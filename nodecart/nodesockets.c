/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include "../baseenv.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include "mpi.h"
#include "nodecart.h"

#ifdef HAVE_UTMPX_H
#include <utmpx.h>
#endif

static int verbose = 0;

/* TEMPORRY: Move to baseenv/topo and share */
static int getSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket);

/* TEMPORARY: Move to baseenv/util and share */
static int getIntFromEnv(const char *envname, int *value);

static int getNsocketNcore(int *nsocket, int *ncore);

/*@ MPIX_GetSocketAndCore - Get the socket and Core on which the calling
 process was running

Output Parameters:
+ nsocket - Number of sockets on the node
. socketrank - Rank of the socket on which the process is running
- rankinsocket - Rank of the core in the socket on which the process is
 running

Return Value:
'0' on success and '-1' on failure

Notes:
The values returned reflect the socket and core that was running the calling
process at some time.  The runtime system or operating system may move the
process to another core or socket. The 'rankinsocket' may not be
consecutive; for example, on a Cray XE6, these might be only even values,
since it is common to allocate processes to core-modules, which consist of
a pair of cores.

Environment Variables:
+ TOPO_CORESPERCHIP - If set to an integer value (as a string), then every
  socket (processor chip) is considered to have this many cores.
- TOPO_CHIPSPERNODE - If set to an integer value (as a string), then every
  node is considered to have this many chips (sockets).

  @*/
int MPIX_GetSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket)
{

    getSocketAndCPU(nsocket, socketrank, rankinsocket);
    return 0;
}

/* Various implementations for different systems */

#ifdef HAVE_SCHED_GETCPU
static int coresPerChip = -1;

/* This should work for general Linux and Blue Waters in particular */
static int getSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket)
{
    int    cpunum, chip=0;

    coresPerChip = -1;
    getNsocketNcore(nsocket, &coresPerChip);
    /* On error or failure, leaves values as -1 */

    /* Attempt to get the total number of physical cores */
#ifdef HAVE_SYSCTLBYNAME
    if (coresPerChip <= 0 && *nsocket > 0) {
	int    rc, totalcores = -1;
	size_t nlen=sizeof(int);
	rc = sysctlbyname("hw.physicalcpu", &totalcores, &nlen, NULL, 0);
	/* Also of possible interest (and in CTL_MACHDEP):
	   machdep.cpu.core.out, machdep.cpu.cores_per_package */
	if (rc == 0 && totalcores > 0)
	    coresPerChip = *nsocket/totalcores;
    }
#endif

    /* Determine where this process is running */
    cpunum = sched_getcpu();
    if (verbose) printf("cpu # = %d\n", cpunum);
    /* On some systems, the cores are numbered consecutively and the same
       number of cores per chip.  In this case, we can determine the chip
       by looking at the cpunum */
    if (cpunum >= 0 && coresPerChip > 0) {
	*socketrank   = cpunum / coresPerChip;
	/* Normalize the cpunum as per chip */
	*rankinsocket = cpunum - *socketrank * coresPerChip;
    }
    return 0;
}
#elif defined(NOTUSED)
/* This could use hwloc, for example */
static int getSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket)
{
    return -1;
}
#else
/* No information available */
static int getSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket)
{
    *nsocket      = -1;
    *socketrank   = -1;
    *rankinsocket = -1;
    return -2;
}
#endif

#include <errno.h>
/* Get an integer value from an environment value.  If the environment
   variable is a valid integer, value is set and the routine returns 0.
   If it is not set, value is unchanged and the routine returns 1.
   If there is an error, a message will be printed and the routine returns
   -1
*/
static int getIntFromEnv(const char *envname, int *value)
{
    char       *s;

    s = getenv(envname);
    if (s && *s) {
	char *endptr;
	/* Use strtol since atoi has no error indicator. strtol is awkward,
	   but it has an error indicator */
	long val;
	errno = 0;
	/* A base of zero permits 0nn for octal, 0x for hex, and other for 10 */
	val = strtol(s,&endptr,0);
	if (errno != 0 || endptr == 0 || *endptr != '\0') {
	    if (verbose) {
		fprintf(stderr, "Invalid value for %s: %s\n", envname, s );
	    }
	    return -1;
	}
	else if (val > 0) {
	    *value = (int)val;
	    return 0;
	}
    }
    return 1;
}

/* Routines to get the numbers of sockets and number of cores on each
   socket, assuming all sockets have the same number of cores.
   Prefers to use hwloc, then read /proc/cpuinfo, then use external
   information in environment variables (i.e., have the user read the
   documentation on the system and set values).
*/
static int nsocketVal=-1, ncoreVal=-1, nninit=0;
#ifdef HAVE_HWLOC_H
#include "hwloc.h"
static hwloc_topology_t hwtopology;
/* This routine inspired by hwloc/doc/hwloc-hello.c */
static int getNsocketNcore(int *nsocket, int *ncore)
{
    int depth;

    if (nninit) {
	*nsocket = nsocketVal;
	*ncore   = ncoreVal;
	return 0;
    }
    hwloc_topology_init(&hwtopology);
    hwloc_topology_load(hwtopology);
    nninit = 1;

    depth = hwloc_get_type_depth(hwtopology, HWLOC_OBJ_SOCKET);
    if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
	/* Can't get values - leave them as -1 */
	hwloc_topology_destroy(hwtopology);
	return 1;
    }
    nsocketVal = hwloc_get_nbobjs_by_depth(hwtopology, depth);

    depth = hwloc_get_type_depth(hwtopology, HWLOC_OBJ_CORE);
    if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
	/* Return leaving ncoreVal set to -1 */
	hwloc_topology_destroy(hwtopology);
	return 1;
    }
    ncoreVal = hwloc_get_nbobjs_by_depth(hwtopology, depth);
    hwloc_topology_destroy(hwtopology);

    *nsocket = nsocketVal;
    *ncore   = ncoreVal;
    return 0;
}
#elif defined(HAVE_PROC_CPUINFO) && 0
/* This is harder than it looks.  This code returns the wrong core
   count on AMD Interlagos, because /proc/cpuinfo returns two groups
   of eight cores for the same phyiscal package.  Another system
   returned physical ids of 0 and 2 (no 1), so the number of distinct
   ids, rather than the largest, must be used to count packages.*/
#include <ctype.h>
static int getNsocketNcore(int *nsocket, int *ncore)
{
    FILE *fp;
    char fieldname[128], value[1024], inputline[1024];
    char *p;
    int  rc, ncores, maxncores=-1, chipno, maxchipno=-1;

    if (nninit) {
	*nsocket = nsocketVal;
	*ncore   = ncoreVal;
	return 0;
    }
    fp = fopen("/proc/cpuinfo", "r");
    while(!feof(fp)) {
	fgets(inputline,1024,fp);
	p = strchr(inputline,':');
	if (p) {
	    char *f=fieldname, *p1=inputline;
	    while (*p1 != ':' && p1 < inputline+1023) { *f++ = *p1++; } *f=0;
	    while (isspace(f[-1]) && f > fieldname) *--f = 0;
	    f = value; p1 = p+1;
	    while (isspace(*p1) && p1 < inputline+1023) p1++;
	    while (*p1 && *p1 != '\n' && p1 < inputline+1023) { *f++ = *p1++; }
	    *f=0;

	    if (strcmp(fieldname, "physical id") == 0) {
		chipno = atoi(value);
		if (chipno > maxchipno) maxchipno = chipno;
	    }
	    else if (strcmp(fieldname, "cpu cores") == 0) {
		ncores = atoi(value);
		if (ncores > maxncores) maxncores = ncores;
	    }
	}
    }
    fclose(fp);
    nninit = 1;
    if (maxchipno != -1) {
	nsocketVal = maxchipno + 1;
	*nsocket   = nsocketVal;
    }
    if (maxncores > 0) {
	ncoreVal  = maxncores;
	*ncore    = ncoreVal;
    }
    return 0;
}
#else
/* Default - use environment variables if available */
static int getNsocketNcore(int *nsocket, int *ncore)
{
    if (nninit) {
	*nsocket = nsocketVal;
	*ncore   = ncoreVal;
	return 0;
    }
    /* Get external configuration */
    if (getIntFromEnv("TOPO_SOCKETSPERNODE", &nsocketVal) == 0) {
	*nsocket = nsocketVal;
    }
    else if (getIntFromEnv("TOPO_CHIPSPERNODE", &nsocketVal) == 0) {
	*nsocket = nsocketVal;
    }
    if (getIntFromEnv("TOPO_CORESPERCHIP", &ncoreVal) == 0) {
	*ncore = ncoreVal;
    }
    nninit = 1;
    return 0;
}
#endif
