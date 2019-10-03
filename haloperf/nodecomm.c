#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "topoinfo.h"

/*D nodecomm - Measure message-passing performance between multicore nodes

Synopsis:
    nodecomm [-v] [-d] [-o] [-o2] [-pm] [-pname name]

+ -v - Turn on verbose mode
. -d - Turn on debug mode
. -o - Test communication/computation overlap with a work array 10 x the
       message size
. -o2 - Test communication/computation overlap with a work array that is
        roughly the square of the message size
. -pm - Print the mapping of processes to nodes and cores.  By default,
        written to nodecomm-mapping.txt .  See -pname to change the output
        file
- -pname name - If -pm is selected, write the mapping to file 'name'.

D*/
#define MAX_DIM 10
#define MAX_SMP_SIZE 64
#define MAX_WORK_SIZE 10000000
/* Getting additional detail about the communication / computation overlap
   requires inserting timing calls within the loop that is used to average
   over multiple runs.  Since this can distort the times, this data is not
   collected by default.  Define GET_OLAP_DETAIL to get this information. */
#ifdef GET_OLAP_DETAIL
#define OLAP_DETAIL(a) a
#else
#define OLAP_DETAIL(a)
#endif

typedef struct {
    int verbose;   /* Provide more information about the operation of code */
    int debug;     /* Use an artificial mesh topology for debugging */
    int overlap;   /* Measure communication overlapped with computation */
    int printMap;  /* Print the mapping of processes to the topology */
    char *mapName; /* Name for the mapping file */
} options_t;

static int verbose = 0;

/* Forward prototypes for functions in this file */
int getOptions(int argc, char *argv[], options_t *options);
int checkConsistentTopo(topoinfo_t *topoinfo);
int findProcessesMatchingMask(int ncoords, const int mycoords[],
			      const int matchcoords[], const int mask[],
			      int ranks[], int nranks);
int findPartnerInCoord( int ncoords, int rootCoords[], int rootRank,
			int myCoords[],
			int idx, int nodeidx, int dist,
			int partnerCoords[], int *partnerRank);
int getNodePingPong(int ppn, const int masterranks[], const int partnerranks[],
		    int minsize, int maxsize);
int getNodePingPongOverlap(int ppn, const int masterranks[],
			   const int partnerranks[],
			   int minsize, int maxsize, int olaptype);
void printiVec(FILE *fp, const char *str, int n, const int vec[]);
int printMapping(FILE *fp, int n, const int mycoords[], MPI_Comm comm);
void ndummy(int, double *);

int main(int argc, char *argv[])
{
    topoinfo_t *topoinfo;
    options_t   options;
    int         required, provided, wrank, wsize, i, level;
    int         mycoords[MAX_DIM], maxcoords[MAX_DIM], ncoords, nodeidx;
    int         mastercoords[MAX_DIM], partnercoords[MAX_DIM];
    int         mask[MAX_DIM];
    int         activemasterranks[MAX_SMP_SIZE], nmasterranks;
    int         activepartnerranks[MAX_SMP_SIZE], npartnerranks;
    int         masterleader, partnerleader;
    int         samenode, prank, ismaster, ispartner, ppn;

    required = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, required, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    /* Get options controling tests */
    if (getOptions(argc, argv, &options)) {
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* Saving verbose in a global variable makes it easier to add verbose
       output. */
    verbose = options.verbose;

    /* Get topology information about job */
    topoDebug(options.debug);
    topoInit(options.verbose, &topoinfo);

    /* Check assumptions about topology (e.g., same number of processes/node) */
    if (checkConsistentTopo(topoinfo)) {
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Determine the communication partners for the tests, based on the
       topology information */
    /*    Basic tests: between two elements (e.g., nodes) with the same
          number, n, of processes.  Do for each level in the topology
	  hierarchy.
          Run with 1..n communication partners.
	  If the topology at a level is a mesh of ndim, use nodes in each
	  coordinate direction.
    */
    /* First, get a tuple describing this process in the topology */
    topoToArray(topoinfo, mycoords, maxcoords, &ncoords, &nodeidx, MAX_DIM);
    if (options.printMap) {
	FILE *fp=0;
	if (wrank == 0) {
	    const char *fname = (const char *)options.mapName;
	    if (!fname) fname = "nodecomm-mapping.txt";
	    fp = fopen(fname, "w");
	    if (!fp) {
		fprintf(stderr, "Could not open %s\n", fname);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
	printMapping(fp, ncoords, mycoords, MPI_COMM_WORLD);
    }
    if (options.verbose > 0) {
	char str[256];
	topoToStr(topoinfo, 1, str, sizeof(str));

	printf("[%d]:%s\n", wrank, str);
	if (options.verbose > 1) {
	    printf("[%d] %d %d\n", wrank, ncoords, nodeidx);

	    for (i=0; i<nodeidx; i++) {
		printf("Node:[%d:%d] %d in %d\n", wrank, i, mycoords[i], maxcoords[i]);
	    }
	    for (i=nodeidx; i<ncoords; i++) {
		printf("Interconnect:[%d:%d] %d in %d\n", wrank, i, mycoords[i], maxcoords[i]);
	    }
	}
    }
    /* First set of tests.  Pingpong within a node.
       One node: The node containing rank 0
       All nodes: (not yet)
    */
    /* Set the master process and tell everyone who it is */
    if (wrank == 0) {
	for (i=0; i<ncoords; i++) mastercoords[i] = mycoords[i];
    }
    MPI_Bcast(mastercoords, ncoords, MPI_INT, 0, MPI_COMM_WORLD);

    /* For each level in the topology, perform the communication tests */
    for (level=0; level<ncoords; level++) {
	MPI_Barrier(MPI_COMM_WORLD);
	/* Pick the partner node.  For level < nodeidx, its the same node.
	   Else, the closest node such that all coords >= nodeidx
	   are the same except for the level'th.

	   We need to skip the cases where there the level'th coordinate
	   is the same for all processes
	*/
	if (level < nodeidx) {
	    samenode = 1;
	    for (i=nodeidx; i<ncoords; i++) {
		partnercoords[i] = mastercoords[i];
		samenode = samenode && (mycoords[i] == mastercoords[i]);
		}
	}
	else {
	    int minmax[2];
	    /* Check that the level'th coordinate is not the same
	       on all processes */
	    minmax[0] = mycoords[level];
	    minmax[1] = -minmax[0];
	    MPI_Allreduce(MPI_IN_PLACE, minmax, 2, MPI_INT, MPI_MAX,
			  MPI_COMM_WORLD);
	    if (minmax[0] == -minmax[1]) {
		if (options.verbose && wrank == 0)
		    printf("Skipping level %d because all coords the same\n",
			   level);
		continue;
	    }
	    samenode = 1;
	    findPartnerInCoord(ncoords, mastercoords, 0, mycoords, level,
			       nodeidx, 1, partnercoords, &prank);
	    for (i=nodeidx; i<ncoords; i++) {
		samenode = samenode && (mycoords[i] == partnercoords[i]);
		}
	}

	/* Run tests with all participating processes */
	for (i=0; i<ncoords; i++)
	    mask[i] = i > level || i >= nodeidx;
	ismaster = 0;
	nmasterranks = findProcessesMatchingMask(ncoords, mycoords,
					   mastercoords, mask,
					   activemasterranks, MAX_SMP_SIZE);
	ispartner = 0;
	npartnerranks = findProcessesMatchingMask(ncoords, mycoords,
					   partnercoords, mask,
					   activepartnerranks, MAX_SMP_SIZE);

	/* At this point, if n{master/partner}ranks > 0,
	   active{master/partner}ranks contains the
	   ranks in MPI_COMM_WORLD of processes that can participate.
	*/
	if (options.verbose) {
	    if (wrank == 0) printf("For level %d\n", level);
	    if (nmasterranks > 0 && wrank == activemasterranks[0]) {
		printiVec(stdout, "activemaster(c)", nmasterranks,
			  activemasterranks);
	    }
	    if (npartnerranks > 0 && wrank == activepartnerranks[0]) {
		printiVec(stdout, "activepartner(c)", npartnerranks,
			  activepartnerranks);
	    }
	}
	/* Set the number of participating processes (ppn) */
	ppn = nmasterranks;
	MPI_Allreduce(MPI_IN_PLACE, &ppn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	/* Pick a master and broadcast the ranks of the processes that
	   will communicate */
	masterleader = 0;
	if (nmasterranks > 0) masterleader = activemasterranks[0];
	MPI_Allreduce(MPI_IN_PLACE, &masterleader, 1, MPI_INT, MPI_MAX,
		      MPI_COMM_WORLD);
	MPI_Bcast(activemasterranks, ppn, MPI_INT, masterleader,
		  MPI_COMM_WORLD);

	partnerleader = 0;
	if (npartnerranks > 0) partnerleader = activepartnerranks[0];
	MPI_Allreduce(MPI_IN_PLACE, &partnerleader, 1, MPI_INT, MPI_MAX,
		      MPI_COMM_WORLD);
	MPI_Bcast(activepartnerranks, ppn, MPI_INT, partnerleader,
		  MPI_COMM_WORLD);

	if (level >= nodeidx) {
	    if (wrank == 0) {
		printf("Pingpong for\n");
		printiVec(stdout,"master:  ", ppn, activemasterranks);
		printiVec(stdout,"partner: ", ppn, activepartnerranks);
	    }
	    if (!options.overlap) {
		getNodePingPong(ppn, activemasterranks, activepartnerranks,
				1, 16*4096);
	    }
	    else {
		getNodePingPongOverlap(ppn, activemasterranks,
				       activepartnerranks, 1, 16*4096,
				       options.overlap);
	    }

	}
    }

    topoFinalize(&topoinfo);
    MPI_Finalize();
    return 0;
}

int getOptions(int argc, char *argv[], options_t *options)
{
    int i;
    /* Set defaults */
    options->verbose  = 0;
    options->debug    = 0;
    options->printMap = 0;
    options->mapName  = 0;
    options->overlap  = 0;

    for (i=1; i<argc; i++) {
	//printf("Processing arg %s\n", argv[i]);
	if      (strcmp(argv[i],"-v") == 0)  options->verbose = 1;
	else if (strcmp(argv[i],"-d") == 0)  options->debug = 1;
	else if (strcmp(argv[i],"-o") == 0)  {
	    options->overlap = 1;
	}
	else if (strcmp(argv[i],"-o2") == 0)  {
	    options->overlap = 2;
	}
	else if (strcmp(argv[i],"-pm") == 0) options->printMap = 1;
	else if (strcmp(argv[i],"-pname") == 0) {
	    i++;
	    if (i < argc)
		options->mapName = strdup(argv[i]);
	    else {
		fprintf(stderr, "-pname missing value\n");
		fflush(stderr);
		return 1;
	    }
	}
	else {
	    fprintf(stderr, "Unrecognized option %s\n", argv[i]);
	    fflush(stderr);
	    return 1;
	}
    }
    return 0;
}

/* This routine performs a sanity check on the topology */
/* FIXME: This should be provided by the topo package */
int checkConsistentTopo(topoinfo_t *topoinfo)
{
    int          tmpval[2];
    int          wrank, i;
    topoentry_t *e;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* All processes have the same number of levels */
    tmpval[0] = topoinfo->numLevels;
    tmpval[1] = -tmpval[0];
    MPI_Allreduce(MPI_IN_PLACE, tmpval, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (tmpval[0] != -tmpval[1]) {
	if (wrank == 0) {
	    fprintf(stderr, "Inconsistent number of levels [%d,%d]\n",
		    -tmpval[1], tmpval[0]);
	}
	return 1;
    }

    /* All processes have the same number of dimensions at each level */
    e = topoinfo->info;
    for (i=0; i<topoinfo->numLevels; i++) {
	tmpval[0] = e->dim;
	tmpval[1] = -tmpval[0];
	MPI_Allreduce(MPI_IN_PLACE, tmpval, 2, MPI_INT, MPI_MAX,
		      MPI_COMM_WORLD);
	if (tmpval[0] != -tmpval[1]) {
	    if (wrank == 0) {
		fprintf(stderr,
		    "Inconsistent number of siblings [%d,%d] at level %d\n",
		    -tmpval[1], tmpval[0], i);
	    }
	    return 1;
	}
	e = e->next;
    }

    /* All is well */
    return 0;
}

/*
 * Find a partner at a specific distance in a specific direction, given
 * the array coordinates.  This is a collective routine that involves
 * communication between processes.
 */
int findPartnerInCoord( int ncoords, int rootCoords[], int rootRank,
			int myCoords[],
			int idx, int nodeidx, int dist,
			int partnerCoords[], int *partnerRank)
{
    int i, mydist, wrank;
    int mindist[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    if (idx >= ncoords) {
	if (wrank == 0)
	    fprintf(stderr,"idx=%d too large (<=%d)\n", idx, ncoords);
	return -1;
    }
    /* For each process, compute a distance from the rootCoords by
       computing a weighted 1-norm distance, with the weights 100x
       in all dimensions other than the idx (and the within-node coords,
       which are ignored) */
    mydist = 0;
    for (i=nodeidx; i<ncoords; i++) {
	int idist = myCoords[i] - rootCoords[i];
	if (idist < 0) idist = -idist;
	if (i != idx) idist *= 100;
	mydist += idist;
    }

    /* Find the rank of a process at the nearest distance */
    /* First, if this process is on the same node (defined by matching
       all coords from nodeidx on, set the distance large so this will
       never be the closest match to the distance */
    if (mydist == 0) mydist = 1000;
    else             mydist = mydist - dist;
    if (mydist < 0) mydist = - mydist;
    if (verbose > 2) printf("[%d] mydist = %d\n", wrank, mydist);
    mindist[0] = mydist;
    mindist[1] = wrank;
    if (verbose > 2) printf("[%d] input mindist (%d,%d)\n", wrank, mindist[0], mindist[1]);
    MPI_Allreduce(MPI_IN_PLACE, mindist, 1, MPI_2INT, MPI_MINLOC,
		  MPI_COMM_WORLD);
    if (verbose > 2) printf("[%d] output mindist (%d,%d)\n", wrank, mindist[0], mindist[1]);

    if (mindist[1] == wrank) {
	for (i=0; i<ncoords; i++) partnerCoords[i] = myCoords[i];
    }
    MPI_Bcast(partnerCoords, ncoords, MPI_INT, mindist[1], MPI_COMM_WORLD);
    *partnerRank = mindist[1];

    return 0;
}

/*
 * Given a set of coordinates and a mask, find all processes such that
 * the coordinates match those in matchcoords for the coords where the
 * mask is 1.  Return their MPI ranks in ranks, an array of size
 * nranks on input.  Return the number of entries found, or 0 if this
 * process does not match the mask.  Return a negative value on error.
 *
 * Since (normally) the master and partner are disjoint, this should be
 * modified to find both sets of processes.  In the comm_split, use a
 * different color for the master and partner.  Can also return which
 * of the two groups the call process belongs (and could generalize to
 * more than two groups)
 */
int findProcessesMatchingMask(int ncoords, const int mycoords[],
			      const int matchcoords[], const int mask[],
			      int ranks[], int nranks)
{
    int i, inmask, wrank, nsize=0, *allranks;
    MPI_Comm newcomm;
    MPI_Group gw, gn;
    inmask = 1;
    for (i=0; i<ncoords; i++) {
	if (mask[i] && mycoords[i] != matchcoords[i]) {
	    inmask = 0; break;
	}
    }
    /* If MPI_Comm_split is efficient and scalable, its the easiest way
       to find the group of processes that have matched the mask
    */
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_split(MPI_COMM_WORLD, inmask, wrank, &newcomm);
    if (inmask) {
	MPI_Comm_size(newcomm, &nsize);
	if (nsize > nranks) return -nsize;
	/* Use Group translate ranks to get the corresponding ranks in
	   MPI_COMM_WORLD for these processes */
	allranks = (int *)malloc(nsize*sizeof(int));
	if (!allranks) return -1;
	for (i=0; i<nsize; i++) allranks[i] = i;
	MPI_Comm_group(MPI_COMM_WORLD, &gw);
	MPI_Comm_group(newcomm, &gn);
	MPI_Group_translate_ranks(gn, nsize, allranks, gw, ranks);
	MPI_Group_free(&gn);
	MPI_Group_free(&gw);
	free(allranks);
    }
    /* This might be a useful communicator; consider an option to return
       it instead */
    MPI_Comm_free(&newcomm);
    return nsize;
}

/*
 * A Basic communication test between all processes in the two arrays
 * of MPI ranks.
 */
int getNodePingPong(int ppn, const int masterranks[], const int partnerranks[],
		    int minsize, int maxsize)
{
    int wrank, nr, master, partner, tag, ntest, i, len, k;
    int ismaster, ispartner;
    int *sbuf, *rbuf;
    double t;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (verbose) {
	/* May want to limit this to "very verbose" */
	printiVec(stdout, "master:  ", ppn, masterranks);
	printiVec(stdout, "partner: ", ppn, partnerranks);
    }
    /* First, determine at what phase, if any, this process is either
       a master or partner */
    ismaster  = ppn;
    ispartner = ppn;
    for (nr=0; nr<ppn; nr++) {
	if (masterranks[nr] == partnerranks[nr]) {
	    fprintf(stderr,"Improper rank arrays at %d\n", nr);
	    MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (masterranks[nr] == wrank) ismaster = nr;
	if (partnerranks[nr] == wrank) ispartner = nr;
    }
    if (verbose) {
	printf("[%d] ismaster=%d, ispartner=%d\n", wrank, ismaster, ispartner);
    }
    sbuf = (int *)malloc(maxsize * sizeof(int));
    rbuf = (int *)malloc(maxsize * sizeof(int));
    for (i=0; i<maxsize; i++) {
	sbuf[i] = i;
	rbuf[i] = -i;
    }

    /* Run the tests over the message sizes */
    for (len=minsize; len<=maxsize; len *= 2) {
	if (wrank == 0) printf("%d\t", len);
	/* FIXME: Need a better way to set ntest */
	ntest = 10000;
	if (len > 10000) ntest = 100;
	for (nr=0; nr<ppn; nr++) {
	    t = -1.0;
	    tag = 1 + nr * ntest;
	    MPI_Barrier(MPI_COMM_WORLD);
	    if (nr >= ismaster) {
		partner = partnerranks[ismaster];
		/*printf("Sending from %d to %d\n", wrank, partner);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    MPI_Send(sbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD);
		    MPI_Recv(rbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD,
			     MPI_STATUS_IGNORE);
		}
		t = (MPI_Wtime() - t)/ntest;
	    }
	    if (nr >= ispartner) {
		master = masterranks[ispartner];
		/*printf("Receiving from %d to %d\n", wrank, master);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    MPI_Recv(rbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD,
			     MPI_STATUS_IGNORE);
		    MPI_Send(sbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD);
		}
		t = (MPI_Wtime() - t)/ntest;
	    }
	    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX,
			  MPI_COMM_WORLD);
	    if (wrank == 0) printf("%.2e\t", t);
	}
	if (len == 0) len = 1;
	if (wrank == 0) printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(rbuf);
    free(sbuf);
    return 0;
}

/*
 * A Basic communication test between all processes in the two arrays
 * of MPI ranks.  This one also measures the overlap of communication
 * over computation
 */
int getNodePingPongOverlap(int ppn, const int masterranks[],
			   const int partnerranks[],
			   int minsize, int maxsize, int olaptype)
{
    int wrank, nr, master, partner, tag, ntest, i, len, k;
    int  workarraysize, wsize;
    double *workarray = 0;
    int ismaster, ispartner;
    int *sbuf, *rbuf;
    double t, tcomm, twork, tocomm; /* comminit, work, commwait */
    MPI_Request req[2];
#ifdef GET_OLAP_DETAIL
    double t1, t2, tolap[3];
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    /*printiVec(stdout, "master:  ", ppn, masterranks);*/
    /*printiVec(stdout, "partner: ", ppn, partnerranks);*/
    /* First, determine at what phase, if any, this process is either
       a master or partner */
    ismaster  = ppn;
    ispartner = ppn;
    for (nr=0; nr<ppn; nr++) {
	if (masterranks[nr] == partnerranks[nr]) {
	    fprintf(stderr,"Improper rank arrays at %d\n", nr);
	    MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (masterranks[nr] == wrank) ismaster = nr;
	if (partnerranks[nr] == wrank) ispartner = nr;
    }
    /*printf("[%d] ismaster=%d, ispartner=%d\n", wrank, ismaster, ispartner);*/
    sbuf = (int *)malloc(maxsize * sizeof(int));
    rbuf = (int *)malloc(maxsize * sizeof(int));
    if (!sbuf || !rbuf) {
	fprintf(stderr, "Unable to allocate sbuf and rbuf\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i=0; i<maxsize; i++) {
	sbuf[i] = i;
	rbuf[i] = -i;
    }
    if (olaptype == 1) {
	workarraysize = 10 * maxsize;
    }
    else {
	if (maxsize < 30000 && maxsize * maxsize < MAX_WORK_SIZE) {
	    workarraysize = maxsize * maxsize;
	}
	else {
	    workarraysize = MAX_WORK_SIZE;
	}
    }
    workarray = (double *)malloc(workarraysize * sizeof(double));
    if (!workarray) {
	fprintf(stderr, "Unable to allocate workarray of size %ld\n",
		(long)workarraysize);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i=0; i<workarraysize; i++)
	workarray[i] = 0;

    /* Run the tests for a range of sizes */
    for (len=minsize; len<=maxsize; len *= 2) {
	/* FIXME: Need a better way to set ntest */
	ntest = 10000;
	if (len > 10000) ntest = 100;
	if (olaptype == 1) {
	    wsize = 10 * len;
	}
	else {
	    if (len < 30000 && len * len < MAX_WORK_SIZE) wsize = len * len;
	    else  wsize = workarraysize;
	}
	if (wrank == 0) printf("%d\t%d\t", len, wsize);
	for (nr=0; nr<ppn; nr++) {
	    tcomm = -1.0;
	    tag = 1 + nr * ntest;
	    /* Communication alone */
	    MPI_Barrier(MPI_COMM_WORLD);
	    if (nr >= ismaster) {
		partner = partnerranks[ismaster];
		/*printf("Sending from %d to %d\n", wrank, partner);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    MPI_Isend(sbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD,
			      &req[0]);
		    MPI_Irecv(rbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD,
			      &req[1]);
		    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
		}
		tcomm = (MPI_Wtime() - t)/ntest;
	    }
	    if (nr >= ispartner) {
		master = masterranks[ispartner];
		/*printf("Receiving from %d to %d\n", wrank, master);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    MPI_Irecv(rbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD,
			      &req[0]);
		    MPI_Isend(sbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD,
			      &req[1]);
		    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
		}
		tcomm = (MPI_Wtime() - t)/ntest;
	    }
	    MPI_Allreduce(MPI_IN_PLACE, &tcomm, 1, MPI_DOUBLE, MPI_MAX,
			  MPI_COMM_WORLD);

	    /* Computation alone for the active processes */
	    twork = -1.0;
	    if (nr >= ismaster || nr >= ispartner) {
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    for (i=0; i<wsize; i++)
			workarray[i] = workarray[i] + 1.0e-6;
		    ndummy(wsize, workarray);
		}
		twork = (MPI_Wtime() - t)/ntest;
	    }
	    MPI_Allreduce(MPI_IN_PLACE, &twork, 1, MPI_DOUBLE, MPI_MAX,
			  MPI_COMM_WORLD);

	    /* Now, do an overlap communication/computation */
	    tocomm = -1.0;
	    OLAP_DETAIL(tolap[0] = 0.0; tolap[1] = 0.0; tolap[2] = 0.0);
	    MPI_Barrier(MPI_COMM_WORLD);
	    if (nr >= ismaster) {
		partner = partnerranks[ismaster];
		/*printf("Sending from %d to %d\n", wrank, partner);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    OLAP_DETAIL(t1 = MPI_Wtime());
		    MPI_Isend(sbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD,
			      &req[0]);
		    MPI_Irecv(rbuf,len,MPI_INT,partner,tag+k,MPI_COMM_WORLD,
			      &req[1]);
		    OLAP_DETAIL(t2 = MPI_Wtime();
				tolap[0] += t2-t1);
		    for (i=0; i<wsize; i++) workarray[i] += 1.0e-6;
		    ndummy(wsize, workarray);
		    OLAP_DETAIL(t1 = MPI_Wtime();
				tolap[1] += t1-t2);
		    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
		    OLAP_DETAIL(tolap[2] += MPI_Wtime() - t1);
		}
		tocomm = (MPI_Wtime() - t)/ntest;
		OLAP_DETAIL(tolap[0] /= ntest;
		tolap[1] /= ntest;
			    tolap[2] /= ntest);
	    }
	    if (nr >= ispartner) {
		master = masterranks[ispartner];
		/*printf("Receiving from %d to %d\n", wrank, master);*/
		t = MPI_Wtime();
		for (k=0; k<ntest; k++) {
		    OLAP_DETAIL(t1 = MPI_Wtime());
		    MPI_Irecv(rbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD,
			      &req[0]);
		    MPI_Isend(sbuf,len,MPI_INT,master,tag+k,MPI_COMM_WORLD,
			      &req[1]);
		    OLAP_DETAIL(t2 = MPI_Wtime();
				tolap[0] += t2-t1);
		    for (i=0; i<wsize; i++) workarray[i] += 1.0e-6;
		    ndummy(wsize, workarray);
		    OLAP_DETAIL(t1 = MPI_Wtime();
				tolap[1] += t1-t2);
		    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
		    OLAP_DETAIL(tolap[2] += MPI_Wtime() - t1);
		}
		tocomm = (MPI_Wtime() - t)/ntest;
		OLAP_DETAIL(tolap[0] /= ntest;
		tolap[1] /= ntest;
			    tolap[2] /= ntest);
	    }

	    MPI_Allreduce(MPI_IN_PLACE, &tocomm, 1, MPI_DOUBLE, MPI_MAX,
			  MPI_COMM_WORLD);
	    OLAP_DETAIL(
	    MPI_Allreduce(MPI_IN_PLACE, tolap, 3, MPI_DOUBLE, MPI_MAX,
			  MPI_COMM_WORLD));

	    if (wrank == 0) {
#ifdef GET_OLAP_DETAIL
		printf("%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t",
		       tcomm, twork, tocomm, tolap[0], tolap[1], tolap[2]);
#else
		printf("%.2e\t%.2e\t%.2e\t", tcomm, twork, tocomm);
#endif
	    }
	}
	if (len == 0) len = 1;
	if (wrank == 0) printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(rbuf);
    free(sbuf);
    return 0;
}


/* Print a vector of integers */
void printiVec(FILE *fp, const char *str, int n, const int vec[])
{
    int i;
    fprintf(fp, "%s", str);
    for (i=0; i<n; i++)
	fprintf(fp,"%s%d", ((i > 0)?", ":""), vec[i]);
    fprintf(fp,"\n");
    fflush(fp);
}

/* Output the mapping of processes to hardware
 *
 * To ensure that the output is correctly ordered, only process 0 in comm
 * writes to the file.
 *
 * The routine is collective and n must have the same value on all
 * processes in comm.
 */
int printMapping(FILE *fp, int n, const int mycoords[], MPI_Comm comm)
{
    int csize, crank, i;
    int *ocoords;

    MPI_Comm_size(comm, &csize);
    MPI_Comm_rank(comm, &crank);
    if (crank == 0) {
	ocoords = (int*)malloc(n * sizeof(int));
	if (!ocoords) return -1;
	fprintf(fp,"%d:", 0);
	printiVec(fp, "", n, mycoords);
	for (i=1; i<csize; i++) {
	    MPI_Recv(ocoords, n, MPI_INT, i, i, comm, MPI_STATUS_IGNORE);
	    fprintf(fp,"%d:", i);
	    printiVec(fp, "", n, ocoords);
	}
	fflush(fp);
	free(ocoords);
    }
    else {
	MPI_Ssend(mycoords, n, MPI_INT, 0, crank, comm);
    }
    return 0;
}
