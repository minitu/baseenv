/* -*- Mode: C; c-basic-offset:4 ; -*- */
/* At least on some systems, _GNU_SOURCE is required to get the prototype for
   sched_getcpu and to get the correct implementation for CPU_SET.
   Further, it needs to be set early; it is too late to set it when
   <sched.h> is loaded. */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "mpi.h"
#include "topoconf.h"
#include "topoinfo.h"
#include "topoimpl.h"

/*static int topoiGetNodeInfo(int, topoinfo_t *);
  static int topoiGetNodeTopoInfo(topoinfo_t *); */
static int getIntFromEnv(const char *envname, int *value);
static void topoiInitFns(void);
static int topoiSetupNodeInfo(topoinfo_t *topo);

static int coresPerChip = -1;
static int chipsPerNode = -1;

static int verbose = 0;

/* Which set of available topology routines to use */
#define MAX_TOPO_FNS 5
#define MAX_NODE_FNS 5
static int nodeFnsNum = -1,
    nodeFnsIdx = -1,
    netFnsNum = -1,
    netFnsIdx = -1;
static topoiNodeFns_t topoNodeFns[MAX_NODE_FNS];
static topoiNetFns_t  topoNetFns[MAX_TOPO_FNS];

/*
 * Provide basic topology information about the calling thread.  See
 * README for information.
 */
/*@ topoInit - Initialize the hardware topology package

Input Parameter:
. isVerbose - If true, provide verbose output to stdout about the operation
  of the routines.

Output Parameter:
. topo_p - On successful exit, points to a 'topoinfo_t' structure.  This
  lists the topology of the system in which the calling process is executing
  starting from processor cores and working out.  See notes below.

Collective:
  This routine is collective over all processes in the job.  It is not
  necessarily synchronizing, but some systems may require collective
  operations in order to return the requested data.

Return value:
0 on success, non-zero on failure.

Environment Variables:
+ TOPO_CORESPERCHIP - If set to an integer value (as a string), then every
  socket (processor chip) is considered to have this many cores.
- TOPO_CHIPSPERNODE - If set to an integer value (as a string), then every
  node is considered to have this many chips (sockets).

Notes:
The topo package is designed for massively parallel systems with regular
networks, such as a torus or mesh.  If the hwloc package is available,
this package can use hwloc to provide some information about the node
topology.  The 'topoinfo_t' structure provides a list of entries, starting
with the finest grain (typically a core, but could be a chip or node if
more detailed information is not available).

Systems Supported:

Note - only IBM, Cray, and HWLOC currently implemented.

+ IBM Blue Gene/Q - Provides interconnect mesh coordinates
. Cray XE6 and XK7 - Provides interconnent mesh coordinates and core number;
  with 'TOPO_CORESPERCHIP' correctly set, also provides socket number
. Miscellaneous MPI programs that have 'MPI_Comm_split_type'.  This
  routine may be used to determine a topology based on processes that can
  share memory.  This is represented as a 'TOPO_NODE'.
- Miscellaneous Unix systems - With hwloc, provides node information.
  if 'sched_getcpu' is available, this is used to find the core number
  (use the environment variable 'TOPO_CORESPERCHIP' to identify to which
  chip in a multi-chip node a core belongs).
  Without hwloc, the code uses the hostname to identify processes that
  share the same hostname, and identifies these as belonging to the same
  'TOPO_NODE'.

Thread safety:
  If called by multiple threads, they must provide the same value for
  'isVerbose'.
@*/
int topoInit(int isVerbose, topoinfo_t **topo_p)
{
  int        err;
  topoinfo_t *topo;
  static volatile int isInitialized=0;

  /* This uses the "double init".  The thread safety of this is suspect,
     but should work most of the time.  Making this fully thread safe
     requires either fixing on a particular thread interface or using
     non-portable memory fence operations. */
  if (isInitialized == 0) {
      verbose = isVerbose;

      /* Get external configuration */
      if (getIntFromEnv("TOPO_CORESPERCHIP", &coresPerChip) < 0)
	  return -1;
      if (getIntFromEnv("TOPO_CHIPSPERNODE", &chipsPerNode) < 0)
	  return -1;

      topoiInitFns();
      isInitialized = 1;
  }

  *topo_p = (topoinfo_t*)malloc(sizeof(topoinfo_t));
  if (!*topo_p) {
      if (verbose) {
	  fprintf(stderr, "Unable to allocate topoinfo_t\n");
      }
      return -1;
  }
  topo                 = *topo_p;
  topo->numLevels      = 0;
  topo->ranks          = 0;
  topo->nodenum        = -1;   /* Unknown node number */
  topo->nnodes         = -1;   /* Unknown number of nodes */
  topo->nodeComm       = MPI_COMM_NULL;
  topo->nodeLeaderComm = MPI_COMM_NULL;
  topo->info           = 0;
  topo->infoTail       = 0;

  /* List from the finest grain out */
  err = topoNodeFns[nodeFnsIdx].getNodeInfo(0, topo);
  if (err != -1)
      err = topoNetFns[netFnsIdx].getNodeTopoInfo(topo);

  /* Setup node information if necessary */
  if (topo->nnodes == -1)
      topoiSetupNodeInfo(topo);

  return err;
}

/*@ topoFinalize - Finalize the topology package

Input/Output Parameter:
. topo - Address of the pointer returned by 'topoInit'.  Set to NULL on exit.

Return value:
0 on success, non-zero on failure.
@*/
int topoFinalize(topoinfo_t **topo)
{
    int err = 0;
    topoentry_t *e, *enext;

    if (*topo) {
	/* Free all elements */
	e = (*topo)->info;
	while (e) {
	    enext = e->next;
	    free(e);
	    e = enext;
	}
	if ((*topo)->nodeComm != MPI_COMM_NULL)
	    MPI_Comm_free(&(*topo)->nodeComm);
	if ((*topo)->nodeLeaderComm != MPI_COMM_NULL)
	    MPI_Comm_free(&(*topo)->nodeLeaderComm);
	if ((*topo)->ranks)
	    free((*topo)->ranks);
	free(*topo);
	*topo = 0;
    }
    return err;
}

/*@
  topoAvailMethods - Return the available methods for determining the topology

Input Parameter:
. maxnum - Size of the arrays 'nodeMethod' and 'netMethod'

Output Parameters:
+ nnode - Number of node methods defined
. nnet - Number of network methods defined
. nodeMethod - Array containing an entry for each node method
- netMethod - Array containing an entry for each network method

Notes:
 This routine returns the available methods to determine the topology
 information within the node and between nodes (the internode network).
 Values are taken from the enum 'topoMethod_t'.  The methods may be used
 in 'topoSetMethods' to specify how the topology is determined when 'topoInit'
 is called.

  @*/
int topoAvailMethods(int maxnum, int *nnode, int *nnet,
		     topoMethod_t nodeMethod[], topoMethod_t netMethod[])
{
    int i, n;

    topoiInitFns();

    *nnode = nodeFnsNum+1;
    *nnet  = netFnsNum+1;
    n = (nodeFnsNum+1 < maxnum) ? nodeFnsNum+1 : maxnum;
    for (i=0; i<n; i++) nodeMethod[i] = topoNodeFns[i].id;
    n = (netFnsNum+1 < maxnum) ? netFnsNum+1 : maxnum;
    for (i=0; i<n; i++) netMethod[i] = topoNetFns[i].id;

    return 0;
}

/*@
  topoSetMethods - Set the methods used to determine the topology

Input Parameters:
+ nodeMethod - The method used to determine the topology of the node (e.g.,
 the number of cores and sockets)
- netMethod - The method used to determine the topology connecting the
 nodes (e.g., a mesh for a Cray XE or IBM BlueGene).

Notes:
 Valid methods can be discovered by calling 'topoAvailMethods'.
  @*/
int topoSetMethods(topoMethod_t nodeMethod, topoMethod_t netMethod)
{
    int err = 0;
    for (nodeFnsIdx=0; nodeFnsIdx <= nodeFnsNum; nodeFnsIdx++) {
	if (topoNodeFns[nodeFnsIdx].id == nodeMethod) break;
    }
    for (netFnsIdx=0; netFnsIdx <= netFnsNum; netFnsIdx++) {
	if (topoNetFns[netFnsIdx].id == netMethod) break;
    }
    if (nodeFnsIdx > nodeFnsNum) {
	nodeFnsIdx = 0;
	err        = 1;
    }
    if (netFnsIdx > netFnsNum) {
	netFnsIdx = 0;
	err       = 1;
    }
    return err;
}

/*@ topoGetMethodsDesc - Return the description string for the currently
  selected topology methods

Output Parameters:
+ nodestr - Describes the method used to determine the node topology
- netstr  - Describes the method used to determine the interconnect topology
  @*/
int topoGetMethodsDesc(const char **nodestr, const char **netstr)
{
    *nodestr = topoNodeFns[nodeFnsIdx].descString;
    *netstr  = topoNetFns[netFnsIdx].descString;
    return 0;
}


/*@
  topoDebug - Select a virtual topology to return

Input Parameter:
.  kind - Kind of topology to return.

Return value:
 0 on success, non-zero on failure.

Notes:
  Normally, the topo routines return information about the real topology
  of the underlying hardware.  However, this can make debugging programs
  that want to make use of these routines difficult, since the tests would
  have to be run on the target hardware system.  This routine allows the
  use of dummy, virtual topologies.  Valid values for 'kind' include\:
.vb
    kind = 1: Use a multiple of 8 processes, provide a Cray XE6-like
              topology
.ve
  @*/
int topoDebug(int kind)
{
    topoiSetDummyTopo(kind);
    /* Set the methods to the debug methods */
    topoSetMethods(TOPO_DEBUG, TOPO_DEBUG);
    return 0;
}

/*@ topoTypeStr - Return a string name for a topology type

Input Parameter:
. t - A 'topoType_t' value

Return Value:
 A string name for the topology type.  On error, returns 'invalid topo type!'
 The known names include 'torus', 'mesh', 'socket', and 'core'.

Thread Safety:
 This routine is thread safe.
  @*/
const char *topoTypeStr(topoType_t t)
{
  const char *p = 0;
  switch (t) {
  case TOPO_UNKNOWN: p = "unknown"; break;
  case TOPO_NODELIST:p = "nodelist"; break;
  case TOPO_TORUS:   p = "torus"; break;
  case TOPO_MESH:    p = "mesh" ; break;
  case TOPO_BUS:     p = "bus";   break;
  case TOPO_TREE:    p = "tree"; break;
  case TOPO_NODE:    p = "node"; break;
  case TOPO_SOCKET:  p = "socket"; break;
  case TOPO_CORE:    p = "core"; break;
  default:
    p = "invalid topo type!";
  }
  return p;
}

/*@ topoPrint - Print basic information about the topology

Input Parameters:
+ fp - FILE pointer for output
. leader - string to prepend each line of output
- topo - pointer to topology info

Return value:
0 on success, non-zero on failure.

Thread Safety:
 This routine is not thread safe.
@*/
int topoPrint(FILE *fp, const char *leader, topoinfo_t *topo)
{
  topoentry_t *e = topo->info;
  while (e) {
    printf("%s%s: %d:", leader, topoTypeStr(e->topoType), e->dim);
    if (e->coords.coords[0] >= 0) {
      int i;
      for (i=0; i<e->dim; i++)
	printf("%d%s", e->coords.coords[i], (i+1<e->dim) ? "," : "");
    }
    printf("\n");
    e = e->next;
  }
  return 0;
}

/*@ topoToStr - Create a one-line description of the topology, suitable for printing

Input Parameters:
+ topo - Pointer to topology information
. incname - if true, include the string name of each element of the topology
- len - size of 'str'

Output Parameters:
. str - Pointer to a 'char' array that will contain the description of the
topology.

Return value:
0 on success, non-zero on failure.

@*/
int topoToStr(topoinfo_t *topo, int incname, char *str, int len)
{
    topoentry_t *e = topo->info;
    int lenleft = len;
    int addcomma = 0;
    while (e && lenleft > 0) {
	if (incname) {
	    const char *name = topoTypeStr(e->topoType);
	    int i, ilen = strlen(name);
	    if (ilen + 3 < lenleft) {
		if (addcomma) {
		    *str++ = ',';
		    lenleft--;
		}
		for (i=0; i<ilen; i++) {
		    *str++ = *name++;
		    lenleft--;
		}
		*str++ = ':';
		lenleft--;
		addcomma = 0;
	    }
	    else break;
	}
	if (e->coords.coords[0] >= 0) {
	    char digits[10], *d;
	    int i, ilen;
	    for (i=0; i<e->dim; i++) {
		snprintf(digits, 10, "%d", e->coords.coords[i]);
		ilen = strlen(digits);
		if (ilen + 2 < lenleft) {
		    if (addcomma) {
			*str++ = ',';
			lenleft--;
		    }
		    d = digits;
		    while (*d) {
			*str++ = *d++;
			lenleft--;
		    }
		    addcomma = 1;
		}
		else break;
	    }
	}
	e = e->next;
    }
    *str++ = 0;
    return 0;
}

/*@
  topoNodeInfoBasic - Return very basic information on the node

  Input Parameter:
. topo - Pointer to topoinfo that was returned by topoInit.

  Output Parameters:
+ corenum - Points to the core number on the node.  -1 if unknown.
. socketnum - Points to the socket (chip) number on the node.  -1 if unknown.
. corespersocket - Number of cores (defined by the system) per socket (chip).
  -1 if unknown.
- socketspernode - Number of sockets (processor chips) per node.  -1 if unknown.

  Return value:
  Zero on success and non-zero on failure
  @*/
int topoNodeInfoBasic(topoinfo_t *topo, int *corenum, int *socketnum,
		      int *corespersocket, int *socketspernode)
{
    topoentry_t *e = topo->info;

    *corenum        = -1;
    *socketnum      = -1;
    *corespersocket = -1;
    *socketspernode = -1;
    while (e) {
	if (e->topoType == TOPO_CORE) {
	    if (e->coords.coords[0] >= 0) *corenum = e->coords.coords[0];
	    if (e->maxcoords.coords[0] > 0)
		*corespersocket = e->maxcoords.coords[0];
	}
	else if (e->topoType == TOPO_SOCKET) {
	    if (e->coords.coords[0] >= 0) *socketnum = e->coords.coords[0];
	    if (e->maxcoords.coords[0] > 0)
		*socketspernode = e->maxcoords.coords[0];
	}
	e = e->next;
    }
    return 0;
}

/*@
  topoNodeEnumeration - Return a numbering of the nodes and the ranks
  of the processes on the same node

Input Parameter:
. topo - Pointer to topoinfo that was returned by topoInit.

Input/Output Parameter:
. nranks - On input, the size of 'noderanks'.  On output, the number of
 processes on this node (including the calling process)

Output Parameters:
+ numnodes - Number of nodes available
. mynodenum - Number of the node on which this process is running
- noderanks - Ranks, in MPI_COMM_WORLD, of processes on the same node as the
 calling process.  In the value of '*nranks' on input is less than the
 total number of ranks, only the first '*nranks' values are returned, and
 '*nranks' is set to the total number of ranks.

Return value:

Notes:
  @*/
int topoNodeEnumeration(topoinfo_t *topo, int *numnodes, int *mynodenum,
			int *nranks, int noderanks[])
{
#if 0
    topoentry_t *e = topo->info;
#endif
    int          i, nr;

    if (topo->nnodes != -1) {
	*numnodes  = topo->nnodes;
	*mynodenum = topo->nodenum;
	nr = topo->nranks;
	if (nr > *nranks) nr = *nranks;
	*nranks    = topo->nranks;
	for (i=0; i<nr; i++) {
	    noderanks[i] = topo->ranks[i];
	}
	return 0;
    }

    *nranks = -1;
    return -1;
#if 0
    *mynodenum        = -1;

    /* Find the node information */
    while (e) {
	/* Skip elements within a node */
	if (e->topoType == TOPO_SOCKET || e->topoType == TOPO_CORE) {
	    e = e->next;
	    continue;
	}

	if (e->topoType == TOPO_NODE) {
	    /* Keep looking */
	}
	}
	else if (e->topoType == TOPO_MESH || e->topoType == TOPO_TORUS) {
	    /* The processes on the same nodes are the ones with the
	       same coordinates in the interconnect.  At this point, we
	       don't have an address, but we can use MPI_Comm_split
	       to determine the same group */
	    int color, wrank, nsize, i, *ranks;
	    MPI_Comm nodecomm, listcomm;

	    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
	    color = e->coords.coords[0];
	    for (i=1; i<e->dim; i++)
		color = e->coords.coords[i] + e->maxcoords[i-1]*color;
	    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &nodecomm);

	    MPI_Comm_rank(nodecomm, &nrank);
	    MPI_Comm_group(MPI_COMM_WORLD, &gw);
	    MPI_Comm_group(nodecomm, &gn);
	    ranks    = (int *)malloc(nsize * sizeof(int));
	    e->ranks = (int *)malloc(nsize * sizeof(int));
	    for (i=0; i<nranks; i++) ranks[i] = i;
	    MPI_Group_translate_ranks(gn, nsize, ranks, gw, &e->ranks);
	    free(ranks);
	    MPI_Group_free(&gn);
	    MPI_Group_free(&gw);

	    color = (nrank == 0) ? 1 : MPI_UNDEFINED;
	    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &listcomm);
	    MPI_Comm_rank(listcomm, );
	    MPI_Comm_size(listcomm, );
	    MPI_Comm_free(&nodecomm);
	}
	else {
	    /* Type BUS, TREE, or UNKNOWN */
	    /* Did not find a node type yet.  Argh. */
	}
	e = e->next;
    }
    /* Did not find a way to determine the number of nodes */
    *nranks = -1;
#endif
}

/*@
  topoMeshCoords - If the topology contains a mesh or torus, return the coordinates of the calling process in that mesh or torus

Input Parameter:
. topo - Pointer to topology information

Input/Output Parameter:
. ndim - Pointer to dimension of mesh or torus.  On input, the size of the
  arrays 'meshcoords' and 'qtorus'.  On output, the size of the mesh or torus.

Output Parameters:
+ meshcoords - The coordinates of the calling process in the mesh.
- qtorus - 1 if this dimension is periodic, 0 otherwise.

Return Value:
 0 on success (a mesh topology element was found and the data returned).
 Nonzero (currently '1') on failure, and '-1' if the topology does not
 contain a mesh or torus.

Notes:
The current implementation does not indicate whether any particular dimension
is periodic, only whether the entire multidimensional mesh is periodic.
In addition, a subset of a torus is still returned as a torus even though the
min and max coordinates in this dimension are not connected directly together.
 @*/
int topoMeshCoords(topoinfo_t *topo, int *ndim, int meshcoords[], int qtorus[])
{
    topoentry_t *e = topo->info;
    int          i, foundMesh = 0;

    while (e) {
	if (e->topoType == TOPO_MESH) {
	    if (e->coords.coords[0] >= 0) {
		if (*ndim < e->dim) return 1;
		*ndim = e->dim;
		for (i=0; i<e->dim; i++) {
		    meshcoords[i] = e->coords.coords[i];
		    qtorus[i]     = 0;
		}
	    }
	    foundMesh = 1;
	    break;
	}
	else if (e->topoType == TOPO_TORUS) {
	    if (e->coords.coords[0] >= 0) {
		if (*ndim < e->dim) return 1;
		*ndim = e->dim;
		for (i=0; i<e->dim; i++) {
		    meshcoords[i] = e->coords.coords[i];
		    qtorus[i]     = 1;
		}
	    }
	    foundMesh = 1;
	    break;
	}
	e = e->next;
    }
    return foundMesh ? 0 : -1;
}

/*@
  topoMeshContainer - If the topology contains a mesh or torus, return the minimum and maximum dimensions of the n-d section that contains all processes

Input Parameter:
. topo - Pointer to topology information

Input/Output Parameter:
. ndim - Pointer to dimension of mesh or torus.  On input, the size of the
  arrays 'mindim', 'maxdim', and 'qtorus'.  On output, the size of the
  mesh or torus.

Output Parameters:
+ mindim - The minimum coordinates of any process in the mesh or torus.
. maxdim - The maximum coordinates of any process in the mesh or torus.
- qtorus - 1 if this dimension is periodic, 0 otherwise.

Return Value:
 0 on success (a mesh topology element was found and the data returned).
 Nonzero (currently '1') on failure, and '-1' if the topology does not
 contain a mesh or torus.

Notes:
The current implementation does not indicate whether any particular dimension
is periodic, only whether the entire multidimensional mesh is periodic.
In addition, a subset of a torus is still returned as a torus even though the
min and max coordinates in this dimension are not connected directly together.
 @*/
int topoMeshContainer(topoinfo_t *topo, int *ndim,
		      int mindim[], int maxdim[], int qtorus[])
{
    topoentry_t *e = topo->info;
    int          i, foundMesh = 0;

    while (e) {
	if (e->topoType == TOPO_MESH) {
	    if (e->coords.coords[0] >= 0) {
		if (*ndim < e->dim) return 1;
		*ndim = e->dim;
		for (i=0; i<e->dim; i++) {
		    mindim[i] = e->mincoords.coords[i];
		    maxdim[i] = e->maxcoords.coords[i];
		    qtorus[i]     = 0;
		}
	    }
	    foundMesh = 1;
	    break;
	}
	else if (e->topoType == TOPO_TORUS) {
	    if (e->coords.coords[0] >= 0) {
		if (*ndim < e->dim) return 1;
		*ndim = e->dim;
		for (i=0; i<e->dim; i++) {
		    mindim[i] = e->mincoords.coords[i];
		    maxdim[i] = e->maxcoords.coords[i];
		    qtorus[i]     = 1;
		}
	    }
	    foundMesh = 1;
	    break;
	}
	e = e->next;
    }
    return foundMesh ? 0 : -1;
}

/*
 Algorithm:
 1. Compute a unique index number for each node, based on the mesh size.
 2. Use this number to compute a communicator for each node.
    Then split, with only the root of the node communicator partipating on
    each node (use MPI_UNDEFINED for color for the other processes on the
    node).
 3. Find the node with minimum coordinate.  That becomes node 0.
 4. While not done: Find next closest node.  Break ties by taking node
    with the minimum unique index
 5. Broadcast node number to the processes on the same node.

 Notes:
 This algorithm is not scalable, and should be used only for modest
 numbers of nodes.
 */
typedef struct { int coords[5], dist, rank, nodenum; } ringdist_t;
int topoGetRingFromMesh(topoinfo_t *topo, int ndim, const int meshcoords[],
			const int qtorus[], int *nodenum)
{
    MPI_Comm nodecomm, rootcomm;
    int      wrank, uidx, i, j, k, nrank, color;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* 1. Compute the unique index */
    uidx = 0;
    for (i=0; i<ndim; i++)
	uidx = uidx * 24 + meshcoords[i];

    /* 2a. Create the node communicator */
    MPI_Comm_split(MPI_COMM_WORLD, uidx, wrank, &nodecomm);

    /* 2b. Create the communicator of node roots */
    MPI_Comm_rank(nodecomm, &nrank);
    color = (nrank == 0) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &rootcomm);

    if (nrank == 0) {
	int mynodenum = -1;
	int rsize, rrank, minuidx, mindist, rmindist, currank, dist;
	int *allcoords;
	ringdist_t *elm;
	MPI_Comm_size(rootcomm, &rsize);
	MPI_Comm_rank(rootcomm, &rrank);

	/* 3. Find the node with minimum uidx */
	MPI_Allreduce(&uidx, &minuidx, 1, MPI_INT, MPI_MIN, rootcomm);
	if (minuidx == uidx) {
	    mynodenum = 0;
	}

	/* 4. Gather the coordinates, sort them, and assign */
	allcoords = (int *)malloc(ndim*rsize*sizeof(int));
	MPI_Allgather((int*)meshcoords, ndim, MPI_INT, allcoords, ndim, MPI_INT,
		      rootcomm);
	/* 4a: Load up the structure describing the entries */
	elm = (ringdist_t *)malloc(rsize * sizeof(ringdist_t));
	mindist  = 100000;
	rmindist = -1;
	for (i=0; i<rsize; i++) {
	    for (j=0; j<ndim; j++)
		elm[i].coords[j] = allcoords[i*ndim+j];
	    elm[i].rank    = i;
	    elm[i].nodenum = -1;
	    dist = 0;
	    for (j=0; j<ndim; j++)
		dist += elm[i].coords[j];
	    elm[i].dist    = dist;
	    if (dist < mindist) {
		mindist  = dist;
		rmindist = i;
	    }
	}
	elm[rmindist].nodenum = 0;
	currank = rmindist;
	for (i=1; i<rsize; i++) {
	    /* Find the next closest process */
	    mindist = 100000;
	    rmindist = -1;
	    for (j=0; j<rsize; j++) {
		if (elm[j].nodenum >= 0) continue;
		dist = 0;
		for (k=0; k<ndim; k++)
		    dist += elm[j].coords[k] - elm[currank].coords[k];
		if (dist < mindist) {
		    mindist  = dist;
		    rmindist = j;
		}
	    }
	    elm[rmindist].nodenum = i;
	}
	*nodenum = elm[rrank].nodenum;
	free(allcoords);
	free(elm);
    }
    MPI_Bcast(nodenum, 1, MPI_INT, 0, nodecomm);

    if (rootcomm != MPI_COMM_NULL) MPI_Comm_free(&rootcomm);
    MPI_Comm_free(&nodecomm);
    return 0;
}

/*@
  topoToArray - Convert a topoinfo structure into a array of coordinates

Input Parameters:
+ topo - Pointer to topology information
- maxdim - Maximum size of arrays 'mycoords' and 'maxcoords'

Output Parameters:
+ mycoords - a nonnegative integer giving a value for each `level` of the
  topology, starting from the smallest (typically a processor core)
. maxcoords - The maximum value of the coordinate for this level
. nlevels - The number of `levels` of the topology.
- nodeidx - The index of the first level of the topology above that of
  a node.

Notes:
This routine provides a way to convert a topoinfo structure into an array of
coordinates that can be used map calling processes onto the physical
topology.  'mycoords[i]' contains a coordinate for the calling process at the
ith level, with the finest level at i=0.  Levels between 0 and nodeidx-1 are
within a single node.  Level 'nodeidx' (if less than 'nlevels') is the first
level of the interconnect between nodes.  If the interconnnect is a mesh,
there will be a separate level for each dimension of the interconnect.
  @*/
int topoToArray(topoinfo_t *topoinfo, int mycoords[], int maxcoords[],
		int *nlevels, int *nodeidx, int maxdim)
{
  topoentry_t *e = topoinfo->info;
  int              i, level=0, foundNode = 0;
  while (e) {
      if (level >= maxdim) {
	  fprintf(stderr, "Too many levels in topology for arrays\n");
	  return 1;
      }
      if (foundNode == 0 &&
	  e->topoType != TOPO_NODE &&
	  e->topoType != TOPO_SOCKET &&
	  e->topoType != TOPO_CORE) {
	  *nodeidx = level;
	  foundNode = 1;
      }
      for (i=0; i<e->dim && level < maxdim; i++) {
	  mycoords[level]  = e->coords.coords[i];
	  maxcoords[level] = e->maxcoords.coords[i];
	  level++;
      }
      e = e->next;
  }
  if (!foundNode) *nodeidx = level;
  *nlevels = level;
  return 0;
}

/* Service routines */
/* Allocate an entry */
topoentry_t *topoiAllocEntry(topoinfo_t *topo)
{
  topoentry_t *e;
  e = (topoentry_t *)malloc(sizeof(topoentry_t));
  if (!e) {
      if (verbose) {
	  fprintf(stderr, "Unable to allocate a topoentry_t\n");
      }
      return NULL;
  }
  e->coords.coords[0]        = -1;
  e->maxcoords.coords[0]     = -1;
  e->mincoords.coords[0]     = -1;
  e->maxtopocoords.coords[0] = -1;
  e->next                    = 0;
  if (topo) {
    if (topo->infoTail) topo->infoTail->next = e;
    if (!topo->info) topo->info = e;
    topo->infoTail = e;
  }
  return e;
}


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

/* Include the system-specific routines for this architecture */
#ifdef HAVE_CRAY_RCA
#include "topocray.inc"
#endif
#ifdef HAVE_BGQ_MPIX
#include "topobgq.inc"
#endif
#ifdef HAVE_MPI_COMM_SPLIT_TYPE
#include "topompisplit.inc"
#endif
#ifdef HAVE_MPI_WITH_NODENAME
#include "topompiname.inc"
#endif

/* Unknown interconnect topology.  Provide no information */
static int topoiGetNodeTopoInfoDEFAULT(topoinfo_t *topo)
{
  return 0;
}


/* Get information about the node (often requires different routines
   than the network topology, so these are setup separately.
   If the code to get the topology also provides the node info, the
   name PROVIDED_GETNODEINFO will be defined

   Error returns are:
   0  - no error, successful return
   -1 - error (such as memory allocation)
   1  - no information available
*/
#ifdef HAVE_HWLOC_H
#include "nodehwloc.inc"
#endif

#ifdef HAVE_SCHED_GETCPU
#include "nodegetcpu.inc"
#endif

#ifdef HAVE_SCHED_GETAFFINITY
#include "nodegetaff.inc"
#endif


static int topoiGetNodeInfoDEFAULT(int isMultithreaded, topoinfo_t *topo)
{
  /* Nothing available */
  return 1;
}

#if 0
 elif defined(USE_MPI_WITH_NODENAME)
 include "nodempiname.inc"
 elif defined(HAVE_MPI_COMM_SPLIT_TYPE)
 include "nodempisplit.inc"
#endif


static void topoiInitFns(void)
{
    static int initialized = 0;
    if (initialized) return;
    initialized = 1;
#ifdef HAVE_CRAY_RCA
    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoCRAY;
    topoNetFns[netFnsNum].descString      = "Cray RCA";
    topoNetFns[netFnsNum].id              = TOPO_CRAY_RCA;
#endif
#ifdef HAVE_BGQ_MPIX
    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoBGQ;
    topoNetFns[netFnsNum].descString      = "IBM BGQ MPIX";
    topoNetFns[netFnsNum].id              = TOPO_BGQ_MPIX;
    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo   = topoiGetNodeInfoBGQ;
    topoNodeFns[nodeFnsNum].descString    = "IBM BGQ MPIX";
    topoNodeFns[nodeFnsNum].id            = TOPO_BGQ_MPIX;
#endif
#ifdef HAVE_MPI_COMM_SPLIT_TYPE
    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoMPISPLIT;
    topoNetFns[netFnsNum].descString      = "Generic with MPI_Comm_split_type";
    topoNetFns[netFnsNum].id              = TOPO_MPISPLIT;
#endif
#ifdef HAVE_MPI_WITH_NODENAME
    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoNODENAME;
    topoNetFns[netFnsNum].descString     = "Generic with nodename";
    topoNetFns[netFnsNum].id              = TOPO_NODENAME;
#endif
    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoDEFAULT;
    topoNetFns[netFnsNum].descString      = "Generic default";
    topoNetFns[netFnsNum].id              = TOPO_DEFAULT;

    netFnsNum++;
    topoNetFns[netFnsNum].getNodeTopoInfo = topoiGetNodeTopoInfoDEBUG;
    topoNetFns[netFnsNum].descString      = "Debug only";
    topoNetFns[netFnsNum].id              = TOPO_DEBUG;

#ifdef HAVE_HWLOC_H
    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo      = topoiGetNodeInfoHWLOC;
    topoNodeFns[nodeFnsNum].descString       = "HWLOC";
    topoNodeFns[nodeFnsNum].id               = TOPO_HWLOC;
#endif

#ifdef HAVE_SCHED_GETCPU
    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo      = topoiGetNodeInfoGETCPU;
    topoNodeFns[nodeFnsNum].descString       = "sched_getcpu";
    topoNodeFns[nodeFnsNum].id               = TOPO_GETCPU;
#endif

#ifdef HAVE_SCHED_GETAFFINITY
    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo      = topoiGetNodeInfoGETAFFINITY;
    topoNodeFns[nodeFnsNum].descString       = "getaffinity";
    topoNodeFns[nodeFnsNum].id               = TOPO_GETAFFINITY;
#endif

    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo      = topoiGetNodeInfoDEFAULT;
    topoNodeFns[nodeFnsNum].descString       = "Generic Node Information";
    topoNodeFns[nodeFnsNum].id               = TOPO_GENERIC;

    nodeFnsNum++;
    topoNodeFns[nodeFnsNum].getNodeInfo      = topoiGetNodeInfoDEBUG;
    topoNodeFns[nodeFnsNum].descString       = "Debug";
    topoNodeFns[nodeFnsNum].id               = TOPO_DEBUG;

    /* The default mechanism is the first one loaded.  There is always
       a default, so index 0 will always be valid */
    nodeFnsIdx = 0;
    netFnsIdx  = 0;
}

/*
 * Internal routine to setup information about the node.
 */
static int topoiSetupNodeInfo(topoinfo_t *topo)
{
    topoentry_t *e;
    int         setupComms = 0, wrank, color, i;
    MPI_Comm    nodecomm, listcomm;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Check whether there is anything to do */
    if (topo->nnodes != -1) return 0;

    e = topo->info;
    while (e) {
	if (e->topoType == TOPO_NODELIST) {
	    /* Use the node number of this node as the color
	       in creating the nodecomm */
	    /* Note, however, that nodelist has usually set topo->nnodes and
	       the other fields */
	    MPI_Comm_split(MPI_COMM_WORLD, e->coords.coords[0], wrank,
			   &nodecomm);

	    topo->nodenum = e->coords.coords[0];
	    topo->nnodes  = e->maxcoords.coords[0];
	    setupComms = 1;
	    break;
	}
	else if (e->topoType == TOPO_TORUS || e->topoType == TOPO_MESH) {
	    color = e->coords.coords[0] - e->mincoords.coords[0];
	    for (i=1; i<e->dim; i++)
		color = e->coords.coords[i] - e->mincoords.coords[i] + 
		    (e->maxcoords.coords[i-1] - e->mincoords.coords[i-1] + 1)*color;
	    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &nodecomm);
	    setupComms = 1;
	    break;
	}
	e = e->next;
    }

    if (setupComms) {
	int *ranks, nrank, nsize, nodenums[2];
	MPI_Group   gn, gw;

	MPI_Comm_rank(nodecomm, &nrank);
	MPI_Comm_size(nodecomm, &nsize);
	fflush(stdout);
	MPI_Comm_group(MPI_COMM_WORLD, &gw);
	MPI_Comm_group(nodecomm, &gn);
	ranks       = (int *)malloc(nsize * sizeof(int));
	topo->ranks = (int *)malloc(nsize * sizeof(int));
	for (i=0; i<nsize; i++) ranks[i] = i;
	MPI_Group_translate_ranks(gn, nsize, ranks, gw, topo->ranks);
	free(ranks);
	MPI_Group_free(&gn);
	MPI_Group_free(&gw);

	color = (nrank == 0) ? 0 : MPI_UNDEFINED;
	MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &listcomm);
	if (color == 0) {
	    MPI_Comm_rank(listcomm, &nodenums[0]);
	    MPI_Comm_size(listcomm, &nodenums[1]);
	}
	MPI_Bcast(nodenums, 2, MPI_INT, 0, nodecomm);
	topo->nranks         = nsize;
	topo->nodenum        = nodenums[0];
	topo->nnodes         = nodenums[1];
	topo->nodeComm       = nodecomm;
	topo->nodeLeaderComm = listcomm;
    }

    return 0;
}
