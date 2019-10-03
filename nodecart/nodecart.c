/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpi.h"
#include "nodecart.h"

/*
 * These routines implement a virtual Cartesian process topology by grouping
 * processes on a single node together.
 *
 * These are intended to illustate an alternative implementation of
 * MPI_Cart_create that while not optimal, is better than current implmentations
 * and is simple and cheap to implement.
 *
 * MPIX_Nodecart_create  - Like MPI_Cart_create
 * MPIX_Nodecart_shift   - Like MPI_Cart_shift
 * MPIX_Nodecart_coords  - Like MPI_Cart_coords
 * MPIX_Nodecart_get     - Like MPI_Cart_get
 * MPIX_Nodecart_rank    - Like MPI_Cart_rank
 * MPIX_Nodecart_dim_get - Like MPI_Cartdim_get
 * MPIX_Nodecart_sub     - Like MPI_Cart_sub
 *
 * No counterpart to MPI_Cart_map has been implemented
 *
 *
 * MPIX_Nodecomm_create - Convenience routine to create communicators for
 *    the node and a communicator of node leaders.
 *
 * Information about communication
 *
 *        MPIX_GetOffNodeCounts - Given a set of ranks, determine how many of
 *              those ranks target the same node as the caller, and how
 *              many are offnode.
 *
 * Convenience routines to simplify example codes.  These create and work
 * with a "communication info structure", comminfo_t, that contains a
 * communicator, ranks of the neighboring processes, and information about
 * the process mesh.  The structure is defined in mpimesh.h
 *
 * Finally, if compiled with -DBUILD_TEST, this file creates a main program
 * that can be used as a simple test.
 */

/* Naming convention
   There are several layers of hierarchy; a uniform naming convention is
   helpful.
   The heirarchy is:
       comm or cartcomm
       node
       socket
       core

   The number of each type of object is given by num<object>, as in numnodes
   (the number of nodes).
   The rank of the calling process in each of these is <object>rank, as
   in noderank (range [0,numnodes)) or socketrank (range [0,numsockets))
   The number of processes within an object is given by <object>size, as
   in nodesize, and the rank of a process within that object is given by
   rankIn<object> (range [0,<object>size)).

   Only the top two levels are required (comm and node). Some systems
   may make information available about the socket and core, and this
   can be used to further improve the process mapping.  These follow
   the same naming scheme.

   A socket is made up of cores (range [0,numcores)) and the rankInsocket
   reflects the core number in the socket on which the process is running.
   Since the core number is sometimes available, this is stored as well.

   Also need to know the core in socket and socket in node.
   Have a separate coordinate array:
      node#/socket#/core#
   if only have nodes, then socket# is always 0, and core# will be within
   nodesize
 */

/* Information about the placement of this process on a node. This
   assumes that the OS does not (often) migrate processes to different
   cores or chips, which is a good assumption on an HPC system but may not
   be as relevant on other systems */
typedef struct {
    /* Extra information about this process on the node. -1 == N/A */
    int nchips,           /* Number of separate processor chips on this node */
	chiprank,         /* rank (in [0,nchips)) of chip */
	rankinchip,       /* rank of core (in [0,ncores)) of core */
	chipnumber,       /* Number of the chip on which this process was
			     running when this field was set */
	corenumber;       /* Number of the core on which this process was
			     running when this field was set */
} procnode_t;

/* A structure with just the information about the node topology, as
   induced by the communicator used to create a new nodecart (or similar)
   communicator. This does *not* include the nodecomm or leadercomm. */
typedef struct {
    int      nnodes,      /* Number of nodes (size of leadercomm where
			     defined) */
	noderank,         /* Rank of this node in leader comm (valid for
			     all processes, not just the one in leadercomm) */
	rankinnode,       /* Rank of this process in the node */
	nodesize,
	minnsize,         /* Minimum size of nodecomm over all nodes */
	maxnsize;         /* Maximum size of nodecomm over all nodes */
} nodetopo_t;

typedef struct {
    MPI_Comm nodecomm,    /* Communicator of processes on the same node */
	leadercomm;       /* Communicator of the processes with rank==0 in
			     nodecomm */
    int      nnodes,      /* Number of nodes (size of leadercomm where
			     defined) */
	noderank;         /* Rank of this node in leader comm (valid for
			     all processes, not just the one in leadercomm) */
} nodeinfo_t;

#define MAX_DIM 5

typedef struct {
    MPI_Comm nodecomm, leadercomm;
    int      interdims[MAX_DIM], intradims[MAX_DIM],
	intercoords[MAX_DIM], intracoords[MAX_DIM],
	dims[MAX_DIM], coords[MAX_DIM],
	periodic[MAX_DIM];
    int ndim, meshsize;   /* Number of dimensions and product of all
			     dimensions (size of the mesh) */
} nodecart_t;

static int nodeinfoKeyval = MPI_KEYVAL_INVALID,
    nodecartKeyval = MPI_KEYVAL_INVALID,
    procnodeKeyval = MPI_KEYVAL_INVALID,
    nodetopoKeyval = MPI_KEYVAL_INVALID;

/* Set to true to make the code use socket and core information in
   determining the coordinates of a process within a node */
static int useSocketInfo = 0;

int MPIX_SetUseSocket(int flag)
{
    int oldval = useSocketInfo;
    useSocketInfo = flag;
    return oldval;
}

/* Define PRIVATE to either null or static; null to allow testing of
   internal functions */
#ifndef PRIVATE
#define PRIVATE static
#endif

static int nodecartDelFn(MPI_Comm comm, int keyval, void *attr, void *estate);
static int nodeinfoDelFn(MPI_Comm comm, int keyval, void *attr, void *estate);
static int procnodeCopyFn(MPI_Comm comm, int keyval, void *estate,
			  void *attr_in, void *attr_out, int *flag);
static int nodetopoCopyFn(MPI_Comm comm, int keyval, void *estate,
			  void *attr_in, void *attr_out, int *flag);
static int procnodeDelFn(MPI_Comm comm, int keyval, void *attr, void *estate);
static int nodetopoDelFn(MPI_Comm comm, int keyval, void *attr, void *estate);
PRIVATE int decompDims(int ndims, const int dims[], int nsize,
		       int intradims[], int interdims[]);
PRIVATE int findInInterval(int val, int offset, int high);
PRIVATE void rankToCoords(int ndims, const int dims[], int rank, int coords[]);
PRIVATE void coordsToRank(int ndims, const int dims[], const int coords[],
			  int *rank);
PRIVATE void rankShift(int ndims, const int dims[], const int coords[],
		       const int periodic[], int rank,
		       int direction, int disp, int *rsource, int *rdest);
PRIVATE int factor(int n, int *nf_ptr, int factors[], int powers[]);
static void decompDimsError(int ndims, const int dims[], nodecart_t *ninfo);

static void GetChipCoreNumber(MPI_Comm comm, int chipnum, int corenum,
			      int *chiprank, int *corerank);

/* Temporary debugging */
static int cvar_nodecart_ppn = 0;
static int cvar_nodecart_verbose = 0;
static FILE *vfp = 0;

/*@ MPIX_Nodecomm_create - Create communicators for intra- and inter-node
 communication

 Input Parameter:
. incomm - Input communicator

 Output Parameters:
+ nodecomm   - Communicator containing all processes on the same node
. leadercomm - A communicator of the leaders on each node, defined as
 the processes with rank 0 in nodecomm.  If the calling process is not
 in rank 0 in nodecomm, sets to 'MPI_COMM_NULL'
. noderank   - Rank of node in 'leadercomm', valid on all processes
- nnodes     - Number of nodes.

 Notes:
 This routine is collective over 'incomm'.
 Information is saved in an attribute, so that information about the
 node and leader communicators can be recovered from an attribute (which
 is private to this module).

 There is no requirement that every node have the same number of processes.

 Optionally, if a suitable API is available, this routine will also determine
 and store information about the number of processor chips on the node of
 the calling process, as well as which chip an which core on that chip the
 calling process was running on when this routine was called.  If this
 information is not available, a value of '-1' is stored for those fields
 in the nodeinfo_t structure.  This data is not returned by this routine,
 but is available to the nodecart routines.
  @*/
int MPIX_Nodecomm_create(MPI_Comm incomm, MPI_Comm *nodecomm,
			 MPI_Comm *leadercomm, int *noderank, int *nnodes)
{
    int         rank, nrank, color;
    int         sz[2];
    nodeinfo_t *ninfo=0;
    procnode_t *pinfo=0, *pinfo2=0;
    nodetopo_t *tinfo=0;

    /* Check to see if we already have this information */
    if (nodeinfoKeyval != MPI_KEYVAL_INVALID) {
	int        flag;
	nodeinfo_t *nninfo;
	MPI_Comm_get_attr(incomm, nodeinfoKeyval, &nninfo, &flag);
	if (flag) {
	    *nodecomm   = nninfo->nodecomm;
	    *leadercomm = nninfo->leadercomm;
	    *nnodes     = nninfo->nnodes;
	    return MPI_SUCCESS;
	}
    }

    /* We do not have this info.  Determine it */
    MPI_Comm_rank(incomm, &rank);

    /* Check first for the special case used for debugging */
    if (cvar_nodecart_ppn > 0) {
	color = rank / cvar_nodecart_ppn;
	MPI_Comm_split(incomm, color, rank, nodecomm);
    }
    else {
	MPI_Comm_split_type(incomm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
			    nodecomm);
    }
    MPI_Comm_rank(*nodecomm, &nrank);
    if (nrank == 0) color = 0;
    else            color = MPI_UNDEFINED;
    MPI_Comm_split(incomm, color, rank, leadercomm);

    /* Return the number of nodes to all processes */
    if (color == 0) {
	MPI_Comm_size(*leadercomm, nnodes);
	MPI_Comm_rank(*leadercomm, noderank);
	sz[0] = *nnodes;
	sz[1] = *noderank;
    }
    /* Note that this creates nnodes concurrent, disjoint broadcasts */
    MPI_Bcast(sz, 2, MPI_INT, 0, *nodecomm);
    if (color != 0) {
	*nnodes   = sz[0];
	*noderank = sz[1];
    }

    /* Save information as an attribute on incomm */

    /* Create the keyval if it doesn't already exist.
       We do not define a copy function because we are returning the same
       nodecomm (and leadercomm).  If we made a copy, we would either need
       to dup those communicators or maintain a reference count.  Either
       is possible, but until there is a significant use case, we choose
       the simpler option and do not provide a copy option */
    if (nodeinfoKeyval == MPI_KEYVAL_INVALID) {
	MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, nodeinfoDelFn,
			       &nodeinfoKeyval, NULL);
    }
    /* Create the structure into which we'll record the info */
    ninfo = (nodeinfo_t *)malloc(sizeof(nodeinfo_t));
    if (!ninfo) {
	fprintf(stderr, "Unable to allocate memory for nodeinfo_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ninfo->nodecomm   = *nodecomm;
    ninfo->leadercomm = *leadercomm;
    ninfo->nnodes     = *nnodes;
    ninfo->noderank   = *noderank;

    /* Temporary - also attach a separate procnode_t attribute that is
       copied to the derivative communicators */
    if (procnodeKeyval == MPI_KEYVAL_INVALID) {
	MPI_Comm_create_keyval(procnodeCopyFn, procnodeDelFn, &procnodeKeyval,
	    NULL);
    }
    pinfo = (procnode_t *)malloc(sizeof(procnode_t));
    if (!pinfo) {
	fprintf(stderr, "Unable to allocate memory for procnode_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* Set the chip info, if possible.  If not, the values are set to -1 */
    MPIX_GetSocketAndCPU(&pinfo->nchips, &pinfo->chipnumber,
			 &pinfo->corenumber);
    /* Both the rankinchip and chiprank need to be compressed to be in the
       the range of [0,size) */
    /* Algorithm:
       By node, gather the core numbers. sort and use index as chiprank.
       Save core# in the structure.
       Ditto for chip
    */
    GetChipCoreNumber(ninfo->nodecomm, pinfo->chipnumber, pinfo->corenumber,
		      &pinfo->chiprank, &pinfo->rankinchip);

    MPI_Comm_set_attr(incomm, procnodeKeyval, pinfo);
    /* To set on another comm, must make a copy, because the attr_destroy
       routine frees the memory */
    pinfo2 = (procnode_t *)malloc(sizeof(procnode_t));
    if (!pinfo2) {
	fprintf(stderr, "Unable to allocate memory for procnode_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    *pinfo2 = *pinfo;
    MPI_Comm_set_attr(ninfo->nodecomm, procnodeKeyval, pinfo2);
    /* Also set on leadercomm? */

    /* Attach the information to the communicator */
    MPI_Comm_set_attr(incomm, nodeinfoKeyval, ninfo);

    /* Temporary - also attach a separate nodetopo_t attribute that is
       copied to the derivative communicators */
    if (nodetopoKeyval == MPI_KEYVAL_INVALID) {
	MPI_Comm_create_keyval(nodetopoCopyFn, nodetopoDelFn, &nodetopoKeyval,
			       NULL);
    }
    tinfo = (nodetopo_t *)malloc(sizeof(nodetopo_t));
    if (!tinfo) {
	fprintf(stderr, "Unable to allocate memory for nodetopo_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    tinfo->nnodes     = ninfo->nnodes;
    tinfo->noderank   = ninfo->noderank;
    /* Check size of each nodecomm - can take the max of (size,-size)
     */
    MPI_Comm_rank(*nodecomm, &tinfo->rankinnode);
    MPI_Comm_size(*nodecomm, &tinfo->nodesize);
    MPI_Comm_size(*nodecomm, &sz[0]);
    sz[1] = -sz[0];
    MPI_Allreduce(MPI_IN_PLACE, sz, 2, MPI_INT, MPI_MAX, incomm);
    tinfo->minnsize   = -sz[1];
    tinfo->maxnsize   = sz[0];
    MPI_Comm_set_attr(incomm, nodetopoKeyval, tinfo);
    /* Also set on nodecomm and leadercomm? */

    return MPI_SUCCESS;
}

/* socketcomm == MPI_COMM_NULL if either no socket info or only one socket */
int MPIX_Socketcomm_create(MPI_Comm nodecomm, MPI_Comm *socketcomm)
{
    procnode_t *pinfo;
    int        flag;

    if (procnodeKeyval == MPI_KEYVAL_INVALID) {
	fprintf(stderr, "No socket information available on nodecomm\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_get_attr(nodecomm, procnodeKeyval, &pinfo, &flag);
    if (!flag || !pinfo) {
	fprintf(stderr, "Unable to access socket info attribute\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (pinfo->nchips < 2) {
	*socketcomm = MPI_COMM_NULL;
	return MPI_SUCCESS;
    }
    MPI_Comm_split(nodecomm, pinfo->chiprank, pinfo->rankinchip, socketcomm);
    return MPI_SUCCESS;
}

/*@ MPIX_Comm_dims_create - Creates a division of processors within a given
 communicator in a Cartesian grid

Input Parameters:
+ comm  - Communicator (handle)
. nnodes - number of nodes in a grid (integer)
- ndims - number of Cartesian dimensions (integer)

Input/Output Parameters:
. dims - integer array of size  'ndims' specifying the number of nodes in each
 dimension. All entries (for 'ndims') are filled in; unlike 'MPI_Dims_create',
 values in 'dims' do not set the size for that dimension.

Notes:
Unlike 'MPI_Dims_create', this is a collective routine over 'comm'.  This
permits the implementation to perform collective communication operations to
determine good decompositions.

@*/
int MPIX_Comm_dims_create(MPI_Comm comm, int nnodes, int ndims, int dims[])
{
    int i, interdims[MAX_DIM], intradims[MAX_DIM], nsize;
    MPI_Comm nodecomm, leadercomm;
    /* Get information about the node.  If the information is not already
       attached to comm, we have to go get it (which requires collective
       communication) */
    nodecomm = MPI_COMM_NULL;
    if (nodeinfoKeyval != MPI_KEYVAL_INVALID) {
	int        flag;
	nodeinfo_t *nninfo;
	MPI_Comm_get_attr(comm, nodeinfoKeyval, &nninfo, &flag);
	if (flag) {
	    nodecomm   = nninfo->nodecomm;
	}
    }
    /* If we didn't get the info, create it */
    if (nodecomm == MPI_COMM_NULL) {
	int err, numnodes, noderank;
	err = MPIX_Nodecomm_create(comm, &nodecomm,
				   &leadercomm, &noderank, &numnodes);
    }

    for (i=0; i<ndims; i++) {
	intradims[i] = 0;
	interdims[i] = 0;
    }
    MPI_Comm_rank(nodecomm, &nsize);
    /* FIXME: Using Dims_create will often produce suboptimal choices.
       It is better to use code such as that in decompDims. */
    if (nsize <= nnodes && (nnodes % nsize) == 0) {
	MPI_Dims_create(nsize, ndims, intradims);
	MPI_Dims_create(nnodes/nsize, ndims, interdims);

	for (i=0; i<ndims; i++) {
	    dims[i] = interdims[i] * intradims[i];
	}
    }
    else {
	/* Error with this node size.  In this case, reset and just use
	   MPI_Dims_create on csize */
	MPI_Dims_create(nnodes, ndims, dims);
    }

    return MPI_SUCCESS;
}

/*@ MPIX_Nodecart_create - Create a Cartesian communicator, using information
  about which nodes are on the same process.

 Input Parameters:
+ comm_old - input communicator (handle)
. ndims - number of dimensions of Cartesian grid (integer)
. dims - integer array of size ndims specifying the number of processes in
  each dimension
. periods - logical array of size ndims specifying whether the grid is
  periodic (true) or not (false) in each dimension
- reorder - ranking may be reordered (true) or not (false) (logical)

Output Parameters:
. comm_cart - communicator with new Cartesian topology (handle)

 Notes:
 Like 'MPI_Cart_create', the values for dims must be provided on input.
 An option to consider is to allow zero values for elements of dims,
 and then let this routine choose the dimensions to best fit the underlying
 physical hardware.

 The routine 'MPIX_Comm_dims_create' may be used to determine good values
 for 'dims'.

 The control variable 'cvar_nodecart_verbose' may be set to a positive value
 to cause output to be generated about this routine''s operation.
  @*/
int MPIX_Nodecart_create(MPI_Comm comm_old, int ndims, const int dims[],
             const int periods[], int reorder, MPI_Comm *comm_cart)
{
    int        nnodes, noderank, nrank, nsize, rr;
    int        i, inrank, insize, tsize, rc;
    nodecart_t *ninfo;

    if (ndims > MAX_DIM) {
	fprintf(stderr, "Maximum number of dimensions is %d\n", MAX_DIM);
	fflush(stderr);
	return MPI_ERR_OTHER;
    }

    /* Create the structure into which we'll record the info */
    ninfo = (nodecart_t *)malloc(sizeof(nodecart_t));
    if (!ninfo) {
	fprintf(stderr, "Unable to allocate memory for nodecart_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Get nodeinfo about comm_old */
    ninfo->nodecomm = MPI_COMM_NULL;
    if (nodeinfoKeyval != MPI_KEYVAL_INVALID) {
	int        flag;
	nodeinfo_t *nninfo;
	MPI_Comm_get_attr(comm_old, nodeinfoKeyval, &nninfo, &flag);
	if (flag) {
	    ninfo->nodecomm   = nninfo->nodecomm;
	    ninfo->leadercomm = nninfo->leadercomm;
	    nnodes            = nninfo->nnodes;
	    noderank          = nninfo->noderank;
	}
    }
    /* If we didn't get the info, create it */
    if (ninfo->nodecomm == MPI_COMM_NULL) {
	int err;
	err = MPIX_Nodecomm_create(comm_old, &ninfo->nodecomm,
				   &ninfo->leadercomm, &noderank, &nnodes);
	if (err) return err;
    }

    MPI_Comm_size(comm_old, &insize);
    if (cvar_nodecart_verbose > 1) {
	MPI_Comm_rank(comm_old, &inrank);
	if (inrank == 0) {
	    if (!vfp) vfp = stdout;
	    fprintf(vfp, "Nnodes = %d\n", nnodes);
	    fflush(vfp);
	}
    }

    /* Create a decomposition of the dims into intranode and internode
       dimensions */
    MPI_Comm_size(ninfo->nodecomm, &nsize);
    rc = decompDims(ndims, dims, nsize, ninfo->intradims, ninfo->interdims);
    if (rc) {
	decompDimsError(ndims, dims, ninfo);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* First, get the coordinate of this process in the decomposition of the
       nodes */
    /* Every process on the node concurrently computes the coordinates
       of this node in the mesh of nodes.  This is faster than having
       one process compute and distribute */
    rankToCoords(ndims, ninfo->interdims, noderank, ninfo->intercoords);
    ninfo->ndim           = ndims;

    /* Get the coordinates of this process on the node */
    MPI_Comm_rank(ninfo->nodecomm, &nrank);
    if (useSocketInfo && procnodeKeyval != MPI_KEYVAL_INVALID) {
	int         socketcoords[MAX_DIM], intracoords[MAX_DIM];
	int         intersocket[MAX_DIM], intrasocket[MAX_DIM];
	procnode_t *pinfo;
	int         socketsize,   /* Number of cores on a socket */
	    flag;
	/* Try to get the socket info */
	MPI_Comm_get_attr(comm_old, procnodeKeyval, &pinfo, &flag);
	if (flag && pinfo && pinfo->nchips > 1) {
	    socketsize = nsize / pinfo->nchips;
	    decompDims(ndims, ninfo->intradims, socketsize,
		       intrasocket, intersocket);
	    if (pinfo->rankinchip >= socketsize) {
		/* Sanity check failed */
		fprintf(stderr, "core rank %d >= cores per socket %d\n",
			pinfo->rankinchip, socketsize);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	    rankToCoords(ndims, intersocket, pinfo->chiprank,   socketcoords);
	    rankToCoords(ndims, intrasocket, pinfo->rankinchip, intracoords);
	    for (i=0; i<ndims; i++) {
		ninfo->intracoords[i] = socketcoords[i]*intrasocket[i] +
		    intracoords[i];
	    }
	}
	else {
	    /* Use noderank-only code */
	    rankToCoords(ndims, ninfo->intradims, nrank, ninfo->intracoords);
	}
    }
    else {
	rankToCoords(ndims, ninfo->intradims, nrank, ninfo->intracoords);
    }

    /* Compute total sizes and coordinates in the full grid */
    tsize = 1;
    for (i=0; i<ndims; i++) {
	ninfo->dims[i]     = ninfo->interdims[i] * ninfo->intradims[i];
	ninfo->coords[i]   = ninfo->intercoords[i] * ninfo->intradims[i] +
	    ninfo->intracoords[i];
	tsize             *= ninfo->dims[i];
	ninfo->periodic[i] = periods[i];
    }
    ninfo->meshsize = tsize;

    /* Check that the decomposition is successful - that the total size
       matches the input communicator */
    if (insize != tsize) {
	fprintf(stderr,
	"ERROR: nodecart requires the same number of processes on each node\n");
	fflush(stderr);
	return MPI_ERR_OTHER;
    }

    /* Now, compute the rank of this process using the two decompositions. */
    coordsToRank(ndims, ninfo->dims, ninfo->coords, &rr);

    if (cvar_nodecart_verbose > 1) {
	if (!vfp) vfp = stdout;
	fprintf(vfp, "rr = %d\n", rr);
	fflush(vfp);
    }

    /* Info on this decomposition can be accessed through the MPIX_Nodecart_xxx
       functions */

    /* Create the new communicator */
    MPI_Comm_split(comm_old, 0, rr, comm_cart);

    /* Attach the information to the communicator */
    if (nodecartKeyval == MPI_KEYVAL_INVALID) {
	MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, nodecartDelFn,
			       &nodecartKeyval, NULL);
    }
    MPI_Comm_set_attr(*comm_cart, nodecartKeyval, ninfo);
    /* Propagate information about the topology to the new communicator */
    if (nodetopoKeyval != MPI_KEYVAL_INVALID) {
	int flag, fflag;
	nodetopo_t *told, *tnew;
	MPI_Comm_get_attr(comm_old, nodetopoKeyval, &told, &flag);
	if (flag) {
	    nodetopoCopyFn(comm_old, nodetopoKeyval, NULL, told, &tnew,
			   &fflag);
	    if (fflag)
		MPI_Comm_set_attr(*comm_cart, nodetopoKeyval, tnew);
	}
    }

    return MPI_SUCCESS;
}

/*@ MPIX_Nodecart_shift - Returns the shifted source and destination
                 ranks, given a shift direction and amount

Input Parameters:
+ comm - communicator with Cartesian structure (handle)
. direction - coordinate dimension of shift (integer)
- disp - displacement (> 0: upwards shift, < 0: downwards shift) (integer)

Output Parameters:
+ rank_source - rank of source process (integer)
- rank_dest - rank of destination process (integer)

Notes:
The 'direction' argument is in the range '[0,n-1]' for an n-dimensional
Cartesian mesh.
@*/
int MPIX_Nodecart_shift(MPI_Comm comm, int direction, int disp,
			int *rank_sources, int *rank_dest)
{
    int        rank, flag;
    nodecart_t *ninfo;

    MPI_Comm_rank(comm, &rank);

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }

    rankShift(ninfo->ndim, ninfo->dims, ninfo->coords, ninfo->periodic,
	      rank, direction, disp, rank_sources, rank_dest);

    if (cvar_nodecart_verbose > 1) {
	if (!vfp) vfp = stdout;
	fprintf(vfp, "coords[%d] = %d, dims = %d, low = %d high = %d\n",
	       direction, ninfo->coords[direction], ninfo->dims[direction],
		*rank_sources, *rank_dest);
	fflush(vfp);
    }

    return MPI_SUCCESS;
}

/*@
MPIX_Nodecart_coords - Determines process coords in Cartesian topology given
                  rank in group

Input Parameters:
+ comm - communicator with Cartesian structure (handle)
. rank - rank of a process within group of 'comm' (integer)
- maxdims - length of vector 'coords' in the calling program (integer)

Output Parameters:
. coords - integer array (of size 'ndims') containing the Cartesian
  coordinates of specified process (integer)
  @*/
int MPIX_Nodecart_coords(MPI_Comm comm, int rank, int maxdims, int coords[])
{
    int        flag;
    nodecart_t *ninfo;

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }

    rankToCoords(ninfo->ndim, ninfo->dims, rank, coords);

    return MPI_SUCCESS;
}

/*@
MPIX_Nodecart_rank - Determines process rank in communicator given Cartesian
                location

Input Parameters:
+ comm - communicator with Cartesian structure (handle)
- coords - integer array (of size 'ndims', the number of dimensions of
    the Cartesian topology associated with 'comm') specifying the Cartesian
  coordinates of a process

Output Parameters:
. rank - rank of specified process (integer)

Notes:
 Out-of-range coordinates are erroneous for non-periodic dimensions.
 @*/
int MPIX_Nodecart_rank(MPI_Comm comm, const int coords[], int *rank)
{
    int        flag;
    nodecart_t *ninfo;

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }

    coordsToRank(ninfo->ndim, ninfo->dims, coords, rank);

    return MPI_SUCCESS;
}

/*@
MPIX_Nodecart_sub - Partitions a communicator into subgroups which
               form lower-dimensional Cartesian subgrids

Input Parameters:
+ comm - communicator with Cartesian structure (handle)
- remain_dims - the  'i'th entry of remain_dims specifies whether the 'i'th
dimension is kept in the subgrid (true) or is dropped (false) (logical
vector)

Output Parameters:
. newcomm - communicator containing the subgrid that includes the calling
process (handle)

Note:
 Needed for the snap CORAL benchmark - they use this instead of
   'MPI_Cart_shift' .
@*/
int MPIX_Nodecart_sub(MPI_Comm comm, const int remain[], MPI_Comm *subcomm)
{
    int        i, flag;
    int        newdims, meshsize, color, k;
    int        noderank, inrank, nrank, rank;
    nodecart_t *ninfo, *ninfonew;

    /* Algorithm: Using the coordinates where remain == FALSE, compute a
       rank using the coords to rank.  Perform comm_split using that
       rank as the color, and the original rank as rank.
       Create the nodeinfo data using the remaining dimensions from the
       input

       Must update nodecomm and leadercomm, since both are likely to
       be changed by the removal of dimensions.

       Note that we can't simply create a communicator with the remaining
       dimensions and then call Nodecart_create on that because we must
       retain the original mapping of all processes
    */

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }

    MPI_Comm_rank(comm, &inrank);
    /* How many remaining dimensions are there? Also compute the color on
       which to split */
    newdims  = 0;
    meshsize = 1;
    color    = 0;
    for (i=0; i<ninfo->ndim; i++) {
	if (remain[i]) {
	    newdims++;
	    meshsize *= ninfo->dims[i];
	}
	else {
	    /* Row major for simplicity */
	    color = color * ninfo->dims[i] + ninfo->coords[i];
	}
    }

    /* Create the subcommunicator */
    MPI_Comm_rank(comm, &inrank);
    MPI_Comm_split(comm, color, inrank, subcomm);

    /* Create the new nodeinfo properties */
    /* Create the structure into which we'll record the info */
    ninfonew = (nodecart_t *)malloc(sizeof(nodecart_t));
    if (!ninfonew) {
	fprintf(stderr, "Unable to allocate memory for nodecart_t\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    ninfonew->ndim = newdims;
    MPI_Comm_rank(ninfo->nodecomm, &noderank);
    MPI_Comm_split(ninfo->nodecomm, color, noderank, &ninfonew->nodecomm);
    /* Need to create the leadercomm for each subcomm, since in most cases,
       the subcomm will have no processes in the original leadercomm */
    k = 0;
    for (i=0; i<ninfo->ndim; i++) {
	if (remain[i]) {
	    ninfonew->interdims[k]   = ninfo->interdims[i];
	    ninfonew->intradims[k]   = ninfo->intradims[i];
	    ninfonew->dims[k]        = ninfo->dims[i];
	    ninfonew->intercoords[k] = ninfo->intercoords[i];
	    ninfonew->intracoords[k] = ninfo->intracoords[i];
	    ninfonew->coords[k]      = ninfo->coords[i];
	    ninfonew->periodic[k]    = ninfo->periodic[i];
	    k++;
	}
    }
    MPI_Comm_rank(ninfonew->nodecomm, &nrank);
    MPI_Comm_rank(*subcomm, &rank);
    if (nrank == 0) color = 0;
    else            color = MPI_UNDEFINED;
    MPI_Comm_split(*subcomm, color, rank, &ninfonew->leadercomm);

    return MPI_SUCCESS;
}

/*@

MPIX_Nodecartdim_get - Retrieves Cartesian topology information
                  associated with a communicator

Input Parameters:
. comm - communicator with Cartesian structure (handle)

Output Parameters:
. ndims - number of dimensions of the Cartesian structure (integer)
@*/
int MPIX_Nodecart_dim_get(MPI_Comm comm, int *ndims)
{
    int        flag;
    nodecart_t *ninfo;

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }
    *ndims = ninfo->ndim;
    return MPI_SUCCESS;
}

/*@

MPIX_Nodecart_get - Retrieves Cartesian topology information associated with a
               communicator

Input Parameters:
+ comm - communicator with Cartesian structure (handle)
- maxdims - length of vectors  'dims', 'periods', and 'coords'
in the calling program (integer)

Output Parameters:
+ dims - number of processes for each Cartesian dimension (array of integer)
. periods - periodicity (true/false) for each Cartesian dimension
(array of logical)
- coords - coordinates of calling process in Cartesian structure
(array of integer)
@*/
int MPIX_Nodecart_get(MPI_Comm comm, int maxdims, int dims[], int periods[],
		      int coords[])
{
    int        flag, i;
    nodecart_t *ninfo;

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* The standard isn't clear what happens if maxdims < ninfo->ndim.
       We abort with an error message */
    if (maxdims < ninfo->ndim) {
	fprintf(stderr, "Nodecart_get: maxdims < ndims for comm\n");
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i=0; i<ninfo->ndim; i++) {
	dims[i]    = ninfo->dims[i];
	periods[i] = ninfo->periodic[i];
	coords[i]  = ninfo->coords[i];
    }
    return MPI_SUCCESS;
}

#if 0
/*@
MPIX_Nodecart_map - Maps process to Cartesian topology information

Input Parameters:
+ comm - input communicator (handle)
. ndims - number of dimensions of Cartesian structure (integer)
. dims - integer array of size 'ndims' specifying the number of processes in
  each coordinate direction
- periods - logical array of size 'ndims' specifying the periodicity
  specification in each coordinate direction

Output Parameters:
. newrank - reordered rank of the calling process; 'MPI_UNDEFINED' if
  calling process does not belong to grid (integer)

  @*/
int MPIX_Nodecart_map(MPI_Comm comm, int ndims, const int dims[],
		      const int periods[], int *newrank)
{
#error 'This routine is not yet implemented'
}
#endif

/*
 * Return the decomposition for each level.  Currently there are only
 * two levels.  dims[0..*ndim-1] is the intra-node dimensions,
 * dims[*ndim...2 **ndim-1] is the inter-node dimensions.

Input Parameter:
. comm - Communicator with 'nodecart' process topology

Output Parameters:
+ nlevel - The number of levels in the dcomposition (will be 2)
. ndim   - The number of dimensions in the decomposition
- dims   - Contains the dimensions of the decomposition on each level

Notes:
This routine is an alpha-test, and the interface might change.
 */
int MPIX_Nodecart_get_decomp(MPI_Comm comm, int *nlevel, int *ndim, int dims[])
{
    int        k, i, flag;
    nodecart_t *ninfo;

    MPI_Comm_get_attr(comm, nodecartKeyval, &ninfo, &flag);
    if (!ninfo || !flag) {
	fprintf(stderr, "Unable to access attribute on comm\n");
        fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1 );
    }

    *nlevel = 2;
    *ndim   = ninfo->ndim;
    k       = 0;
    for (i=0; i<ninfo->ndim; i++)
	dims[k++] = ninfo->intradims[i];
    for (i=0; i<ninfo->ndim; i++)
	dims[k++] = ninfo->interdims[i];

    return MPI_SUCCESS;
}

/* ------------------------------------------------------------------- */
/* Support routines */
/* Given a communicator and a list of target ranks, determine the number of:
 * On node and off node ranks.
 * Must also provide a communicator nodecomm
 * that contains all processes on the node of the caller.
 *
 * This code considers MPI_PROC_NULL as off-node. 
 */
int MPIX_GetOffNodeCounts(MPI_Comm comm, int nr, const int ranks[],
			  MPI_Comm nodecomm,
			  int *onnode, int *offnode)
{
    int *nranks;
    int on, off, i;
    MPI_Group cgroup, ngroup;

    nranks = (int *)malloc(nr * sizeof(int));
    if (!nranks) {
	fprintf(stderr, "Unable to allocate %d words\n", nr);
        fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_group(comm, &cgroup);
    MPI_Comm_group(nodecomm, &ngroup);
    MPI_Group_translate_ranks(cgroup, nr, ranks, ngroup, nranks);
    on  = 0;
    off = 0;
    for (i=0; i<nr; i++) {
	if (nranks[i] == MPI_PROC_NULL || nranks[i] == MPI_UNDEFINED) off++;
	else on ++;
    }
    /* For the number of distinct nodes, we need to map the ranks to
       different nodes */

    free(nranks);
    MPI_Group_free(&cgroup);
    MPI_Group_free(&ngroup);

    *onnode  = on;
    *offnode = off;

    return 0;
}

/* */
int MPIX_GetOffNodeCounts_X(MPI_Comm comm, int nr, const int ranks[],
			    int ncomm, MPI_Comm subcomm[],
			    int *onnode, int *offnode)
{
    int *nranks;
    int on, off, i, j;
    MPI_Group cgroup, ngroup;

    nranks = (int *)malloc(nr * sizeof(int));
    if (!nranks) {
	fprintf(stderr, "Unable to allocate %d words\n", nr);
        fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_group(comm, &cgroup);
    for (j=0; j<ncomm; j++) {
	MPI_Comm_group(subcomm[j], &ngroup);
	MPI_Group_translate_ranks(cgroup, nr, ranks, ngroup, nranks);
	on  = 0;
	off = 0;
	for (i=0; i<nr; i++) {
	    if (nranks[i] == MPI_PROC_NULL || nranks[i] == MPI_UNDEFINED) off++;
	    else on ++;
	}
	onnode[j]  = on;
	offnode[j] = off;
	/* For the number of distinct nodes, we need to map the ranks to
	   different nodes */
	MPI_Group_free(&ngroup);
    }
    free(nranks);

    MPI_Group_free(&cgroup);

    return 0;
}

/* Returns how many of the ranks are on the socket, on the node, but not
   the same socket, or off the node */
int MPIX_GetNestedCounts(MPI_Comm comm, MPI_Comm nodecomm, MPI_Comm socketcomm,
			 int nr, const int ranks[],
			 int *onsock, int *onnode, int *offnode)
{
    int *nranks, *sranks;
    int i;
    MPI_Group cgroup, ngroup, sgroup;

    nranks = (int *)malloc(2*nr * sizeof(int));
    if (!nranks) {
	fprintf(stderr, "Unable to allocate %d words\n", nr);
        fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    sranks = nranks + nr;

    MPI_Comm_group(comm, &cgroup);
    MPI_Comm_group(nodecomm, &ngroup);
    MPI_Comm_group(socketcomm, &sgroup);

    /* First, which are on the node */
    MPI_Group_translate_ranks(cgroup, nr, ranks, ngroup, nranks);
    MPI_Group_translate_ranks(cgroup, nr, ranks, sgroup, sranks);
    *onsock  = 0;
    *onnode  = 0;
    *offnode = 0;
    for (i=0; i<nr; i++) {
	if (nranks[i] == MPI_PROC_NULL || nranks[i] == MPI_UNDEFINED)
	    *offnode = *offnode + 1;
	else {
	    if (sranks[i] == MPI_PROC_NULL || sranks[i] == MPI_UNDEFINED) {
		*onnode = *onnode + 1;
	    }
	    else
		*onsock = *onsock + 1;
	}
    }

    free(nranks);
    MPI_Group_free(&cgroup);
    MPI_Group_free(&ngroup);
    MPI_Group_free(&sgroup);

    return 0;
}

/* The routines rankToCoords and coordsToRank can be compiled for either
   row-major or column-major, with the defines USE_ROW_MAJOR (the default)
   or USE_COLUMN_MAJOR */

PRIVATE void rankToCoords(int ndims, const int dims[], int rank, int coords[])
{
    int i, s;
    s = dims[0];
    for (i=1; i<ndims; i++)
	s *= dims[i];
#ifdef USE_COLUMN_MAJOR
    for (i=ndims-1; i>=0; i--) {
	s = s / dims[i];
	coords[i] = rank / s;
	rank = rank - coords[i] * s;
    }
#else
    for (i=0; i<ndims; i++) {
	s = s / dims[i];
	coords[i] = rank / s;
	rank = rank - coords[i] * s;
    }
#endif
}

PRIVATE void coordsToRank(int ndims, const int dims[], const int coords[],
			  int *rank)
{
    int i, s, r;
    s = dims[0];
    for (i=1; i<ndims; i++)
	s *= dims[i];
#ifdef USE_COLUMN_MAJOR
    r = coords[ndims-1];
    for (i=ndims-2; i>=0; i--) {
        r = r * dims[i] + coords[i];
    }
#else
    r = coords[0];
    for (i=1; i<ndims; i++) {
	r = r * dims[i] + coords[i];
    }
#endif
    *rank = r;
}

PRIVATE int findInInterval(int val, int offset, int high)
{
    while (offset + val < 0)
	val += high;
    while (offset + val >= high)
	val -= high;
    return val;
}

PRIVATE void rankShift(int ndims, const int dims[], const int coords[],
		       const int periodic[], int rank,
		       int direction, int disp, int *rsource, int *rdest)
{
    int rfrom, rto;
    int offset, i;

    offset = 1;
#ifdef USE_COLUMN_MAJOR
    for (i=0; i<direction; i++) offset *= dims[i];
#else
    for (i=ndims-1; i>direction; i--) offset *= dims[i];
#endif

    rfrom = -disp;
    rto   = disp;
    if (periodic[direction]) {
	/* Allow disp to be negative, so must make both in the range
	   [0,dims[direction]-1] */
	rfrom = findInInterval(rfrom, coords[direction], dims[direction]);
	rto   = findInInterval(rto,   coords[direction], dims[direction]);
	rfrom = rank + rfrom * offset;
	rto   = rank + rto * offset;
    }
    else {
	if (rfrom + coords[direction] < 0 ||
	    rfrom + coords[direction] >= dims[direction])
	    rfrom =  MPI_PROC_NULL;
	else
	    rfrom = rank + rfrom * offset;
	if (rto + coords[direction] < 0 ||
	    rto + coords[direction] >= dims[direction])
	    rto = MPI_PROC_NULL;
	else
	    rto = rank + rto * offset;
    }
    *rsource = rfrom;
    *rdest   = rto;
}

/* Given "inval" from each process in comm, order them and return the rank
   of this process in the list */

static int chiprankcompare(const void *a, const void *b);
static int chiprankcompare(const void *a, const void *b)
{
    int *aa = (int *)a;
    int *bb = (int *)b;
    int r;

    /* Compare the chip # first */
    r = *bb - *aa;
    if (r != 0) return r;

    /* Same chip, so return the core # on that chip */
    return bb[1] - aa[1];
}
static void GetChipCoreNumber(MPI_Comm comm, int chipnum, int corenum,
			      int *chiprank, int *corerank)
{
    int csize, i, chipr, corer, chipid, coreid;
    int *allvals, invals[2];
    MPI_Comm_size(comm, &csize);
    allvals = (int *)malloc(2*csize*sizeof(int));
    invals[0] = chipnum;
    invals[1] = corenum;
    MPI_Allgather(invals, 2, MPI_INT, allvals, 2, MPI_INT, comm);
    qsort(allvals, csize, 2*sizeof(int), chiprankcompare);

    chipr = -1;
    corer = -1;
    chipid = -10;
    coreid = -10;
    for (i=0; i<csize; i++) {
	if (allvals[2*i] != chipid) {
	    chipr ++;
	    chipid = allvals[2*i];
	    corer = 0;
	}
        else corer++;

	if (chipid == chipnum && allvals[2*i+1] == corenum) {
	    *chiprank = chipr;
	    *corerank = corer;
	    break;
	}
    }

    /* Sanity check */
    if (i == csize) {
	fprintf(stderr, "Could not find input value %d,%d in allvals!\n",
		chipnum, corenum);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    free(allvals);
}

/* ------------------------------------------------------------------------ */

/* Temporary routine to provide for testing.  This could be implemented
   using CVARs or removed */
void MPIX_Nodecart_cvar_set(const char *name, int value)
{
    if (strcmp(name, "ppn") == 0)
	cvar_nodecart_ppn     = value;
    else if (strcmp(name, "debug") == 0)
	cvar_nodecart_verbose = value;
    else {
	fprintf(stderr, "Unrecognized cvar %s\n", name); fflush(stderr);
    }
}

/* ----------------------------------------------------------------------- */
/* Internal functions */

/* Attribute functions */
static int nodeinfoDelFn(MPI_Comm comm, int keyval, void *attr, void *estate)
{
    nodeinfo_t *ninfo = (nodeinfo_t *)attr;

    if (!ninfo) return MPI_ERR_OTHER;

    if (ninfo->nodecomm != MPI_COMM_NULL)
	MPI_Comm_free(&ninfo->nodecomm);
    if (ninfo->leadercomm != MPI_COMM_NULL)
	MPI_Comm_free(&ninfo->leadercomm);
    free(ninfo);

    return 0;
}

/* Attribute functions */
static int nodecartDelFn(MPI_Comm comm, int keyval, void *attr, void *estate)
{
    nodecart_t *ninfo = (nodecart_t *)attr;

    if (!ninfo) return MPI_ERR_OTHER;

    /* For now, nodecomm and leadercomm are just references.  An MPI
       implementation should increment and decrement their reference
       counts.  It is not necessary to Dup them. (or is it?) */
    free(ninfo);

    return 0;
}

static int procnodeCopyFn(MPI_Comm comm, int keyval, void *estate,
			  void *attr_in, void *attr_out, int *flag)
{
    procnode_t *pnew, *pold =(procnode_t*)attr_in;
    /* Make a copy of the values in pold */
    pnew = (procnode_t *)malloc(sizeof(procnode_t));
    pnew->nchips            = pold->nchips;
    pnew->chiprank          = pold->chiprank;
    pnew->rankinchip        = pold->rankinchip;
    *(procnode_t**)attr_out = pnew;
    *flag = 1;
    return MPI_SUCCESS;
}
static int procnodeDelFn(MPI_Comm comm, int keyval, void *attr, void *estate)
{
    free(attr);
    return MPI_SUCCESS;
}
static int nodetopoCopyFn(MPI_Comm comm, int keyval, void *estate,
			  void *attr_in, void *attr_out, int *flag)
{
    nodetopo_t *tnew, *told =(nodetopo_t*)attr_in;
    /* Make a copy of the values in told */
    tnew             = (nodetopo_t *)malloc(sizeof(nodetopo_t));
    tnew->nnodes     = told->nnodes;
    tnew->noderank   = told->noderank;
    tnew->rankinnode = told->rankinnode;
    tnew->nodesize   = told->nodesize;
    *(nodetopo_t**)attr_out = tnew;
    *flag = 1;
    return MPI_SUCCESS;
}
static int nodetopoDelFn(MPI_Comm comm, int keyval, void *attr, void *estate)
{
    free(attr);
    return MPI_SUCCESS;
}

/* This routine is used to factor the number of processes on a single node,
   and so assumes that the number of processes on a single node is no
   greater than about 1k.  Adding a few more factors would push that up
   to over 16K. */
/* Worst case need at most 6 unique factors to reach over 2310, so the
   factor arrays need not be large */
#define MAX_FACTORS 7
static int smallprimes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
	37, 41, 43, 47, 53, 59, 61, -1 };
PRIVATE int factor(int n, int *nf_ptr, int factors[], int powers[])
{
    int i, nf = 0;

    /* Lots of optimizations possible.  E.g., once smallprimes > sqrt(n),
       you are done, since you've already removed the small
     */
    for (i=0; smallprimes[i] > 0; i++) {
	int f = smallprimes[i];
	if ((n % f) == 0) {
	    factors[nf] = f;
	    powers[nf]  = 1;
	    n = n / f;
	    while ((n % f) == 0) {
		powers[nf]++;
		n = n / f;
	    }
	    nf++;
	    if (n == 1) break;
	}
    }
    if (n > 1) {
	/* Last value is prime > sqrt(n) */
	factors[nf] = n;
	powers[nf]  = 1;
	nf++;
    }
    *nf_ptr = nf;
    return 0;
}
/* This is a heuristic but simple algorithm:
   1) Factor the node size (nsize)
   2) Starting with the *largest* prime factor, find the largest remaing
      dim that it divides exactly, and apply it there
   3) Continue until all factors are used
  This does not guarantee a "best" decomposition (which also depends on your
  definition of best", but is easy to implement and in common cases, does
  produce an optimal decomposition.
*/
PRIVATE int decompDims(int ndims, const int dims[], int nsize,
		       int intradims[], int interdims[])
{
    int f[MAX_FACTORS], p[MAX_FACTORS], nf;
    int i, k, fac, remain, ridx, remaindims[10];

    factor(nsize, &nf, f, p);

    for (i=0; i<ndims; i++) {
	intradims[i]  = 1;
	remaindims[i] = dims[i];
    }

    nf = nf-1;
    while (nf >= 0) {
	fac = f[nf];
	p[nf]--;
	if (p[nf] == 0) nf--;
	remain = nsize;
	ridx   = -1;
	for (k=0; k<ndims; k++) {
	    if ( (remaindims[k] % fac) == 0) {
		if (intradims[k] < remain) {
		    remain = intradims[k];
		    ridx   = k;
		}
	    }
	}
	if (ridx < 0) {
	    return -1;
	}
	remaindims[ridx] = remaindims[ridx]/fac;
	intradims[ridx]  *= fac;
    }
    for (k=0; k<ndims; k++)
	interdims[k] = dims[k] / intradims[k];
    return 0;
}

static void decompDimsError(int ndims, const int dims[], nodecart_t *ninfo)
{
    int i;
    /* PANIC: Error in result. */
	fprintf(stderr, "Unable to get node decomposition! ndims=%d\n",
		ndims);
	fprintf(stderr, "\tdims = [");
	for (i=0; i<ndims; i++)
	    fprintf(stderr, "%d%c", dims[i], (i == ndims-1) ? ']':',');
	fprintf(stderr, "\n\tintradims = [");
	for (i=0; i<ndims; i++)
	    fprintf(stderr, "%d%c", ninfo->intradims[i],
		    (i == ndims-1) ? ']':',');
	fprintf(stderr, "\n\tinterdims = [");
	for (i=0; i<ndims; i++)
	    fprintf(stderr, "%d%c", ninfo->interdims[i],
		    (i == ndims-1) ? ']':',');
	fprintf(stderr, "\n");
	fflush(stderr);
}

/*@
  MPIX_NodeSubset - Create a subset of a given communicator with a specified
  number of nodes

Input Parameters:
+ incomm - Communicator from which to form the subset
. leadercomm - Communicator of node leaders for 'incomm'
. nodecomm - Communicator of processes on this node
- nnodes - Number of nodes to retain in the subset

Output Parameter:
. outcommPtr - Pointer to the created subset communicator.

Notes:
  Given an incomm, with leadercomm of the nodes and nodecomm for the
  processes on the same node as the calling process, create a subset
  communicator with nnodes nodes.

@*/
int MPIX_NodeSubset(MPI_Comm incomm, MPI_Comm leadercomm, MPI_Comm nodecomm,
		    int nnodes, MPI_Comm *outcommPtr)
{
    int color, inrank;

    /* 1. Set the split color to 1 for nodes to include on leadercomm
       Also check that nnodes is in range */
    if (leadercomm != MPI_COMM_NULL) {
	int lrank, lsize;
	MPI_Comm_rank(leadercomm, &lrank);
	MPI_Comm_size(leadercomm, &lsize);
	if (lsize < nnodes) {
	    if (lrank == 0) {
		fprintf(stderr,
			"NodeSubset: number of nodes %d > available %d\n",
			nnodes, lsize);
                fflush(stderr);
	    }
	    *outcommPtr = MPI_COMM_NULL;
	    return MPI_ERR_OTHER;
	}
	color = MPI_UNDEFINED;
	if (lrank < nnodes) color = 1;
    }
    MPI_Bcast(&color, 1, MPI_INT, 0, nodecomm);

    /* 2. Split on color */
    MPI_Comm_rank(incomm, &inrank);
    MPI_Comm_split(incomm, color, inrank, outcommPtr);

    return MPI_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* Print the multi-level process decomposition sizes to fp */
int MPIX_PrintNodeDecomp(FILE *fp, MPI_Comm comm)
{
    /* Only need the process decomposition the first time */
    int i, j, k, nlevel, ndim, dims[16];
    MPIX_Nodecart_get_decomp(comm, &nlevel, &ndim, dims);
    k = 0;
    for (i=0; i<nlevel; i++) {
	fprintf(fp,"[");
	for (j=0; j<ndim; j++)
	    fprintf(fp, "%d%c", dims[k++], (j<ndim-1)?'x':']');
    }
    if (i < nlevel - 1) fprintf(fp, " x ");
    fprintf(fp, "\n");
    fflush(fp);
    return MPI_SUCCESS;
}

int MPIX_PrintNodeCommCounts(FILE *fp, MPI_Comm comm, int nr, const int ranks[],
			     MPI_Comm nodecomm)
{
    int ncnt[2], mincnt[2], maxcnt[2], sumcnt[2];
    int nsum[2], nmincnt[2], nmaxcnt[2], nsumcnt[2];
    int wrank, csize;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(comm, &csize);

    MPIX_GetOffNodeCounts(comm, nr, ranks,
	                  nodecomm, &ncnt[0]/*onnode*/, &ncnt[1]/*offnode*/);
    if (cvar_nodecart_verbose) {
	/* In verbose case, add this information to the output file */
	fprintf(fp, "%d: (nr=%d) on=%d off=%d\n", wrank, nr, ncnt[0], ncnt[1]);
	fflush(fp);
    }

    /* Compute the min, max, and average by process */
    MPI_Allreduce(ncnt, mincnt, 2, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(ncnt, maxcnt, 2, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(ncnt, sumcnt, 2, MPI_INT, MPI_SUM, comm);
    /* Compute the number of on and off node for the *node* */
    MPI_Allreduce(ncnt, nsum, 2, MPI_INT, MPI_SUM, nodecomm);
    /* Compute the min/max/average by node */
    MPI_Allreduce(nsum, nmaxcnt, 2, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(nsum, nmincnt, 2, MPI_INT, MPI_MIN, comm);
    /* nsumcnt should be sumcnt *ppn */
    MPI_Allreduce(nsum, nsumcnt, 2, MPI_INT, MPI_SUM, comm);
    if (wrank == 0) {
	int resultlen;
	char label[MPI_MAX_OBJECT_NAME];
	MPI_Comm_get_name(comm, label, &resultlen);
	/* Note label is the null string if no name was set */
	fprintf(fp,
		"Communicator %s (min/max/avg) on-node %d/%d/%.2f",
		label, mincnt[0], maxcnt[0],
		((double)sumcnt[0])/csize);
	fprintf(fp, " off-node %d/%d/%.2f\n",
		mincnt[1], maxcnt[1], ((double)sumcnt[1])/csize);
	fprintf(fp,
		"\tby node: on-node %d/%d/%.2f off-node %d/%d/%.2f\n",
		nmincnt[0], nmaxcnt[0], ((double)nsumcnt[0])/csize,
		nmincnt[1], nmaxcnt[1], ((double)nsumcnt[1])/csize);
	fflush(fp);
    }

    return MPI_SUCCESS;
}

int MPIX_PrintNodeCommCounts_X(FILE *fp,
			       MPI_Comm comm, int nr, const int ranks[],
			       MPI_Comm nodecomm, MPI_Comm socketcomm)
{
    int ncnt[3], mincnt[3], maxcnt[3], sumcnt[3];
    int nsum[3], nmincnt[3], nmaxcnt[3], nsumcnt[3];
    int wrank, csize, i;
    int onnode[2], offnode[2], onsock;
    MPI_Comm subcomms[2];
    int resultlen;
    char label[MPI_MAX_OBJECT_NAME];

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(comm, &csize);

    MPI_Comm_get_name(comm, label, &resultlen);
    /* Note label is the null string if no name was set */

    subcomms[0] = nodecomm;
    subcomms[1] = socketcomm;
    MPIX_GetOffNodeCounts_X(comm, nr, ranks, 2, subcomms,
			    onnode, offnode);
    MPIX_GetNestedCounts(comm, nodecomm, socketcomm, nr, ranks,
			 &onsock, &onnode[0], &offnode[0]);

    if (cvar_nodecart_verbose) {
	/* In verbose case, add this information to the output file */
	for (i=0; i<2; i++) {
	    fprintf(fp, "%d: level %d (nr=%d) on=%d off=%d\n",
		    wrank, i, nr, onnode[i], offnode[i]);
	}
	fflush(fp);
    }

    for (i=0; i<2; i++) {
	/* Compute the min, max, and average by process */
	ncnt[0] = onnode[i];
	ncnt[1] = offnode[i];
	MPI_Allreduce(ncnt, mincnt, 2, MPI_INT, MPI_MIN, comm);
	MPI_Allreduce(ncnt, maxcnt, 2, MPI_INT, MPI_MAX, comm);
	MPI_Allreduce(ncnt, sumcnt, 2, MPI_INT, MPI_SUM, comm);
	/* Compute the number of on and off node for the *node* */
	MPI_Allreduce(ncnt, nsum, 2, MPI_INT, MPI_SUM, nodecomm);
	/* Compute the min/max/average by node */
	MPI_Allreduce(nsum, nmaxcnt, 2, MPI_INT, MPI_MAX, comm);
	MPI_Allreduce(nsum, nmincnt, 2, MPI_INT, MPI_MIN, comm);
	/* nsumcnt should be sumcnt *ppn */
	MPI_Allreduce(nsum, nsumcnt, 2, MPI_INT, MPI_SUM, comm);
	if (wrank == 0) {
	    const char *sublabel;

	    if (i == 0) sublabel = "node";
	    else if (i == 1) sublabel = "socket";
	    else sublabel = "unknown";
	    fprintf(fp,
		    "Communicator %s (min/max/avg) on-%s %d/%d/%.2f",
		    label, sublabel, mincnt[0], maxcnt[0],
		    ((double)sumcnt[0])/csize);
	    fprintf(fp, " off-%s %d/%d/%.2f\n", sublabel,
		    mincnt[1], maxcnt[1], ((double)sumcnt[1])/csize);
	    fprintf(fp,
		    "\tby node: on-%s %d/%d/%.2f off-%s %d/%d/%.2f\n",
		    sublabel,
		    nmincnt[0], nmaxcnt[0], ((double)nsumcnt[0])/csize,
		    sublabel,
		    nmincnt[1], nmaxcnt[1], ((double)nsumcnt[1])/csize);
	    fflush(fp);
	}
    }

    ncnt[0] = onsock;
    ncnt[1] = onnode[0];
    ncnt[2] = offnode[0];

    MPI_Allreduce(ncnt, mincnt, 3, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(ncnt, maxcnt, 3, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(ncnt, sumcnt, 3, MPI_INT, MPI_SUM, comm);
    /* Compute the number of on and off node for the *node* */
    MPI_Allreduce(ncnt, nsum, 3, MPI_INT, MPI_SUM, nodecomm);
    /* Compute the min/max/average by node */
    MPI_Allreduce(nsum, nmaxcnt, 3, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(nsum, nmincnt, 3, MPI_INT, MPI_MIN, comm);
    /* nsumcnt should be sumcnt *ppn */
    MPI_Allreduce(nsum, nsumcnt, 3, MPI_INT, MPI_SUM, comm);
    if (wrank == 0) {
	fprintf(fp,
		"Communicator %s (min/max/avg) on-sock %d/%d/%.2f\n",
		label, mincnt[0], maxcnt[0],
		((double)sumcnt[0])/csize);
	fprintf(fp,
		"Communicator %s (min/max/avg) on-node(not sock) %d/%d/%.2f\n",
		label, mincnt[1], maxcnt[1],
		((double)sumcnt[1])/csize);
	fprintf(fp,
		"Communicator %s (min/max/avg) off-node %d/%d/%.2f\n",
		label, mincnt[2], maxcnt[2],
		((double)sumcnt[2])/csize);
	fprintf(fp,
		"\tby node: on-sock %d/%d/%.2f on-node %d/%d/%.2f off-node %d/%d/%.2f\n",
		nmincnt[0], nmaxcnt[0], ((double)nsumcnt[0])/csize,
		nmincnt[1], nmaxcnt[1], ((double)nsumcnt[1])/csize,
		nmincnt[2], nmaxcnt[2], ((double)nsumcnt[2])/csize);
	fflush(fp);
    }

    return MPI_SUCCESS;
}

/*
 * Temporary routine to output information about the process location
 * in the topology
 */
int MPIX_GetNodeTopoInfo(MPI_Comm comm, int *nnodes, int *noderank,
			 int *nodesize, int *rankonnode,
			 int *nchips, int *chiprank, int *rankinchip)
{
    int        flag = 0;
    nodetopo_t *nninfo=0;
    procnode_t *pinfo=0;

    if (procnodeKeyval != MPI_KEYVAL_INVALID) {
	MPI_Comm_get_attr(comm, procnodeKeyval, &pinfo, &flag);
	if (flag && pinfo) {
	    *nchips     = pinfo->nchips;
	    *chiprank   = pinfo->chiprank;
	    *rankinchip = pinfo->rankinchip;
	}
    }
    if (nodetopoKeyval != MPI_KEYVAL_INVALID) {
	MPI_Comm_get_attr(comm, nodetopoKeyval, &nninfo, &flag);
	if (!flag || !nninfo) {
	    int rank;
	    MPI_Comm_rank(comm, &rank);
	    if (rank == 0) {
		fprintf(stderr, "No nodetopo information on communicator\n");
		fflush(stderr);
	    }
	}
	else {
	    *nnodes     = nninfo->nnodes;
	    *noderank   = nninfo->noderank;
	    *nodesize   = nninfo->nodesize;
	    *rankonnode = nninfo->rankinnode;
	}
    }

    return MPI_SUCCESS;
}

/* ----------------------------------------------------------------------- */
#ifdef DEBUG_DECOMP
/* This is used as an in-file test program to test the basic operation
   of the nodecart routines */
typedef struct { int ndims, nsize, dims[MAX_DIM]; } test_t;

test_t tests[] = { {2, 2, 4, 3},
		   {2, 2, 3, 4},
		   {2, 8, 16, 16},
		   {3, 2, 4, 3, 5},
		   {3, 3, 4, 3, 5},
		   {3, 6, 4, 3, 5},
		   {3, 16, 8, 8, 4},
		   {3, 24, 25, 25, 24},
		   {3, 24, 2*2*5, 5*5, 3*2*5},
		   {3, 24, 2*2*5, 3*5, 2*5*5},
		   {0, 0, 0} };

int main(int argc, char **argv)
{
    int i, j, rc, intradims[MAX_DIM], interdims[MAX_DIM];

    MPI_Init(0,0);

    for (i=0; tests[i].ndims > 0; i++) {
	for (j=0; j<tests[j].ndims; j++) {
	    interdims[j] = 0;
	    interdims[j] = 0;
	}
	rc = decompDims(tests[i].ndims, tests[i].dims, tests[i].nsize,
		       intradims, interdims);
	if (rc) {
	    fprintf(stderr, "decompDims for test %d returned %d\n", i, rc);
	}
	printf("nsize = %d\n", tests[i].nsize);
	if (tests[i].ndims == 2) {
	    printf("dims: [%d,%d], interdims [%d,%d], intradims [%d,%d]\n",
		   tests[i].dims[0], tests[i].dims[1],
		   interdims[0], interdims[1],
		   intradims[0], intradims[1]);
	}
	else if (tests[i].ndims == 3) {
	    printf("dims: [%d,%d,%d], interdims [%d,%d,%d], intradims [%d,%d,%d]\n",
		   tests[i].dims[0], tests[i].dims[1], tests[i].dims[2],
		   interdims[0], interdims[1], interdims[2],
		   intradims[0], intradims[1], intradims[2]);
	}
	else {
	    printf("dims size of %d not supported!\n", tests[i].ndims);
	}
	fflush(stdout);
    }

    MPI_Finalize();

    return 0;
}
#endif /* DEBUG_DECOMP */
