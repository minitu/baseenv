/* */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "topoinfo.h"
/*
 * This file contains routines to compute distances and volumes of communication
 * at each level of the topology.
 */

/*
 * Gather data about the partners in communication.  This is a collective
 * call over the communicator
 */
/*@
  topodistInit - Initialize the routines to provide data about communication
 topology

Input Parameters:
+ comm - MPI Communicator
. nsend - Number of processes that this process sends to
. sendrank - Array of ranks of processes that this process sends to
. nrecv - Number of processes that this process receives from
. recvrank - Array of ranks of processes that this process receives from
- ti - Pointer to valid (initialized) 'topoinfo_t'

Output Parameters:
. td - Pointer to valid and initialized 'topodist_t' structure

Return Value:

  @*/
int topodistInit(MPI_Comm comm, int nsend, int sendrank[],
		 int nrecv, int recvrank[],
		 topoinfo_t *ti, topodist_t **td)
{
    topoarray_t *targetarray;
    topodist_t  *dinfo;
    MPI_Request *rq;
    int         i;
    MPI_Datatype arrtype;

    dinfo           = (topodist_t *)malloc(sizeof(topodist_t));

    /* Get my information as an array */
    topoToArray(ti, dinfo->myarray.mycoords, dinfo->myarray.maxcoords,
		&dinfo->myarray.nlevels,
		&dinfo->myarray.nodeidx, MAX_TOPO_TOTAL_DIM);
    MPI_Comm_rank(MPI_COMM_WORLD, &dinfo->myarray.rank);

    MPI_Type_contiguous(3+2*MAX_TOPO_TOTAL_DIM, MPI_INT, &arrtype);
    MPI_Type_commit(&arrtype);

    /* Allocate space for partner topoinfo */
    targetarray = (topoarray_t *)malloc(nsend * sizeof(topoarray_t));
    dinfo->ntargets = nsend;
    dinfo->tarray   = targetarray;
    /* Irecv to send ranks */
    rq = (MPI_Request *)malloc((nrecv + nsend)* sizeof(MPI_Request));
    for (i=0; i<nsend; i++) {
	MPI_Irecv(&targetarray[i], 1, arrtype, sendrank[i], 0, comm, &rq[i]);
    }
    /* Isend to recv ranks */
    for (i=0; i<nrecv; i++) {
	MPI_Isend(&dinfo->myarray, 1, arrtype, recvrank[i], 0, comm,
		  &rq[i+nsend]);
    }
    /* Waitall */
    MPI_Waitall(nrecv+nsend, rq, MPI_STATUSES_IGNORE);
    free(rq);

    MPI_Type_free(&arrtype);
    *td = dinfo;

    return 0;
}

/*@
  topoMeshHopDistance - Computing the distance from the calling process to the
  process with the given rank

Input Parameters:
+ dt - Pointer to a valid topodist_t structure
- rank - Rank of the process to determine the distance to.  See notes.

Output Parameter:
. dist - distance in hops to the process with rank 'rank'

Return value:
 '0' on success, non-zero otherwise.  '-1' if the topology does not include
 a mesh or torus.

Notes:
 'rank' must belong to one of the send ranks used to initialize the topodist_t
 'dt'.  Thus, this routine is `not` collective and may be used by any process
 at any time once 'dt' has been initialized.

 The distance is computing in `hops`: this is the 1-norm or Manhatten distance
 between the calling process and the process with 'rank' in a mesh or torus
 topology.
  @*/
int topoMeshHopDistance(topodist_t *dt, int rank, int *dist)
{
    int i, j, d;

    /* Note: This code currently *assumes* a mesh for the interconnect
       past the node. */
    for (i=0; i<dt->ntargets; i++) {
	if (dt->tarray[i].rank == rank) {
	    if (dt->myarray.nodeidx != dt->tarray[i].nodeidx) {
		return -1;
	    }
	    if (dt->myarray.nlevels != dt->tarray[i].nlevels) {
		return -1;
	    }
	    d = 0;
	    for (j=dt->myarray.nodeidx; j<dt->myarray.nlevels; j++) {
		int d1 = dt->myarray.mycoords[j] - dt->tarray[i].mycoords[j];
		if (d1 < 0) d1 = -d1;
		d += d1;
	    }
	    *dist = d;
	    return 0;
	}
    }

    return -1;
}

void topodistFree(topodist_t *dt)
{
    free(dt->tarray);
    free(dt);
}

/*
 * Notes: I want to be able to compute and print out:
 *   For each NODE, the amount of data and number of targets sent, received
 *   For each CHIP, ditto
 *
 * So, need a routine that returns these values for each level of the topology
 *
 *   For each rank pair, compute the distance (hops).  Hops may be weighted.
 */

/*@
  topodistNodeCommInfo - Provide some information about off-node communication

  Input Parameters:
+ td - Topodist information returned by topodistInit
- comm - Communicator of processes to use.  This routine is collective over
 this communicator.  The most common choice of 'comm' is 'MPI_COMM_WORLD'

  Output Parameters:
+ partnerOnNode - Number of processes that this process communicates with
  that are on the same node
. partnerOffNode - As above, but for partners that are off node
. nodecomm_p     - A communicator containing only the processes that share
  the node with the calling process.
- totalOffNode - The total number of processes on this node that are
  communicating to one or more processes on a different node

  @*/
int topodistNodeCommInfo(topodist_t *td, MPI_Comm comm,
			 int *partnerOnNode, int *partnerOffNode,
			 MPI_Comm *nodecomm_p, int *totalOffNode)
{
    int nOffNode,        /* Number of partners off of the node for each
			    process */
	nOnNode,         /* Number of partners on the node for each process */
	i, j,
	sameNode;        /* */
    int color, rank, nrank, anyoff;
    MPI_Comm nodecomm;

    nOffNode = 0;
    nOnNode  = 0;

    /* For each partner, compute on/offnode information */
    for (i=0; i<td->ntargets; i++) {
	/* are the topoinfo values comparable? */
	sameNode = 0;
	if (td->tarray[i].nodeidx == td->myarray.nodeidx &&
	    td->tarray[i].nlevels == td->myarray.nlevels) {
	    /* Check the topo indices above the node */
	    sameNode = 1;
	    for (j=td->myarray.nodeidx; j<td->myarray.nlevels; j++) {
		if (td->tarray[i].mycoords[j] != td->myarray.mycoords[j])
		    sameNode = 0;
	    }
	}
	if (sameNode) {
	    nOnNode++;
	}
	else {
	    nOffNode++;
	}
    }

    *partnerOnNode  = nOnNode;
    *partnerOffNode = nOffNode;

    /* Find all of information by node, as defined above */
    MPI_Comm_rank(comm, &rank);
    color = 0;
    /* Check for the case where no interconnect network is known */
    if (td->myarray.nodeidx < td->myarray.nlevels)
	color = td->myarray.mycoords[td->myarray.nodeidx];
    for (i=td->myarray.nodeidx+1; i<td->myarray.nlevels; i++)
	color = td->myarray.mycoords[i] + td->myarray.maxcoords[i-1]*color;
    MPI_Comm_split(comm, color, rank, &nodecomm);

    anyoff = (nOffNode > 0);
    MPI_Allreduce(&anyoff, totalOffNode, 1, MPI_INT, MPI_SUM, nodecomm);
    MPI_Comm_rank(nodecomm, &nrank);
    *nodecomm_p = nodecomm;

    return 0;
}
