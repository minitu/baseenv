/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "topoinfo.h"
/* #include "meshtopo.h" */
/* Add prototypes for some of the topology information */
#ifdef HAVE_FINDCLIQUES_H
#include "findcliques.h"
#endif

/*
   This file contains routines to compute some of the properties of
   a mapping of MPI processes to processors.

   Operations supported:

   1. Given an array of ranks representing a set of bidirectional
   exchanges (e.g., for a halo exchange) and a "clique number" for the
   calling process, for each process compute the number of
   communication partners in the same clique.  This is a collective
   call.

   2. ToDo: version of (1), but uses a tuple representing the location of
   the calling process in a hierarchy.  This may include FP-core, chip, node,
   interconnect-node, etc.

*/

/*
 * Algorithm:
 * Create a buffer into which to receive the clique number of each partner
 * Post Irecvs, perform sends, waitall (ensures safe even if eager threshold
 * is 0 bytes).
 * Compare resulting array wih myClique.
 *
 * Question: Also determine the number of different cliques?
 */

/*
 * Given partner ranks and the cliqueNum of this process, determine the
 * cliqueNum of each of the partner processes.  Place those values in the
 * output array cliqueArray.  This is a collective routine; in addition,
 * the ranks array must represent bidirectional communication (e.g.,
 * communication used for exchanges).
 */
int TopoFindPartnerCliquenum( int nranks, const int ranks[], MPI_Comm comm,
			      int myClique, int cliqueArray[] )
{
    MPI_Request *rArray = 0;
    int i;

    rArray      = (MPI_Request *)malloc( nranks * sizeof(MPI_Request) );

    for (i=0; i<nranks; i++)
	MPI_Irecv( cliqueArray + i, 1, MPI_INT, ranks[i], 0, comm, rArray+i );
    for (i=0; i<nranks; i++)
	MPI_Send( &myClique, 1, MPI_INT, ranks[i], 0, comm );
    MPI_Waitall( nranks, rArray, MPI_STATUSES_IGNORE );

    free( rArray );

    return 0;
}

#if 0
/*
 * Given partner ranks and the topology information for of this
 * process, determine the topology information each of the partner processes.
 * Place those values in the output array cliqueArray.  This is a
 * collective routine; in addition, the ranks array must represent
 * bidirectional communication (e.g., communication used for
 * exchanges).
 */
int TopoFindPartnerTopoinfo( int nranks, const int ranks[], MPI_Comm comm,
			     TopoDesc myDesc, TopoDesc partnerInfo[] )
{
    MPI_Request *rArray = 0;
    int i;

    rArray      = (MPI_Request *)malloc( nranks * sizeof(MPI_Request) );

    for (i=0; i<nranks; i++)
	MPI_Irecv( partnerInfo[i].coords, 3, MPI_INT, ranks[i], 0,
		   comm, rArray+i );
    for (i=0; i<nranks; i++)
	MPI_Send( myDesc.coords, 3, MPI_INT, ranks[i], 0, comm );
    MPI_Waitall( nranks, rArray, MPI_STATUSES_IGNORE );

    free( rArray );

    return 0;
}

#endif

/* Determine how many partner processes are in the same clique as the
 * calling process.  This is a collective call, with the same restrictions
 */
int TopoComputeRanksSameClique(int nranks, const int ranks[], MPI_Comm comm,
			       int myClique, int *nSameClique )
{
    int *cliqueArray = 0;
    int i, match;

    cliqueArray = (int *)malloc( nranks * sizeof(int) );

    if (TopoFindPartnerCliquenum( nranks, ranks, comm,
				  myClique, cliqueArray )) {
	fprintf( stderr, "Error returned from TopoFindPartnerCliquenum\n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    match = 0;
    for (i=0; i<nranks; i++)
	if (cliqueArray[i] == myClique) match++;

    *nSameClique = match;

    return 0;
}

#if 0
/* */
int TopoMeshDistances( int nranks, const int ranks[], MPI_Comm comm,
		       int dist[] )
{
    TopoDesc myDesc, *partnerDesc;
    int i, j, d, d1;

    myDesc.kind = TOPO_MESH;
    myDesc.ndim = TopoMyMeshDim();
    partnerDesc = (TopoDesc *)malloc( nranks * sizeof(TopoDesc) );
    TopoMyMeshCoords(myDesc.ndim, &myDesc.coords[0]);

    TopoFindPartnerTopoinfo( nranks, ranks, comm, myDesc, partnerDesc );

    for (i=0; i<nranks; i++) {
	d = 0;
	for (j=0; j<myDesc.ndim; j++) {
	    d1 = myDesc.coords[j] - partnerDesc[i].coords[j];
	    if (d1 < 0) d1 = -d1;
	    d += d1;
	}
	dist[i] = d;
    }
    free( partnerDesc );

    return 0;
}

/*
 * Approximate mapping of processes to a
 */

/*
 * Temporary code to acquire information about specific systems
 */
#ifdef HAVE_CRAY_RCA
/* Based on code from A Bahtele and N. Jain, UIUC */
#include <pmi.h>
#include <rca_lib.h>

/* Check for consistency in dimension information */
#if TOPO_MAX_DIM < 4
#error TOPO_MAX_DIM must be at least 4
#endif
static rca_mesh_coord_t TopoMeshCoords;
static int myNid;
static int minCoords[TOPO_MAX_DIM], maxCoords[TOPO_MAX_DIM];

int TopoInitMesh(int *ndims)
{
    int wrank, rc;
    int loc[4];

    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
    PMI_Get_nid( wrank, &myNid );
    rc = rca_get_meshcoord( myNid, &TopoMeshCoords );
    if (rc == -1) return 1;

    /* To get the range of coords, need to min/max the values */
    loc[0] = TopoMeshCoords.mesh_x;
    loc[1] = TopoMeshCoords.mesh_y;
    loc[2] = TopoMeshCoords.mesh_z;
    /* There are two nodes per Gemini NIC.  They have consecutive
       NID values, so we just take NID mod 2. */
    loc[3] = myNID % 2;
    MPI_Allreduce( loc, minCoords, 4, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( loc, maxCoords, 4, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

    *ndims = 3; /* The Cray XE6/XT7 is a 3-D torus */
    return 0;
}
int TopoMyMeshDim(void)
{
    return 4;
}

int TopoMyMeshCoords( int ndims, int coords[] )
{
    coords[0] = TopoMeshCoords.mesh_x;
    coords[1] = TopoMeshCoords.mesh_y;
    coords[2] = TopoMeshCoords.mesh_z;

    return 0;
}

int TopoMeshContainer( int ndims, int mindims[3], int maxdims[3] )
{
    int i, maxdim;
    maxdim = (ndims > 4) ? 4 : ndims;
    for (i=0; i<maxdim; i++) {
	mindims[i] = minCoords[i];
	maxdims[i] = maxCoords[i];
    }

    return 0;
}
int TopoMeshMaxDims( int ndims, int dims[] )
{
    rca_mesh_coord_t dims;
    /* rca_get_max_dimension gives the size of the *machines* mesh, not
       the jobs mesh.  For that, we need to compute min and max of the
       coordinates. */
    rca_get_max_dimension( &dims );
    dims[0] = dims.mesh_x+1;
    dims[1] = dims.mesh_y+1;
    dims[2] = dims.mesh_z+1;
    if (ndims > 3)
	dims[3] = 2;   /* Two nodes per NIC */

    return 0;
}
#elif defined(HAVE_BGQ_MPIX)
/* See 6.3.2 in the BG/Q Redbook */
#include "mpix.h"

#define MAX_TORUS_DIMS 5
static MPIX_Hardware_t TopoMeshCoords;
static int meshNdims;
static int minCoords[MAX_TORUS_DIMS], maxCoords[MAX_TORUS_DIMS];
static int myCoords[MAX_TORUS_DIMS+1];

int TopoInitMesh(int *ndims)
{
    int wrank, rc;
    int loc[5];

    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
    /* This horrible code is because the name of this routine changed
       for no reason sometime after the documentation was released. */
#ifdef HAVE_MPIX_HARDWARE
    MPIX_Hardware( &TopoMeshCoords );
#else
    MPIX_Init_hw( &TopoMeshCoords );
#endif
    MPIX_Torus_ndims( &meshNdims );

    /* The torus coords are also available in TopoMeshCoords */
    MPIX_Rank2torus( wrank, myCoords );
    /* There is also Torus2Rank( int coords[], int *rank ) to convert
       in the other direction */
    MPI_Allreduce( myCoords, minCoords, meshNdims, MPI_INT, MPI_MIN,
		   MPI_COMM_WORLD );
    MPI_Allreduce( myCoords, maxCoords, meshNdims, MPI_INT, MPI_MAX,
		   MPI_COMM_WORLD );

    *ndims = 5;  /* BG/Q is a 5-d torus */
    return 0;
}
int TopoMyMeshDim(void)
{
    return 5;
}
int TopoMyMeshCoords( int ndims, int coords[] )
{
    int i;

    if (ndims > 5) ndims = 5;
    for (i=0; i<ndims; i++) coords[i] = myCoords[i];

    return 0;
}

int TopoMeshContainer( int ndims, int mindims[], int maxdims[] )
{
    int i;
    if (ndims > meshNdims) ndims = meshNdims;
    for (i=0; i<meshNdims; i++) {
	mindims[i] = minCoords[i];
	maxdims[i] = maxCoords[i];
    }
    return 0;
}

int TopoMeshMaxDims(int ndims, int dims[])
{
    /* Fields in the MPIX_Hardware_t include torus size in each dimension.
       Note that there are more than three dimensions */
    /*
    *xdim = dims.mesh_x+1;
    *ydim = dims.mesh_y+1;
    *zdim = dims.mesh_z+1;
    */
    return 0;
}
#else
/* No topology information is available.  Provide stubs that return
   failure.  The halo code will use failure from TopoInitMesh to
   bypass the topology-related computations. */
int TopoInitMesh(int *ndims)
{
    *ndims = 0;
    return 1;  /* Indicate no information available */
}
int TopoMyMeshDim(void)
{
    return 0;
}
int TopoMyMeshCoords(int ndims, int coords[])
{
    return 1;
}
int TopoMeshContainer(int ndims, int mindims[], int maxdims[])
{
    return 1;
}
int TopoMeshMaxDims(int ndims, int dims[])
{
    return 1;
}

#endif /* Select on Cray/BG/other */

#endif
