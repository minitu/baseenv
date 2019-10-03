#include <stdio.h>
#include "mpi.h"
#include "nodecart.h"

/* This is a test program for the generation of a 2-level decomposition */

int decompDims(int ndims, const int dims[], int nsize,
		       int intradims[], int interdims[]);
void rankToCoords(int ndims, const int dims[], int rank, int coords[]);
void coordsToRank(int ndims, const int dims[], const int coords[],
			  int *rank);

#define MAX_DIM 3

int main(int argc, char **argv)
{
    int i, f, dims[MAX_DIM], intradims[MAX_DIM], interdims[MAX_DIM],
	ndims, nodesize;
    int socketdims[MAX_DIM], nsockets;
    int intrasocket[MAX_DIM], intersocket[MAX_DIM];
    int nodecoords[MAX_DIM], byrank[MAX_DIM];
    int totsize, rank, socksize;

    MPI_Init(0,0);

    ndims = 2;
    for (i=0; i<ndims; i++) {
	dims[i]        = 0;
	intradims[i]   = 0;
	interdims[i]   = 0;
        socketdims[i]  = 0;
	intrasocket[i] = 0;
	intersocket[i] = 0;
    }

    nsockets = 2;
    nodesize = 8;  /* 8 processes on node, 4 per socket */
    socksize = nodesize / nsockets;

#if 0
    int nsize, nodedims[MAX_DIM];
    /* Try using dims_create to make a nice breakdown of the sockets.
       But note that this has the same problems as before - there may
       be user constraints on the breakdown */
    MPI_Dims_create(nsockets, ndims, socketdims);
    nodedims[0] = 0;
    nodedims[1] = 0;

    dims[0] = 4;
    dims[1] = 3;
    nsize   = 2;

    /* Decompose the node across the socket and between socket */
    f = decompDims(ndims, nodedims, nsockets, socketdims, intradims);
    /* Decompose the */
    f += decompDims(ndims, dims, nsize, intradims, interdims);
    printf("perim = %d, nsize = %d\n", f, nsize);
    printf("dims: [%d,%d], interdims [%d,%d], intradims [%d,%d]\n",
	   dims[0], dims[1], interdims[0], interdims[1],
	   intradims[0], intradims[1]); fflush(stdout);
#endif
    totsize = 8*nodesize;
    /* Create an initial dims */
    MPI_Dims_create(totsize, ndims, dims);
    /* First, decompose to the node */
    f = decompDims(ndims, dims, nodesize, intradims, interdims);
    /* Now, decompose the intradims (dims on the node) into sockets */
    f += decompDims(ndims, intradims, nodesize/nsockets, intrasocket,
		    intersocket);
    /* Location hierarchy is now:
       interdims
       intersocket
       intrasocket
    */
    printf("nodesize = %d, nsockets = %d\n", nodesize, nsockets);
    printf("Decomp is\ndims =\t[%d,%d]\n", dims[0], dims[1]);
    printf("interdims =\t[%d,%d]\n", interdims[0], interdims[1]);
    printf("intersock =\t[%d,%d]\n", intersocket[0], intersocket[1]);
    printf("intrasock =\t[%d,%d]\n", intrasocket[0], intrasocket[1]);

    /* Also need to do ranktocoords and coordstorank and check the mapping */
    for (rank=0; rank<totsize; rank++) {
	/* convert each rank into the global coordinates, and the
	   coordinates in each level of the hierarchy (node, socket, core) */
	int noderank = rank / nodesize;
	int rankinnode = rank - noderank * nodesize;
	int srank    = rankinnode / socksize;
	int crank    = rank % socksize;
	int intercoords[MAX_DIM], socketcoords[MAX_DIM], intracoords[MAX_DIM];
	rankToCoords(ndims, interdims,   noderank, intercoords);

	rankToCoords(ndims, intersocket, srank,    socketcoords);
	rankToCoords(ndims, intrasocket, crank,    intracoords);
	printf("rank %d: [%d,%d] - [%d,%d] - [%d,%d]\n", rank,
	       intercoords[0], intercoords[1],
	       socketcoords[0], socketcoords[1],
	       intracoords[0], intracoords[1]);

	/* Compute the nodecoords from the socket and intrasocket coords */
	for (i=0; i<ndims; i++) {
	    nodecoords[i] = socketcoords[i]*intrasocket[i] + intracoords[i];
	}
	rankToCoords(ndims, intradims, rankinnode, byrank);
	printf("rank %d: coords in node [%d,%d]; from node rank, [%d,%d]\n",
	       rank, nodecoords[0], nodecoords[1],
	       byrank[0], byrank[1]);

	/* Now, redo this with a different set of socket and core #s */
	srank    = rankinnode % nsockets;
	crank    = rankinnode / nsockets;
	rankToCoords(ndims, intersocket, srank,    socketcoords);
	rankToCoords(ndims, intrasocket, crank,    intracoords);
	printf("rank(rr) %d: [%d,%d] - [%d,%d] - [%d,%d]\n", rank,
	       intercoords[0], intercoords[1],
	       socketcoords[0], socketcoords[1],
	       intracoords[0], intracoords[1]);

	/* Compute the nodecoords from the socket and intrasocket coords */
	for (i=0; i<ndims; i++) {
	    nodecoords[i] = socketcoords[i]*intrasocket[i] + intracoords[i];
	}
	rankToCoords(ndims, intradims, rankinnode, byrank);
	printf("rank(rr) %d: coords in node [%d,%d]; from node rank, [%d,%d]\n",
	       rank, nodecoords[0], nodecoords[1],
	       byrank[0], byrank[1]);
    }
    /* Now need to go backwards.  Given the coords on the node (we have socket
       and core), compute the rank for the process in the desired global grid
       (but using an arbitrary node location for the process) */
    /* Lets assume that the node# is rank/nodesize, but that the process
       are assigned in round robin to the sockets.  We want our updated mapping
       to take processes on the same socket as neighbors as much as possible */

#if 0
    /* The idea is to take the intradims, find compatible socketdims, and
       coredims is created from that */
    /* */
    for (i=0; i<ndims; i++) {
	socketdims[i] = 0;
	coredims[i]   = 0;
    }
    MPI_Dims_create(ndims, nsockets, socketdims);
    MPI_Dims_create(ndims, ncores, coredims);

    for (rank=0; rank<totsize; rank++) {
	/* Enumeration of this rank across node, socket, core */
	int noderank = rank / nodesize;
	int srank = (rank - noderank * nodesize) % nsockets;
	int crank = (rank - noderank * nodesize) / nsockets;
	int socketcoords[MAX_DIM], corecoords[MAX_DIM];

	/* What we'd like to do is compute an alternate rank so that there is
	   less inter-socket communication */

	rankToCoords(2, socketdims, srank, socketcoords);
	rankToCoords(2, coredims,   crank, corecoords);
    }
#endif

    return 0;
}
