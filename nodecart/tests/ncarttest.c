/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "nodecart.h"

int main(int argc, char *argv[])
{
    int wrank, wsize, north, south, east, west;
    MPI_Comm cartcomm, ncartcomm;
    int dims[3], coords[3], periods[3], debug=0, i;
    int crank, ncrank, ranks[5], wranks[5];
    MPI_Comm nodecomm, leadercomm;
    MPI_Group gworld, gcart, gnodecart;
    int nnodes, noderank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Look for debugging and verbose arguments */
    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-debug") == 0) debug = 1;
	else {
	    if (wrank == 0) {
		fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD,1);
	    }
	}
    }

    /* Temp: For testing on laptop, use debug option to set the node size */
    if (debug) {
	if (wsize > 32 && (wsize % 16) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 16);
	}
	else if (wsize >= 12 && (wsize % 6) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 6);
	    MPIX_Nodecart_cvar_set("debug", 2);
	}
	else if (wsize >=8 && (wsize % 4) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 4);
	}
	else if ((wsize & 0x1) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 2);
	}
    }

    MPI_Comm_group(MPI_COMM_WORLD, &gworld);

    /* Get the node size */
    MPIX_Nodecomm_create(MPI_COMM_WORLD, &nodecomm, &leadercomm,
			 &noderank, &nnodes);

    if (wrank == 0) {
	int nsize;
	MPI_Comm_size(nodecomm, &nsize);
	printf("SMP: nodes = %d, nodesize = %d\n", nnodes, nsize);
	fflush(stdout);
    }

    /* Print the processes in the same node for each leader */
    if (leadercomm != MPI_COMM_NULL) {
	int j, nsize, *nranks, *nwranks, nodenum;
	MPI_Group ngroup;

	MPI_Comm_size(nodecomm, &nsize);
	MPI_Comm_group(nodecomm, &ngroup);
	nranks = (int *)malloc(2*nsize*sizeof(int));
	nwranks = nranks + nsize;
	for (j=0; j<nsize; j++) nranks[j] = j;

	MPI_Group_translate_ranks(ngroup, nsize, nranks, gworld, nwranks);

	MPI_Comm_rank(leadercomm, &nodenum);

	/* Simple sequentialization nodes */
	for (i=0; i<nnodes; i++) {
	    if (i == nodenum) {
		printf("%d: process on node %d:", wrank, nodenum);
		for (j=0; j<nsize; j++) printf("%d ", nwranks[j]);
		printf("\n");fflush(stdout);
	    }
	}
	free(nranks);
	MPI_Group_free(&ngroup);
    }
    {
	int numnodes, noderk, nodesize, rankonnode, nchips, chiprank,
	    rankonchip;
	int err;
	err = MPIX_GetNodeTopoInfo(MPI_COMM_WORLD, &numnodes, &noderk,
				   &nodesize, &rankonnode, &nchips,
				   &chiprank, &rankonchip);
	if (err == MPI_SUCCESS) {
	    printf("%d: Nodes=%d, noderank=%d, nodesize=%d, rinnode=%d, nchips=%d, chiprank=%d, rchip=%d\n",
		   wrank, numnodes, noderk, nodesize, rankonnode,
		   nchips, chiprank, rankonchip);
	    fflush(stdout);
	}
    }

    /* MPI Cartesian routine */
    for (i=0; i<2; i++) {
	dims[i] = 0;
	periods[i] = 0;
    }
    MPI_Dims_create(wsize, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cartcomm);
    MPI_Comm_set_name(cartcomm, "cartcomm");
    MPI_Cart_shift(cartcomm, 0, 1, &west, &east);
    MPI_Cart_shift(cartcomm, 1, 1, &north, &south);
    MPI_Comm_rank(cartcomm, &crank);
    MPI_Cart_coords(cartcomm, crank, 2, coords);
    MPI_Comm_group(cartcomm, &gcart);

    ranks[0] = west;
    ranks[1] = east;
    ranks[2] = north;
    ranks[3] = south;
    ranks[4] = crank;
    MPI_Group_translate_ranks(gcart, 5, ranks, gworld, wranks);
    if (wrank == 0) printf("Cartdecomp [%d x %d]\n", dims[0], dims[1]);
    printf("%d(%d): (%d,%d):%d:%d:%d:%d\n", wrank, crank, coords[0], coords[1],
	   west, east, north, south);
    printf("%d(%d):%d: (%d,%d):%d:%d:%d:%d\n", wrank, crank, wranks[4],
	   coords[0], coords[1], wranks[0], wranks[1], wranks[2], wranks[3]);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    MPIX_PrintNodeCommCounts(stdout, cartcomm, 4, ranks, nodecomm);

    /* MPIX Nodecart routine */
    MPIX_Nodecart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &ncartcomm);
    MPI_Comm_set_name(ncartcomm, "ncartcomm");
    MPIX_Nodecart_shift(ncartcomm, 0, 1, &west, &east);
    MPIX_Nodecart_shift(ncartcomm, 1, 1, &north, &south);
    { int remain[2], sz;
	MPI_Comm ncart1;

	remain[0] = 0; remain[1] = 1;
        MPIX_Nodecart_sub(ncartcomm, remain, &ncart1);
	MPI_Comm_size(ncart1, &sz);
	printf("Sizeof sub is %d\n", sz);fflush(stdout);
        MPI_Comm_free(&ncart1);
	remain[1] = 0; remain[0] = 1;
        MPIX_Nodecart_sub(ncartcomm, remain, &ncart1);
        MPI_Comm_free(&ncart1);
    }
    MPI_Comm_rank(ncartcomm, &ncrank);
    MPIX_Nodecart_coords(ncartcomm, ncrank, 2, coords);
    MPI_Comm_group(ncartcomm, &gnodecart);

    ranks[0] = west;
    ranks[1] = east;
    ranks[2] = north;
    ranks[3] = south;
    ranks[4] = ncrank;
    MPI_Group_translate_ranks(gnodecart, 5, ranks, gworld, wranks);
    printf("%d(%d): (%d,%d):%d:%d:%d:%d\n", wrank, ncrank, coords[0], coords[1],
	   west, east, north, south);
    printf("%d(%d):%d: (%d,%d):%d:%d:%d:%d\n", wrank, ncrank, wranks[4],
	   coords[0], coords[1], wranks[0], wranks[1], wranks[2], wranks[3]);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    MPIX_PrintNodeCommCounts(stdout, ncartcomm, 4, ranks, nodecomm);

    /* Get node info and check on quality of mappings */

    /* Free communicators */
    MPI_Comm_free(&cartcomm);
    MPI_Comm_free(&ncartcomm);
    /* Free groups */
    MPI_Group_free(&gworld);
    MPI_Group_free(&gcart);
    MPI_Group_free(&gnodecart);

    MPI_Finalize();
    return 0;
}


