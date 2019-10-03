#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "nodecart.h"

void printLayout(FILE *fp, MPI_Comm outcomm, MPI_Comm comm, const char *label);

int main(int argc, char *argv[])
{
    MPI_Comm nodecomm, leadercomm, ncartcomm;
    int      i, dims[2], periods[2];
    int      noderank, nnodes, wsize, wrank;
    FILE     *fp=0;
    char     *defaultOut = "mapping.txt";
    char     *outfile    = 0;

    MPI_Init(0,0);

    if (argc == 2) outfile = argv[1];

    MPIX_Nodecomm_create(MPI_COMM_WORLD, &nodecomm, &leadercomm, &noderank,
			 &nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    if (wrank == 0) {
	if (!outfile) outfile = defaultOut;
	fp = fopen(outfile, "w");
    }

    for (i=0; i<2; i++) {
	dims[i]    = 0;
	periods[i] = 0;
    }
    MPI_Dims_create(wsize, 2, dims);
    MPIX_Nodecart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &ncartcomm);
    /* Some information is created by Nodecart_create, so even for
       MPI_COMM_WORLD, we wait until we have run Nodecart_create */
    printLayout(fp, MPI_COMM_WORLD, MPI_COMM_WORLD,
		"Layout for MPI_COMM_WORLD");
    printLayout(fp, MPI_COMM_WORLD, ncartcomm,
		"Layout for ncart_create in WORLD order");
    printLayout(fp, ncartcomm, ncartcomm,
		"Layout for ncart_create in cartcomm order");

    if (wrank == 0) {
	fclose(fp);
    }

    MPI_Comm_free(&ncartcomm);

    MPI_Finalize();
    return 0;
}

/* Output information on the layout of processes in comm to the file
   fp. Add a label if label is non-NULL. The process with rank 0 in
   outcomm writes all output, and output is in rank order for outcomm.
   use MPI_COMM_WORLD for outcomm to get rank order in world, use comm
   to get rank order in the communicator comm.
   Comm must have the nodeinfo attribute set.
*/
void printLayout(FILE *fp, MPI_Comm outcomm, MPI_Comm comm, const char *label)
{
    int      vals[9];
    int      noderank, nnodes, nodesize, rankonnode, nchips, chiprank,
	rankinchip, worldrank;
    int      wrank, wsize, rc;

    MPI_Comm_rank(outcomm, &wrank);
    MPI_Comm_size(outcomm, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);

    rc = MPIX_GetNodeTopoInfo(comm, &nnodes, &noderank, &nodesize,
			      &rankonnode, &nchips, &chiprank, &rankinchip);

    MPI_Allreduce(MPI_IN_PLACE, &rc, 1, MPI_INT, MPI_MAX, outcomm);
    if (rc != 0) return;

    vals[0] = nnodes;
    vals[1] = noderank;
    vals[2] = nodesize;
    vals[3] = rankonnode;
    vals[4] = nchips;
    vals[5] = chiprank;
    vals[6] = rankinchip;
    MPI_Comm_rank(comm, &vals[7]);
    MPI_Comm_rank(MPI_COMM_WORLD, &vals[8]);

    if (wrank == 0) {
	int  i, rvals[9];
	if (label) fprintf(fp, "%s\n", label);
	fprintf(fp, "i\twrank\trank\tnnodes\tnrank\tnsize\tronnode\tnchips\trchip\tronchip\n");
	fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 0, vals[8], vals[7],
		vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6]);
	for (i=1; i<wsize; i++) {
	    MPI_Recv(rvals, 9, MPI_INT, i, 0, outcomm,
		     MPI_STATUS_IGNORE);
	    fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, rvals[8],
		    rvals[7],
		    rvals[0], rvals[1], rvals[2], rvals[3], rvals[4],
		    rvals[5], rvals[6]);
	}
	fflush(fp);
    }
    else {
	MPI_Send(vals, 9, MPI_INT, 0, 0, outcomm);
    }
}
