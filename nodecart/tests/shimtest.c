#include <stdio.h>
#include "mpi.h"
#include "ncartshim.h"

/* This is a test of the nodecart shim file.  If */
int main(int argc, char **argv)
{
    int wrank, wsize, dims[2], periods[2], i;
    int west, east, north, south, ndim=2, crank, coords[2];
    int createcalls, subcalls;
    MPI_Comm cartcomm;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    for (i=0; i<ndim; i++) {
	dims[i]    = 0;
	periods[i] = 0;
    }
    MPI_Dims_create(wsize, ndim, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 1, &cartcomm);
    MPI_Cart_shift(cartcomm, 0, 1, &west, &east);
    MPI_Cart_shift(cartcomm, 1, 1, &north, &south);
    MPI_Comm_rank(cartcomm, &crank);
    MPI_Cart_coords(cartcomm, crank, ndim, coords);

    /* This is a special routine provided by the shim library to make it
       easier to confirm that the shim library has been correctly applied
       and that the routine interception has occured */
    MPIX_CartShim_info(&createcalls, &subcalls);
    if (wrank == 0) {
	if (createcalls != 1 && subcalls != 0) {
	    printf("Test FAILED: createcalls = %d (should be 1) and subcalls = %d (should be 0)\n",
		   createcalls, subcalls);
	}
	else {
	    printf("Test PASSED\n");
	}
	fflush(stdout);
    }

    MPI_Comm_free(&cartcomm);
    MPI_Finalize();
    return 0;
}
