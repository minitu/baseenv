#include <stdio.h>
#include "mpi.h"
#include "nodecart.h"
#include "ncartshim.h"

/* This file contains interface routines that exploit the MPI profiling
   interface to allow use of the MPIX_Nodecart routines.  Because this
   shim should *only* be used when the intent is to replace the versions of
   the MPI_Cart routines from the MPI implementation */

/* Track use of these routines. */
static int cartCreateCalls = 0;
static int cartSubCalls = 0;

int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[],
		    const int periods[],
                    int reorder, MPI_Comm *comm_cart)
{
    int err;
    err = MPIX_Nodecart_create(comm_old, ndims, dims, periods, reorder,
			       comm_cart);
    cartCreateCalls++;
    return err;
}

int MPI_Cartdim_get(MPI_Comm comm, int *ndims)
{
    int err;
    err = MPIX_Nodecart_dim_get(comm, ndims);
    return err;
}

int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[],
		 int coords[])
{
    int err;
    err = MPIX_Nodecart_get(comm, maxdims, dims, periods, coords);
    return err;
}

int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank)
{
    int err;
    err = MPIX_Nodecart_rank(comm, coords, rank);
    return err;
}

int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[])
{
    int err;
    err = MPIX_Nodecart_coords(comm, rank, maxdims, coords);
    return err;
}

int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
		   int *rank_source, int *rank_dest)
{
    int err;
    err = MPIX_Nodecart_shift(comm, direction, disp, rank_source, rank_dest);
    return err;
}

int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm)
{
    int err;
    err = MPIX_Nodecart_sub(comm, remain_dims, newcomm);
    cartSubCalls++;
    return err;
}

int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[],
		 const int periods[], int *newrank)
{
    fprintf(stderr, "MPI_Cart_map not implemented for nodecart\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
    return MPI_ERR_OTHER;
}

/* This routine makes it easier to check that the shim routines were
   invoked.
*/
int MPIX_CartShim_info(int *createcalls, int *subcalls)
{
    *createcalls = cartCreateCalls;
    *subcalls    = cartSubCalls;
    return MPI_SUCCESS;
}
