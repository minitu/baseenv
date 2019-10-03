/*
 * Simplified Fortran interface for the MPIX_Nodecart routines
 */
#include "../baseenv.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpi.h"
#include "nodecart.h"

#define MPI_CART_CREATE FC_FUNC_(mpi_cart_create,MPI_CART_CREATE)
#define MPI_CART_SHIFT FC_FUNC_(mpi_cart_shift,MPI_CART_SHIFT)
#define MPI_CART_COORDS FC_FUNC_(mpi_cart_coords,MPIC_NODECART_COORDS)
#define MPI_CART_RANK FC_FUNC_(mpi_cart_rank,MPI_CART_RANK)
#define MPI_CART_SUB FC_FUNC_(mpi_cart_sub,MPI_CART_SUB)
#define MPIX_CARTSHIM_INFO FC_FUNC_(mpix_cartshim_info,MPIX_CARTSHIM_INFO)

/* Track use of these routines. */
static int cartCreateCalls = 0;
static int cartSubCalls = 0;

#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif

/* Define the prototoypes for strict compiles */
void MPI_CART_CREATE(MPI_Fint *incomm, MPI_Fint *ndims, MPI_Fint dims[],
		     MPI_Fint periods[], MPI_Fint *reorder,
		     MPI_Fint *comm_cart, MPI_Fint *ierr);
void MPI_CART_SHIFT(MPI_Fint *comm, MPI_Fint *direction, MPI_Fint *disp,
		    MPI_Fint *rank_sources, MPI_Fint *rank_dest,
		    MPI_Fint *ierr);
void MPI_CART_COORDS(MPI_Fint *comm, MPI_Fint *rank, MPI_Fint *maxdims,
		     MPI_Fint coords[], MPI_Fint *ierr);
void MPI_CART_RANK(MPI_Fint *comm, MPI_Fint coords[], MPI_Fint *rank,
		   MPI_Fint *ierr);
void MPI_CART_SUB(MPI_Fint *comm, MPI_Fint remain[], MPI_Fint *subcomm,
		  MPI_Fint *ierr);
int MPIX_CARTSHIM_INFO(MPI_Fint *createcalls, MPI_Fint *subcalls);
/* End of prototypes */

/* Beginning of shims */

void MPI_CART_CREATE(MPI_Fint *incomm, MPI_Fint *ndims, MPI_Fint dims[],
		     MPI_Fint periods[], MPI_Fint *reorder,
		     MPI_Fint *comm_cart, MPI_Fint *ierr)
{
    int      err;
    MPI_Comm c_comm_cart;
    /* FIXME: Assume sizeof(MPI_Fint) == sizeof(int) for the arrays
       dims and periods */
    err = MPIX_Nodecart_create(MPI_Comm_f2c(*incomm), *ndims, dims,
			       periods, *reorder, &c_comm_cart);
    if (err) { *ierr = (MPI_Fint)err; return; }
    else *ierr = MPI_SUCCESS;
    *comm_cart = MPI_Comm_c2f(c_comm_cart);
}
void MPI_CART_SHIFT(MPI_Fint *comm, MPI_Fint *direction, MPI_Fint *disp,
		    MPI_Fint *rank_sources, MPI_Fint *rank_dest,
		    MPI_Fint *ierr)
{
    int c_source, c_dest, err;
    err = MPIX_Nodecart_shift(MPI_Comm_f2c(*comm), *direction, *disp,
			      &c_source, &c_dest);
    *rank_sources = c_source;
    *rank_dest    = c_dest;
    *ierr         = err;
}
void MPI_CART_COORDS(MPI_Fint *comm, MPI_Fint *rank, MPI_Fint *maxdims,
		     MPI_Fint coords[], MPI_Fint *ierr)
{
    int err;
    err   = MPIX_Nodecart_coords(MPI_Comm_f2c(*comm), *rank, *maxdims, coords);
    *ierr = err;
}
void MPI_CART_RANK(MPI_Fint *comm, MPI_Fint coords[], MPI_Fint *rank,
		   MPI_Fint *ierr)
{
    int err;
    err = MPIX_Nodecart_rank(MPI_Comm_f2c(*comm), coords, rank);
    *ierr = err;
}
void MPI_CART_SUB(MPI_Fint *comm, MPI_Fint remain[], MPI_Fint *subcomm,
		  MPI_Fint *ierr)
{
    int      err;
    MPI_Comm subcomm_p;
    err = MPIX_Nodecart_sub(MPI_Comm_f2c(*comm), remain, &subcomm_p);
    if (err == MPI_SUCCESS) {
	*subcomm = MPI_Comm_c2f(subcomm_p);
    }
    *ierr = err;
}

/* This routine makes it easier to check that the shim routines were
   invoked.
*/
int MPIX_CARTSHIM_INFO(MPI_Fint *createcalls, MPI_Fint *subcalls)
{
    *createcalls = cartCreateCalls;
    *subcalls    = cartSubCalls;
    return MPI_SUCCESS;
}
