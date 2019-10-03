/*
 * Simplified Fortran interface for the MPIX_Nodecart routines
 */
#include "../baseenv.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpi.h"
#include "nodecart.h"

#define MPIX_NODECOMM_CREATE FC_FUNC_(mpix_nodecomm_create,MPIX_NODECOMM_CREATE)
#define MPIX_NODECART_CREATE FC_FUNC_(mpix_nodecart_create,MPIX_NODECART_CREATE)
#define MPIX_NODECART_SHIFT FC_FUNC_(mpix_nodecart_shift,MPIX_NODECART_SHIFT)
#define MPIX_NODECART_COORDS FC_FUNC_(mpix_nodecart_coords,MPIC_NODECART_COORDS)
#define MPIX_COMM_DIMS_CREATE FC_FUNC_(mpix_comm_dims_create,MPIX_COMM_DIMS_CREATE)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void MPIX_NODECOMM_CREATE(MPI_Fint *incomm, MPI_Fint *nodecomm,
      MPI_Fint *leadercomm, MPI_Fint *noderank, MPI_Fint *nnodes,
      MPI_Fint *ierr)
{
    MPI_Comm c_nodecomm, c_leadercomm;
    int      c_noderank, c_nnodes;
    int      err;
    err = MPIX_Nodecomm_create(MPI_Comm_f2c(*incomm),
			       &c_nodecomm, &c_leadercomm, &c_noderank,
			       &c_nnodes);
    if (err) { *ierr = (MPI_Fint)err; return; }
    else *ierr = MPI_SUCCESS;
    *noderank = (MPI_Fint)c_noderank;
    *nnodes   = (MPI_Fint)c_nnodes;
    *nodecomm = (MPI_Fint)MPI_Comm_c2f(c_nodecomm);
    if (c_leadercomm != MPI_COMM_NULL)
	*leadercomm = (MPI_Fint)MPI_Comm_c2f(c_leadercomm);
    else
	*leadercomm = MPI_COMM_NULL;
}
void MPIX_NODECART_CREATE(MPI_Fint *incomm, MPI_Fint *ndims, MPI_Fint dims[],
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
void MPIX_NODECART_SHIFT(MPI_Fint *comm, MPI_Fint *direction, MPI_Fint *disp,
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
void MPIX_NODECART_COORDS(MPI_Fint *comm, MPI_Fint *rank, MPI_Fint *maxdims,
			  MPI_Fint coords[], MPI_Fint *ierr)
{
    int err;
    err   = MPIX_Nodecart_coords(MPI_Comm_f2c(*comm), *rank, *maxdims, coords);
    *ierr = err;
}
void MPIX_COMM_DIMS_CREATE(MPI_Fint *comm, MPI_Fint *nsize, MPI_Fint *ndims,
			   MPI_Fint dims[], MPI_Fint *ierr)
{
    int err;
    err = MPIX_Comm_dims_create(MPI_Comm_f2c(*comm), *nsize, *ndims, dims);
    *ierr = err;
}

