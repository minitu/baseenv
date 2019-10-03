/*
 * Nonblocking version with datatypes and communication/computation overlap
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi_ddt_overlap.c in the Advanced MPI tutorial */
double stencil_ddt_ov(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    MPI_Comm comm=pi->comm;
    double heat; /* total heat in system */
    MPI_Datatype north_south_type;
    MPI_Datatype east_west_type;
    double t, tt;

    /* create north-south datatype */
    t = MPI_Wtime();
    MPI_Type_contiguous(bx, MPI_DOUBLE, &north_south_type);
    MPI_Type_commit(&north_south_type);
    /* create east-west type */
    MPI_Type_vector(by,1,bx+2,MPI_DOUBLE, &east_west_type);
    MPI_Type_commit(&east_west_type);
    st->commInit = MPI_Wtime() - t;

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);
    tt = MPI_Wtime();

    for(iter=0; iter<si->niters; ++iter) {

	/* exchange data with neighbors */
	MPI_Request reqs[8];
	t = MPI_Wtime();
	MPI_Isend(&aold[ind(1,1)] /* north */, 1, north_south_type, pi->north, 9, comm, &reqs[0]);
	MPI_Isend(&aold[ind(1,by)] /* south */, 1, north_south_type, pi->south, 9, comm, &reqs[1]);
	MPI_Isend(&aold[ind(bx,1)] /* east */, 1, east_west_type, pi->east, 9, comm, &reqs[2]);
	MPI_Isend(&aold[ind(1,1)] /* west */, 1, east_west_type, pi->west, 9, comm, &reqs[3]);
	MPI_Irecv(&aold[ind(1,0)] /* north */, 1, north_south_type, pi->north, 9, comm, &reqs[4]);
	MPI_Irecv(&aold[ind(1,by+1)] /* south */, 1, north_south_type, pi->south, 9, comm, &reqs[5]);
	MPI_Irecv(&aold[ind(bx+1,1)] /* west */, 1, east_west_type, pi->east, 9, comm, &reqs[6]);
	MPI_Irecv(&aold[ind(0,1)] /* east */, 1, east_west_type, pi->west, 9, comm, &reqs[7]);
	st->commStart += MPI_Wtime() - t;

	/* update inner grid points */
	heat = 0.0;
	t = MPI_Wtime();
	for(j=2; j<by; ++j) {
	    for(i=2; i<bx; ++i) {
		anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		heat += anew[ind(i,j)];
#endif
	    }
	}
	st->compInterior += MPI_Wtime() - t;

	/* wait for communication to complete */
	t = MPI_Wtime();
	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	st->commComplete += MPI_Wtime() - t;

	/* update outer grid points */
	t = MPI_Wtime();
	for(j=1; j < by+1; j+=by-1) {
	    for(i=2; i<bx; ++i) { /* north, south -- two elements less per row (first and last) to avoid "double computation" in next loop! */
		anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		heat += anew[ind(i,j)];
#endif
	    }
	}

	/* update outer grid points */
	for(j=1; j < by+1; ++j) {
	    for(i=1; i<bx+1; i+=bx-1) { /* east, west -- full columns */
		anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		heat += anew[ind(i,j)];
#endif
	    }
	}
	st->compBndy += MPI_Wtime() - t;

	/* swap arrays */
	tmp=anew; anew=aold; aold=tmp;
    }
    /* Ensure all processes are done */
    MPI_Barrier(comm);
    st->total = MPI_Wtime() - tt;

    /* get final heat in the system */
#if !COMPUTE_HEAT_EACH_ITERATION
    heat = 0.0;
    for(j=1; j<by+1; ++j) {
	for(i=1; i<bx+1; ++i) {
	    heat += aold[ind(i,j)];
	}
    }
#endif
    MPI_Allreduce(MPI_IN_PLACE, &heat, 1, MPI_DOUBLE, MPI_SUM, comm);

return heat;
}
