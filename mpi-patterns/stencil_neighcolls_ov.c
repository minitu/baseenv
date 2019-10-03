/*
 * Use nonblocking neighbor collectives
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi_carttopo_neighcolls.c in the Advanced MPI tutorial.
   Note that this requires a cartesian topology */
double stencil_neighcolls_ov(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    double t, tt;

    double heat; /* total heat in system */

    double *restrict sbuf = (double*)calloc(1,2*bx*sizeof(double)+2*by*sizeof(double)); /* send buffer (west, east, north, south) */
    double *restrict rbuf = (double*)calloc(1,2*bx*sizeof(double)+2*by*sizeof(double)); /* receive buffer (w, e, n, s) */
    int counts[4]; /* = {by, by, bx, bx}; */
    int displs[4]; /* = {0, by, 2*by, 2*by+bx}; */

    counts[0] = by; counts[1] = by; counts[2] = bx; counts[3] = bx;
    displs[0] = 0;  displs[1] = by; displs[2] = 2*by; displs[3] = 2*by+bx;

    /* Start all processes at roughly the same time */
    MPI_Barrier(pi->comm);

    tt = MPI_Wtime();
    for(iter=0; iter<si->niters; ++iter) {
	MPI_Request req;
	/* exchange data with neighbors */
	t = MPI_Wtime();
	for(i=0; i<by; ++i) sbuf[i] = aold[ind(1,i+1)]; /* pack west */
	for(i=0; i<by; ++i) sbuf[by+i] = aold[ind(bx,i+1)]; /* pack east */
	for(i=0; i<bx; ++i) sbuf[2*by+i] = aold[ind(i+1,1)]; /* pack north */
	for(i=0; i<bx; ++i) sbuf[2*by+bx+i] = aold[ind(i+1,by)]; /* pack south */
	st->commPack += MPI_Wtime() - t;
	t = MPI_Wtime();
	MPI_Ineighbor_alltoallv(sbuf, counts, displs, MPI_DOUBLE, rbuf, counts, displs, MPI_DOUBLE, pi->comm, &req);
	st->commStart += MPI_Wtime() - t;

	/* update inner grid points */
	t = MPI_Wtime();
	heat = 0.0;
	for(j=2; j<by; ++j) {
	    for(i=2; i<bx; ++i) {
		anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		heat += anew[ind(i,j)];
#endif
	    }
	}
	st->compInterior += MPI_Wtime() - t;

	t = MPI_Wtime();
	MPI_Wait(&req, MPI_STATUS_IGNORE);
	st->commComplete += MPI_Wtime() - t;

	t = MPI_Wtime();
	for(i=0; i<by; ++i) aold[ind(0,i+1)] = rbuf[i]; /* unpack loop */
	for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbuf[by+i]; /* unpack loop */
	for(i=0; i<bx; ++i) aold[ind(i+1,0)] = rbuf[2*by+i]; /* unpack loop */
	for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbuf[2*by+bx+i]; /* unpack loop */
	st->commUnpack += MPI_Wtime() - t;

	t = MPI_Wtime();
	/* update outer grid points */
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
    MPI_Barrier(pi->comm);
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
    MPI_Allreduce(MPI_IN_PLACE, &heat, 1, MPI_DOUBLE, MPI_SUM, pi->comm);

    free(sbuf); free(rbuf);

    return heat;
}
