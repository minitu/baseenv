/*
 * Nonblocking version with datatypes and RMA
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi_rma.c in the Advanced MPI tutorial */
/* Update these to use Win_allocate */
void stencil_rma_init(procInfo *pi, stencilInfo *si, stencilTime *st)
{
}
void stencil_rma_free(procInfo *pi, stencilInfo *si, stencilTime *st)
{
}

double stencil_rma(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    MPI_Comm comm=pi->comm;
    MPI_Win  win;
    double t, tt;
    MPI_Aint size = (bx+2)*(by+2);
    double *commbuf;
    double *restrict sbufnorth, *restrict sbufsouth, *restrict sbufeast,
	*restrict sbufwest;   /* send buffers */
    double *restrict rbufnorth, *restrict rbufsouth, *restrict rbufeast,
	*restrict rbufwest;   /* receive buffers index into window */
    double heat; /* total heat in system */
#ifdef DEBUG_HALO
    static int dbgcnt = 0;
#endif

    /* Create the windows from preallocated memory */

    t = MPI_Wtime();

    /* Requires that anew = aold + size */
    MPI_Win_allocate((2*bx+2*by)*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &commbuf, &win);

    st->commInit = MPI_Wtime() - t;
    if (aold + size != anew) {
	fprintf(stderr, "aold + size + %p; anew = %p\n", aold + size, anew);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    sbufnorth = (double*)calloc(1,bx*sizeof(double)); /* send buffers */
    sbufsouth = (double*)calloc(1,bx*sizeof(double));
    sbufeast  = (double*)calloc(1,by*sizeof(double));
    sbufwest  = (double*)calloc(1,by*sizeof(double));

    rbufnorth = (double*)&commbuf[0]; /* receive buffers index into window */
    rbufsouth = (double*)&commbuf[bx];
    rbufeast  = (double*)&commbuf[2*bx];
    rbufwest  = (double*)&commbuf[(2*bx+by)];

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);
    tt = MPI_Wtime();
    for(iter=0; iter<si->niters; ++iter) {

	/* exchange data with neighbors */
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NRD-FBefore");
	}
#endif
	MPI_Win_fence(0, win);

	t = MPI_Wtime();
        if (pi->north != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) sbufnorth[i] = aold[ind(i+1,1)]; /* pack loop - last valid region */
        if (pi->south != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) sbufsouth[i] = aold[ind(i+1,by)]; /* pack loop */
        if (pi->east != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) sbufeast[i] = aold[ind(bx,i+1)]; /* pack loop */
        if (pi->west != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) sbufwest[i] = aold[ind(1,i+1)]; /* pack loop */

	st->commPack += MPI_Wtime() - t;
	t = MPI_Wtime();

/* To avoid an apparent bug in Cray MPI, we don't call MPI_Put if the target
   process is MPI_PROC_NULL (i.e., no process) */
	if (pi->north != MPI_PROC_NULL)
	    MPI_Put(sbufnorth, bx, MPI_DOUBLE, pi->north, bx, bx, MPI_DOUBLE, win);
	if (pi->south != MPI_PROC_NULL)
	    MPI_Put(sbufsouth, bx, MPI_DOUBLE, pi->south, 0, bx, MPI_DOUBLE, win);
	if (pi->east != MPI_PROC_NULL)
	    MPI_Put(sbufeast, by, MPI_DOUBLE, pi->east, 2*bx+by, by, MPI_DOUBLE, win);
	if (pi->west != MPI_PROC_NULL)
	    MPI_Put(sbufwest, by, MPI_DOUBLE, pi->west, 2*bx, by, MPI_DOUBLE, win);

	st->commStart += MPI_Wtime() - t;
	t = MPI_Wtime();

	MPI_Win_fence(0, win);

        st->commComplete += MPI_Wtime() - t;
	t = MPI_Wtime();

#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
            printSoln(si, aold, "NRD-FAfter");
        }
	dbgcnt++;
#endif

        if (pi->north != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) aold[ind(i+1,0)] = rbufnorth[i]; /* unpack loop - into ghost cells */
        if (pi->south != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbufsouth[i]; /* unpack loop */
        if (pi->east != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbufeast[i]; /* unpack loop */
        if (pi->west != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) aold[ind(0,i+1)] = rbufwest[i]; /* unpack loop */
	st->commUnpack += MPI_Wtime() - t;

	/* update grid points */
	heat = 0.0;
	t = MPI_Wtime();
	for(j=1; j<by+1; ++j) {
	    for(i=1; i<bx+1; ++i) {
		anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		heat += anew[ind(i,j)];
#endif
	    }
	}
	st->compUpdate += MPI_Wtime() - t;

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

    MPI_Win_free(&win);

    return heat;
}
