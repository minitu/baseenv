/*
 * Simple nonblocking version with no overlap, using persistent routines
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi.c in the Advanced MPI tutorial */
double stencil_pnb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    /* allocate communication buffers */
    double *restrict sbufnorth = (double*)calloc(1,bx*sizeof(double)); /* send buffers */
    double *restrict sbufsouth = (double*)calloc(1,bx*sizeof(double));
    double *restrict sbufeast = (double*)calloc(1,by*sizeof(double));
    double *restrict sbufwest = (double*)calloc(1,by*sizeof(double));
    double *restrict rbufnorth = (double*)calloc(1,bx*sizeof(double)); /* receive buffers */
    double *restrict rbufsouth = (double*)calloc(1,bx*sizeof(double));
    double *restrict rbufeast = (double*)calloc(1,by*sizeof(double));
    double *restrict rbufwest = (double*)calloc(1,by*sizeof(double));
    double heat=0; /* total heat in system */
    double t, tt;    /* Temp used to get times */
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    MPI_Comm comm=pi->comm;
    MPI_Request reqs[8];
#ifdef DEBUG_HALO
    static int dbgcnt = 0;
#endif

    t = MPI_Wtime();
    MPI_Send_init(sbufnorth, bx, MPI_DOUBLE, pi->north, 9, comm, &reqs[0]);
    MPI_Send_init(sbufsouth, bx, MPI_DOUBLE, pi->south, 9, comm, &reqs[1]);
    MPI_Send_init(sbufeast, by, MPI_DOUBLE, pi->east, 9, comm, &reqs[2]);
    MPI_Send_init(sbufwest, by, MPI_DOUBLE, pi->west, 9, comm, &reqs[3]);
    MPI_Recv_init(rbufnorth, bx, MPI_DOUBLE, pi->north, 9, comm, &reqs[4]);
    MPI_Recv_init(rbufsouth, bx, MPI_DOUBLE, pi->south, 9, comm, &reqs[5]);
    MPI_Recv_init(rbufeast, by, MPI_DOUBLE, pi->east, 9, comm, &reqs[6]);
    MPI_Recv_init(rbufwest, by, MPI_DOUBLE, pi->west, 9, comm, &reqs[7]);
    st->commInit = MPI_Wtime() - t;

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);

    tt = MPI_Wtime();
    for(iter=0; iter<si->niters; ++iter) {

	/* exchange data with neighbors */
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NP-UBefore");
	}
#endif
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
        MPI_Startall(8,reqs);
	st->commStart += MPI_Wtime() - t;

	t = MPI_Wtime();
	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	st->commComplete += MPI_Wtime() - t;

	t = MPI_Wtime();
        if (pi->north != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) aold[ind(i+1,0)] = rbufnorth[i]; /* unpack loop - into ghost cells */
        if (pi->south != MPI_PROC_NULL)
	    for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbufsouth[i]; /* unpack loop */
        if (pi->east != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbufeast[i]; /* unpack loop */
        if (pi->west != MPI_PROC_NULL)
	    for(i=0; i<by; ++i) aold[ind(0,i+1)] = rbufwest[i]; /* unpack loop */
	st->commUnpack += MPI_Wtime() - t;
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NP-UAfter");
        }
	dbgcnt++;
#endif

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

#if !COMPUTE_HEAT_EACH_ITERATION
    heat = 0.0;
    for(j=1; j<by+1; ++j) {
	for(i=1; i<bx+1; ++i) {
	    heat += aold[ind(i,j)];
	}
    }
#endif
    MPI_Allreduce(MPI_IN_PLACE, &heat, 1, MPI_DOUBLE, MPI_SUM, comm);

    for (i=0; i<8; i++) MPI_Request_free(&reqs[i]);
    free(sbufnorth); free(sbufsouth); free(sbufeast); free(sbufwest);
    free(rbufnorth); free(rbufsouth); free(rbufeast); free(rbufwest);

    return heat;
}
