/*
 * Simple version with nonblocking Irecv and blocking send with no overlap
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi.c in the Advanced MPI tutorial */
double stencil_bnb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int bz = si->bz;
    int i, j, k, p, iter;
    /* allocate communication buffers */
    double *restrict sbufnorth = (double*)calloc(1,bx*bz*sizeof(double)); /* send buffers */
    double *restrict sbufsouth = (double*)calloc(1,bx*bz*sizeof(double));
    double *restrict sbufeast = (double*)calloc(1,by*bz*sizeof(double));
    double *restrict sbufwest = (double*)calloc(1,by*bz*sizeof(double));
    double *restrict sbuftop  = (double*)calloc(1,bx*by*sizeof(double));
    double *restrict sbufbottom = (double*)calloc(1,bx*by*sizeof(double));
    double *restrict rbufnorth = (double*)calloc(1,bx*bz*sizeof(double)); /* receive buffers */
    double *restrict rbufsouth = (double*)calloc(1,bx*bz*sizeof(double));
    double *restrict rbufeast = (double*)calloc(1,by*bz*sizeof(double));
    double *restrict rbufwest = (double*)calloc(1,by*bz*sizeof(double));
    double *restrict rbuftop  = (double*)calloc(1,bx*by*sizeof(double));
    double *restrict rbufbottom = (double*)calloc(1,bx*by*sizeof(double));
    double heat; /* total heat in system */
    double t, tt;    /* Temp used to get times */
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    MPI_Comm comm=pi->comm;

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);

    tt = MPI_Wtime();
    for(iter=0; iter<si->niters; ++iter) {

	/* exchange data with neighbors */
	MPI_Request reqs[6];
	t = MPI_Wtime();
	/* Pack loops for each face */
        if (pi->north != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<bx; ++i) sbufnorth[p++] = aold[ind(i+1,1,k+1)];
	}
        if (pi->south != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<bx; ++i) sbufsouth[p++] = aold[ind(i+1,by,k+1)];
	}
        if (pi->east != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<by; ++i) sbufeast[p++] = aold[ind(bx,i+1,k+1)];
	}
        if (pi->west != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<by; ++i) sbufwest[p++] = aold[ind(1,i+1,k+1)];
	}
        if (pi->top != MPI_PROC_NULL) {
	    p=0;
	    for (j=0; j<bx; ++j)
		for(i=0; i<bx; ++i) sbuftop[p++] = aold[ind(i+1,j+1,1)];
	}
        if (pi->bottom != MPI_PROC_NULL) {
	    p=0;
	    for (j=0; j<by; ++j)
		for(i=0; i<bx; ++i) sbufbottom[p++] = aold[ind(i+1,j+1,bz)];
	}
	st->commPack += MPI_Wtime() - t;
	t = MPI_Wtime();
	MPI_Irecv(rbufnorth, bx*bz, MPI_DOUBLE, pi->north, 9, comm, &reqs[0]);
	MPI_Irecv(rbufsouth, bx*bz, MPI_DOUBLE, pi->south, 9, comm, &reqs[1]);
	MPI_Irecv(rbufeast, by*bz, MPI_DOUBLE, pi->east, 9, comm, &reqs[2]);
	MPI_Irecv(rbufwest, by*bz, MPI_DOUBLE, pi->west, 9, comm, &reqs[3]);
	MPI_Irecv(rbuftop, bx*by, MPI_DOUBLE, pi->top, 9, comm, &reqs[4]);
	MPI_Irecv(rbufbottom, bx*by, MPI_DOUBLE, pi->bottom, 9, comm, &reqs[5]);
	MPI_Send(sbufnorth, bx*bz, MPI_DOUBLE, pi->north, 9, comm);
	MPI_Send(sbufsouth, bx*bz, MPI_DOUBLE, pi->south, 9, comm);
	MPI_Send(sbufeast, by*bz, MPI_DOUBLE, pi->east, 9, comm);
	MPI_Send(sbufwest, by*bz, MPI_DOUBLE, pi->west, 9, comm);
	MPI_Send(sbuftop, bx*by, MPI_DOUBLE, pi->top, 9, comm);
	MPI_Send(sbufbottom, bx*by, MPI_DOUBLE, pi->bottom, 9, comm);
	st->commStart += MPI_Wtime() - t;

	t = MPI_Wtime();
	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
	st->commComplete += MPI_Wtime() - t;

	t = MPI_Wtime();
        if (pi->north != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<bx; ++i) aold[ind(i+1,0,k+1)] = sbufnorth[p++];
	}
        if (pi->south != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<bx; ++i)  aold[ind(i+1,by+1,k+1)] = sbufsouth[p++];
	}
        if (pi->east != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<by; ++i) aold[ind(bx+1,i+1,k+1)] = sbufeast[p++];
	}
        if (pi->west != MPI_PROC_NULL) {
	    p=0;
	    for (k=0; k<bz; ++k)
		for(i=0; i<by; ++i) aold[ind(0,i+1,k+1)] = sbufwest[p++];
	}
        if (pi->top != MPI_PROC_NULL) {
	    p=0;
	    for (j=0; j<bx; ++j)
		for(i=0; i<bx; ++i) aold[ind(i+1,j+1,0)] = sbuftop[p++];
	}
        if (pi->bottom != MPI_PROC_NULL) {
	    p=0;
	    for (j=0; j<by; ++j)
		for(i=0; i<bx; ++i) aold[ind(i+1,j+1,bz+1)] = sbufbottom[p++];
	}
	st->commUnpack += MPI_Wtime() - t;

	/* update grid points */
	heat = 0.0;
	t = MPI_Wtime();
	for(k=1; k<bz+1; ++k) {
	    for(j=1; j<by+1; ++j) {
		for(i=1; i<bx+1; ++i) {
		    anew[ind(i,j,k)] = aold[ind(i,j,k)]/2.0 + (aold[ind(i-1,j,k)] + aold[ind(i+1,j,k)] + aold[ind(i,j-1,k)] + aold[ind(i,j+1,k)] + aold[ind(i,j,k-1)] + aold[ind(i,j,k+1)])/6.0/2.0;
#if COMPUTE_HEAT_EACH_ITERATION
		    heat += anew[ind(i,j,k)];
#endif
		}
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
    for(k=1; k<bz+1; ++k) {
	for(j=1; j<by+1; ++j) {
	    for(i=1; i<bx+1; ++i) {
		heat += aold[ind(i,j,k)];
	    }
	}
    }
#endif
    MPI_Allreduce(MPI_IN_PLACE, &heat, 1, MPI_DOUBLE, MPI_SUM, comm);

    free(sbufnorth); free(sbufsouth); free(sbufeast); free(sbufwest);
    free(sbuftop);   free(sbufbottom);
    free(rbufnorth); free(rbufsouth); free(rbufeast); free(rbufwest);
    free(rbuftop);   free(rbufbottom);

    return heat;
}
