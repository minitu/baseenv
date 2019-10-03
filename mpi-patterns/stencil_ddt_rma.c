/*
 * Nonblocking version with datatypes and RMA
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* See stencil_mpi_ddt_rma.c in the Advanced MPI tutorial */
/* Update these to use Win_allocate */
void stencil_ddt_rma_init(procInfo *pi, stencilInfo *si, stencilTime *st)
{
}
void stencil_ddt_rma_free(procInfo *pi, stencilInfo *si, stencilTime *st)
{
}

double stencil_ddt_rma(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    double *restrict aold = si->aold, *restrict anew=si->anew, *tmp;
    double *abase;
    MPI_Comm comm=pi->comm;
    MPI_Win  win;
    MPI_Datatype north_south_type;
    MPI_Datatype east_west_type;
    double t, tt;
    MPI_Aint size = (bx+2)*(by+2);
    double heat; /* total heat in system */
#ifdef DEBUG_HALO
    static int dbgcnt = 0;
#endif

    /* create north-south datatype */
    t = MPI_Wtime();
    MPI_Type_contiguous(bx, MPI_DOUBLE, &north_south_type);
    MPI_Type_commit(&north_south_type);
    /* create east-west type */
    MPI_Type_vector(by, 1, bx+2, MPI_DOUBLE, &east_west_type);
    MPI_Type_commit(&east_west_type);
    /* Create the windows from preallocated memory */
    abase = aold;   /* used to check offset calculation */
    /* Requires that anew = aold + size */
    MPI_Win_create(aold, 2*size*sizeof(double), sizeof(double), MPI_INFO_NULL,
		   comm, &win);
    st->commInit = MPI_Wtime() - t;
    if (aold + size != anew) {
	fprintf(stderr, "aold + size + %p; anew = %p\n", aold + size, anew);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);
    tt = MPI_Wtime();
    for(iter=0; iter<si->niters; ++iter) {
	MPI_Aint offset = size*((iter)%2);
	/* exchange data with neighbors */
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NRD-FBefore");
	}
#endif
	t = MPI_Wtime();
	MPI_Win_fence(0, win);
#ifdef DEBUG_HALO
	if (abase + offset != aold) {
	    myAbort(MPI_COMM_WORLD, 1, "Incorrect calculation of offset\n");
	}
#endif

/* To avoid an apparent bug in Cray MPI, we don't call MPI_Put if the target
   process is MPI_PROC_NULL (i.e., no process) */
#ifdef USE_CONTIG_TYPE
	if (pi->north != MPI_PROC_NULL)
	    MPI_Put(&aold[ind(1,1)], 1, north_south_type, pi->north, ind(1,by+1)+offset, 1, north_south_type, win);
	if (pi->south != MPI_PROC_NULL)
	    MPI_Put(&aold[ind(1,by)], 1, north_south_type, pi->south, ind(1,0)+offset, 1, north_south_type, win);
#else
	MPI_Put(&aold[ind(1,1)], bx, MPI_DOUBLE, pi->north, ind(1,by+1)+offset, bx, MPI_DOUBLE, win);
	MPI_Put(&aold[ind(1,by)], bx, MPI_DOUBLE, pi->south, ind(1,0)+offset, bx, MPI_DOUBLE, win);
#endif
	if (pi->east != MPI_PROC_NULL)
	    MPI_Put(&aold[ind(bx,1)], 1, east_west_type, pi->east, ind(0,1)+offset, 1, east_west_type, win);
	if (pi->west != MPI_PROC_NULL)
	    MPI_Put(&aold[ind(1,1)], 1, east_west_type, pi->west, ind(bx+1,1)+offset, 1, east_west_type, win);

	st->commStart += MPI_Wtime() - t;
	t = MPI_Wtime();
	MPI_Win_fence(0, win);
	st->commComplete += MPI_Wtime() - t;
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
            printSoln(si, aold, "NRD-FAfter");
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

    MPI_Type_free(&east_west_type);
    MPI_Type_free(&north_south_type);
    MPI_Win_free(&win);

    return heat;
}
