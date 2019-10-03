/*
 * Nonblocking version with datatypes and RMA
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

/* We use file-local variables for other data needed to implement the
   communication approach */
static MPI_Win win;
static MPI_Comm shmcomm;
static double *restrict mem;
static size_t size;
/* Set to MPI_UNDEFINED if not accessible through shared memory */
static int northSM, southSM, eastSM, westSM;

/* See stencil_mpi_shmem.c in the Advanced MPI tutorial.  That version does
   not include the internode communication */
double stencil_shmem_nb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st)
{
    int bx = si->bx;
    int by = si->by;
    int i, j, iter;
    MPI_Comm comm=pi->comm;
    double heat; /* total heat in system */

    double *tmp;
    double *restrict anew= mem+size;
    double *restrict aold= mem;

    double *restrict northptr, *restrict southptr,
	*restrict eastptr, *restrict westptr;
    double *restrict northptr2, *restrict southptr2,
	*restrict eastptr2, *restrict westptr2;

    /* Allocate space for all edges, even if not needed for communication.
       We could be more careful here and only allocate as needed */
    double *restrict sbufnorth = (double*)calloc(1,bx*sizeof(double)); /* send buffers */
    double *restrict sbufsouth = (double*)calloc(1,bx*sizeof(double));
    double *restrict sbufeast = (double*)calloc(1,by*sizeof(double));
    double *restrict sbufwest = (double*)calloc(1,by*sizeof(double));
    double *restrict rbufnorth = (double*)calloc(1,bx*sizeof(double)); /* receive buffers */
    double *restrict rbufsouth = (double*)calloc(1,bx*sizeof(double));
    double *restrict rbufeast = (double*)calloc(1,by*sizeof(double));
    double *restrict rbufwest = (double*)calloc(1,by*sizeof(double));
    MPI_Request reqs[8];
    MPI_Aint sz;
    int dsp_unit;
    double t, tt;
#ifdef DEBUG_HALO
    static int dbgcnt = 0;
#endif

    /* Because we're using a private aold and anew, we must initialize here */
    initMeshBase(pd, bx, by, aold, anew);

    if (northSM != MPI_UNDEFINED) {
	MPI_Win_shared_query(win, northSM, &sz, &dsp_unit,
				 (void*)&northptr);
	northptr2 = northptr+size;
    }
    if (southSM != MPI_UNDEFINED) {
	MPI_Win_shared_query(win, southSM, &sz, &dsp_unit,
				 (void*)&southptr);
	southptr2 = southptr+size;
    }
    if (eastSM != MPI_UNDEFINED) {
	MPI_Win_shared_query(win, eastSM, &sz, &dsp_unit,
				 (void*)&eastptr);
	eastptr2 = eastptr+size;
    }
    if (westSM != MPI_UNDEFINED) {
	MPI_Win_shared_query(win, westSM, &sz, &dsp_unit,
				 (void*)&westptr);
	westptr2 = westptr+size;
    }

    /* Preload reqs with REQUEST_NULL; we'll only update the ones that
       require communication */
    for (i=0; i<8; i++) reqs[i] = MPI_REQUEST_NULL;

    /* Start all processes at roughly the same time */
    MPI_Barrier(comm);

    tt = MPI_Wtime();
    MPI_Win_lock_all(0, win);
    for(iter=0; iter<si->niters; ++iter) {
	MPI_Win_sync(win);
	MPI_Barrier(shmcomm);
	/* exchange data with neighbors */

#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NS-UBefore");
	}
#endif
	t = MPI_Wtime();
	if (northSM == MPI_UNDEFINED) {
	    if (pi->north != MPI_PROC_NULL) {
		for(i=0; i<bx; ++i) sbufnorth[i] = aold[ind(i+1,1)];
		MPI_Isend(sbufnorth, bx, MPI_DOUBLE, pi->north, 9, comm, &reqs[0]);
		MPI_Irecv(rbufnorth, bx, MPI_DOUBLE, pi->north, 9, comm, &reqs[4]);
	    }
	}
	else if(northSM != MPI_PROC_NULL) {
	    for(i=0; i<bx; ++i) aold[ind(i+1,0)] = northptr[ind(i+1,by)]; /* pack loop - last valid region */
	}
	if (southSM == MPI_UNDEFINED) {
	    if (pi->south != MPI_PROC_NULL) {
		for(i=0; i<bx; ++i) sbufsouth[i] = aold[ind(i+1,by)];
		MPI_Isend(sbufsouth, bx, MPI_DOUBLE, pi->south, 9, comm, &reqs[1]);
		MPI_Irecv(rbufsouth, bx, MPI_DOUBLE, pi->south, 9, comm, &reqs[5]);
	    }
	}
	else if(southSM != MPI_PROC_NULL) {
	    for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = southptr[ind(i+1,1)]; /* pack loop */
	}
	if (eastSM == MPI_UNDEFINED) {
	    if (pi->east != MPI_PROC_NULL) {
		for(i=0; i<by; ++i) sbufeast[i] = aold[ind(bx,i+1)];
		MPI_Isend(sbufeast, by, MPI_DOUBLE, pi->east, 9, comm, &reqs[2]);
		MPI_Irecv(rbufeast, by, MPI_DOUBLE, pi->east, 9, comm, &reqs[6]);
	    }
	}
	else if(eastSM != MPI_PROC_NULL) {
	    for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = eastptr[ind(1,i+1)]; /* pack loop */
	}
	if (westSM == MPI_UNDEFINED) {
	    if (pi->west != MPI_PROC_NULL) {
		for(i=0; i<by; ++i) sbufwest[i] = aold[ind(1,i+1)];
		MPI_Isend(sbufwest, by, MPI_DOUBLE, pi->west, 9, comm, &reqs[3]);
		MPI_Irecv(rbufwest, by, MPI_DOUBLE, pi->west, 9, comm, &reqs[7]);
	    }
	}
	else if(westSM != MPI_PROC_NULL) {
	    for(i=0; i<by; ++i) aold[ind(0,i+1)] = westptr[ind(bx,i+1)]; /* pack loop */
	}
	st->commStart += MPI_Wtime() - t;

	t = MPI_Wtime();
	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	st->commComplete += MPI_Wtime() - t;

	/* Unback any data received with MPI_Irecv */
	t = MPI_Wtime();
	if (northSM == MPI_UNDEFINED && pi->north != MPI_PROC_NULL) {
	    for(i=0; i<bx; ++i) aold[ind(i+1,0)] = rbufnorth[i]; /* unpack loop */
	}
	if (southSM == MPI_UNDEFINED && pi->south != MPI_PROC_NULL) {
	    for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbufsouth[i]; /* unpack loop */
	}
	if (eastSM == MPI_UNDEFINED && pi->east != MPI_PROC_NULL) {
	    for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbufeast[i]; /* unpack loop */
	}
	if (westSM == MPI_UNDEFINED && pi->west != MPI_PROC_NULL) {
	    for(i=0; i<by; ++i) aold[ind(0,i+1)] = rbufwest[i]; /* unpack loop */
	}
	st->commUnpack += MPI_Wtime() - t;
#ifdef DEBUG_HALO
	if (dbgcnt >= debugFirst && dbgcnt <= debugLast) {
	    printSoln(si, aold, "NS-UAfter");
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
	tmp=northptr; northptr=northptr2; northptr2=tmp;
	tmp=southptr; southptr=southptr2; southptr2=tmp;
	tmp=eastptr; eastptr=eastptr2; eastptr2=tmp;
	tmp=westptr; westptr=westptr2; westptr2=tmp;
    }
    MPI_Win_unlock_all(win);
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

    free(sbufnorth); free(sbufsouth); free(sbufeast); free(sbufwest);
    free(rbufnorth); free(rbufsouth); free(rbufeast); free(rbufwest);

    return heat;
}

void stencil_shmem_nb_init(procInfo *pi, stencilInfo *si, stencilTime *st)
{
    MPI_Group commgroup, shmemgroup;
    int       nbrs[4], shmnbrs[4];
    size = (si->bx+2)*(si->by+2); /* process-local grid (including halos (thus +2)) */
    MPI_Comm_split_type(pi->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
			&shmcomm);
    /* Sanity check */
    if (shmcomm == MPI_COMM_NULL) {
	fprintf(stderr, "Unable to get COMM_TYPE_SHARED\n");
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Win_allocate_shared(2*size*sizeof(double), 1, MPI_INFO_NULL,
			    shmcomm, (void *)&mem, &win);

    /* We don't use si->aold and si->anew, since that memory
       is allocated by the default routine */
    /* Determine which directions are within shared memory */
    MPI_Comm_group(pi->comm, &commgroup);
    MPI_Comm_group(shmcomm, &shmemgroup);
    nbrs[0] = pi->north;
    nbrs[1] = pi->south;
    nbrs[2] = pi->east;
    nbrs[3] = pi->west;
    /* Input ranks should all be within commgroup; shmnbrs will
     be MPI_PROC_NULL if not within shared memory group */
    MPI_Group_translate_ranks(commgroup, 4, nbrs, shmemgroup, shmnbrs);
    northSM = shmnbrs[0];
    southSM = shmnbrs[1];
    eastSM  = shmnbrs[2];
    westSM  = shmnbrs[3];
    MPI_Group_free(&commgroup);
    MPI_Group_free(&shmemgroup);
}
void stencil_shmem_nb_free(procInfo *pi, stencilInfo *si, stencilTime *st)
{
    MPI_Comm_free(&shmcomm);
    MPI_Win_free(&win);
    /* Memory allocated in Win_allocate_shared is freed with the window */
}
