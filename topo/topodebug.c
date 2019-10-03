/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include <stdio.h>
#include "mpi.h"
#include "topoconf.h"
#include "topoinfo.h"
#include "topoimpl.h"

static int debugTopoKind=1;

int topoiSetDummyTopo(int kind)
{
    debugTopoKind = kind;
    return 0;
}

int topoiGetNodeTopoInfoDEBUG(topoinfo_t *topo)
{
    int         wsize, wrank;
    topoentry_t *e=0;

    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (debugTopoKind == 1) {
	int ppn;
	/* Create a topology similar to a Cray XE6.
	   The mapping is this:
	   px=py=pz=2 (for now)
	   if wsize % 8 != 0,
	   procs/node = wsize/(px*py*pz)
	   procs/chip = procs/node/2
	   Simple linear assignment of ranks to processes.
	 */
	ppn = wsize / 8;
	/* printf( "ppn = %d, wsize = %d\n", ppn, wsize); */
	if ( (wsize % 8) != 0 || ppn == 0) {
	    ppn = 1;
	}
	wrank   = wrank / ppn;     /* wrank is now wrt the nodes */

	e = topoiAllocEntry(topo);
	e->dim = 4;
	e->topoType = TOPO_TORUS;
	/* Fake the mesh coordinates for this rank */
	e->coords.coords[3] = wrank % 2;  wrank /= 2;
	e->coords.coords[0] = wrank % 2;  wrank /= 2;
	e->coords.coords[1] = wrank % 2;  wrank /= 2;
	e->coords.coords[2] = wrank;

	/* Blue Waters is a 24 x 24 x 24 torus */
	e->maxtopocoords.coords[0] = 24;
	e->maxtopocoords.coords[1] = 24;
	e->maxtopocoords.coords[2] = 24;
	e->maxtopocoords.coords[3] = 2;  /* 2 nodes per NIC */

	/* Setup the min and max coords */
	MPI_Allreduce(e->coords.coords, e->mincoords.coords, 4, MPI_INT,
		      MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(e->coords.coords, e->maxcoords.coords, 4, MPI_INT,
		      MPI_MAX, MPI_COMM_WORLD);
    }
    else {
	/* Unknown debug topology kind */
	return 1;
    }
    return 0;
}

int topoiGetNodeInfoDEBUG(int isMultithreaded, topoinfo_t *topo)
{
    int wsize, wrank, corenum=0, ppchip, chip;
    topoentry_t *e=0;

    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (debugTopoKind == 1) {
	int ppn;
	/* Create a topology similar to a Cray XE6.
	   The mapping is this:
	   px=py=pz=2 (for now)
	   if wsize % 8 != 0,
	   procs/node = wsize/(px*py*pz)
	   procs/chip = procs/node/2
	   Simple linear assignment of ranks to processes.
	 */
	ppn = wsize / 8;
	/* printf( "ppn = %d, wsize = %d\n", ppn, wsize); */
	if ( (wsize % 8) != 0 || ppn == 0) {
	    ppn = 1;
	}
	/* If an even number of processes per node, alternate them among
	   the 2 chips */
	if ((ppn % 2) == 0) {
	    chip    = wrank % 2;
	    corenum = (wrank % ppn) / 2;
	    ppchip  = ppn / 2;
	}
	else {
	    chip    = -1;
	    corenum = wrank % ppn;
	    ppchip  = ppn;
	}
	wrank   = wrank / ppn;     /* wrank is now wrt the nodes */
	if (corenum >= 0) {
	    e = topoiAllocEntry(topo);
	    if (!e) return 1;
	    e->dim = 1;
	    e->maxcoords.coords[0] = ppchip;
	    e->topoType            = TOPO_CORE;
	    e->coords.coords[0]    = corenum;
	    e->coords.coords[1]    = -1;
	}
	if (chip >= 0) {
	    e = topoiAllocEntry(topo);
	    if (!e) return 1;
	    e->dim = 1;
	    e->maxcoords.coords[0] = 2; /* 2 chips per node */
	    e->topoType = TOPO_SOCKET;
	    e->coords.coords[0] = chip;
	    e->coords.coords[1] = -1;
	    corenum = (wrank % ppn) / 2;
	    ppchip  = ppn / 2;
	}
    }
    else {
	/* Unknown debug topology kind */
	return 1;
    }
    return 0;
}


#if 0
/*
 * This internal routine is used to help debug both the topo package itself
 * and programs that are using the package.  Based on the 'kind' parameter,
 * it creates a topology (rather than providing the true topology of the
 * system)
 */
int topoiSetupDummyTopo(int kind, topoinfo_t *topo)
{
    int wsize, wrank, corenum=0, ppchip, chip;
    topoentry_t *e=0;

    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (kind == 1) {
	int ppn;
	/* Create a topology similar to a Cray XE6.
	   The mapping is this:
	   px=py=pz=2 (for now)
	   if wsize % 8 != 0,
	   procs/node = wsize/(px*py*pz)
	   procs/chip = procs/node/2
	   Simple linear assignment of ranks to processes.
	 */
	ppn = wsize / 8;
	/* printf( "ppn = %d, wsize = %d\n", ppn, wsize); */
	if ( (wsize % 8) != 0 || ppn == 0) {
	    ppn = 1;
	}
	/* If an even number of processes per node, alternate them among
	   the 2 chips */
	if ((ppn % 2) == 0) {
	    chip    = wrank % 2;
	    corenum = (wrank % ppn) / 2;
	    ppchip  = ppn / 2;
	}
	else {
	    chip    = -1;
	    corenum = wrank % ppn;
	    ppchip  = ppn;
	}
	wrank   = wrank / ppn;     /* wrank is now wrt the nodes */
	if (corenum >= 0) {
	    e = topoiAllocEntry(topo);
	    if (!e) return 1;
	    e->dim = 1;
	    e->maxcoords.coords[0] = ppchip;
	    e->topoType            = TOPO_CORE;
	    e->coords.coords[0]    = corenum;
	    e->coords.coords[1]    = -1;
	}
	if (chip >= 0) {
	    e = topoiAllocEntry(topo);
	    if (!e) return 1;
	    e->dim = 1;
	    e->maxcoords.coords[0] = 2; /* 2 chips per node */
	    e->topoType = TOPO_SOCKET;
	    e->coords.coords[0] = chip;
	    e->coords.coords[1] = -1;
	    corenum = (wrank % ppn) / 2;
	    ppchip  = ppn / 2;
	}

	e = topoiAllocEntry(topo);
	e->dim = 4;
	e->topoType = TOPO_TORUS;
	/* Fake the mesh coordinates for this rank */
	e->coords.coords[3] = wrank % 2;  wrank /= 2;
	e->coords.coords[0] = wrank % 2;  wrank /= 2;
	e->coords.coords[1] = wrank % 2;  wrank /= 2;
	e->coords.coords[2] = wrank;

	/* Blue Waters is a 24 x 24 x 24 torus */
	e->maxtopocoords.coords[0] = 24;
	e->maxtopocoords.coords[1] = 24;
	e->maxtopocoords.coords[2] = 24;
	e->maxtopocoords.coords[3] = 2;  /* 2 nodes per NIC */
    }
    else {
	/* Unknown debug topology kind */
	return 1;
    }
    return 0;
}
#endif
