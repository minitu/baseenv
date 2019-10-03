/* -*- Mode: C; c-basic-offset:4 ; -*- */
/* See 6.3.2 in the BG/Q Redbook */
#include "mpix.h"

#define MAX_TORUS_DIMS 5
static int haveMPIXTopo = 0;
static MPIX_Hardware_t TopoMeshCoords;
static int meshNdims;
static int myCoords[MAX_TORUS_DIMS+1];

static int topoiGetNodeTopoInfoBGQ(topoinfo_t *topo);

static int topoiGetNodeTopoInfoBGQ(topoinfo_t *topo)
{
    int wrank, rc;
    int loc[5];
    topoentry_t *e=0;

    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
    /* This horrible code is because the name of this routine changed
       for no reason sometime after the documentation was released. */
#ifdef HAVE_MPIX_HARDWARE
    if (!haveMPIXTopo) {
        MPIX_Hardware( &TopoMeshCoords );
        haveMPIXTopo = 1;
    }
#else
    if (!haveMPIXTopo) {
        MPIX_Init_hw( &TopoMeshCoords );
        haveMPIXTopo = 1;
    }
#endif
    MPIX_Torus_ndims( &meshNdims );

    e = topoiAllocEntry(topo);
    e->dim = meshNdims;
    e->topoType = TOPO_TORUS;

    /* The torus coords are also available in TopoMeshCoords */
    MPIX_Rank2torus( wrank, e->coords.coords );

    /* There is also Torus2Rank( int coords[], int *rank ) to convert
       in the other direction */
    MPI_Allreduce(e->coords.coords, e->mincoords.coords, meshNdims, MPI_INT,
	MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce(e->coords.coords, e->maxcoords.coords, meshNdims, MPI_INT,
	MPI_MAX, MPI_COMM_WORLD );

    return 0;
}

/* Indicate that this file also provides the node info */
#define PROVIDED_GETNODEINFO 1

/* Use information from the MPIX_Init_hw to get info on the core */
static int topoiGetNodeInfoBGQ(int isMultithreaded, topoinfo_t *topo)
{
    int cpunum;
    topoentry_t *e;

#ifdef HAVE_MPIX_HARDWARE
    if (!haveMPIXTopo) {
        MPIX_Hardware( &TopoMeshCoords );
        haveMPIXTopo = 1;
    }
#else
    if (!haveMPIXTopo) {
        MPIX_Init_hw( &TopoMeshCoords );
        haveMPIXTopo = 1;
    }
#endif

    cpunum = TopoMeshCoords.coreID;
    if (verbose) printf("CPU number = %d\n", cpunum);
    if (cpunum >= 0) {
	e = topoiAllocEntry(topo);
	if (!e) return -1;
	e->dim = 1;
	e->maxcoords.coords[0] = TopoMeshCoords.ppn;
	e->topoType = TOPO_CORE;
	e->coords.coords[0] = cpunum;
	e->coords.coords[1] = -1;
    }
    return 0;
}

