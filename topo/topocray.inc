/* -*- Mode: C; c-basic-offset:4 ; -*- */
/* Based on code from A Bahtele and N. Jain, UIUC */
#include <pmi.h>
#include <rca_lib.h>

static rca_mesh_coord_t TopoMeshCoords;
static int myNid;

/* Check for consistency in dimension information */
#if MAX_TOPO_DIM < 4
#error MAX_TOPO_DIM must be at least 4 in topoinfo.h
#endif

static int topoiGetNodeTopoInfoCRAY(topoinfo_t *topo);

static int topoiGetNodeTopoInfoCRAY(topoinfo_t *topo)
{
    int wrank, rc, i;
    int loc[4];
    rca_mesh_coord_t dims;
    topoentry_t *e=0;

    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
    PMI_Get_nid( wrank, &myNid );
    rc = rca_get_meshcoord( myNid, &TopoMeshCoords );
    if (rc == -1) return -1;

    e = topoiAllocEntry(topo);
    e->dim = 4;
    e->topoType = TOPO_TORUS;

    /* To get the range of coords, need to min/max the values */
    loc[0] = TopoMeshCoords.mesh_x;
    loc[1] = TopoMeshCoords.mesh_y;
    loc[2] = TopoMeshCoords.mesh_z;
    /* There are two nodes per Gemini NIC.  They have consecutive
       NID values, so we just take NID mod 2. */
    loc[3] = myNid % 2;
    MPI_Allreduce( loc, e->mincoords.coords, 4, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( loc, e->maxcoords.coords, 4, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    for (i=0; i<4; i++)
      e->coords.coords[i] = loc[i];

    rca_get_max_dimension( &dims );
    e->maxtopocoords.coords[0] = dims.mesh_x + 1;
    e->maxtopocoords.coords[1] = dims.mesh_y + 1;
    e->maxtopocoords.coords[2] = dims.mesh_z + 1;
    e->maxtopocoords.coords[3] = 2;  /* 2 nodes per NIC */

    return 0;
}

