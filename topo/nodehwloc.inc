/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include "hwloc.h"
/* Could use one process per node and share information.  For that,
   need to run the nodetopo information first */
static int topoiGetNodeInfoHWLOC(int isMultithreaded, topoinfo_t *topo);

static int topoiGetNodeInfoHWLOC(int isMultithreaded, topoinfo_t *topo)
{
  hwloc_topology_t topology;
  hwloc_cpuset_t   cpuset;
  hwloc_thread_t   threadid;
  hwloc_pid_t      pid;
  hwloc_obj_type_t obtype;
  int              topodepth, err=0;
  unsigned         depth;

  err = hwloc_topology_init(&topology);
  err = hwloc_topology_load(topology);
  topodepth = hwloc_topology_get_depth(topology);
  /* Extract information about the node topology */
  if (verbose) printf( "depth = %d\n", topodepth );
  for (depth=0; depth<topodepth; depth++) {
    int i, nobj;
    char string[256];
    topoentry_t *e=0;
    if (verbose) printf("Objects at level %u\n", depth);
    obtype = hwloc_get_depth_type(topology,depth);
    nobj   = hwloc_get_nbobjs_by_depth(topology,depth);
    /*
     * For the purposes of this function, we want to pay attention to
     *   HWLOC_OBJ_NODE
     *   HWLOC_OBJ_SOCKET
     *   HWLOC_OBJ_CORE
     *   (maybe) HWLOC_OBJ_PU
     *
     *   HWLOC_OBJ_CACHE may also be interesting
     *
     * The rule here should be to return the number of distinct
     * processing levels (not affinitiy levels, which hwloc is doing).
     * Thus, we can skip the cache and other objects.
     * ToDo: Nodes with mixed processor types may need special handling,
     * e.g., nodes with CPU and GPU chips.
     */
    switch (obtype) {
    case HWLOC_OBJ_NODE:
      if (verbose) printf( "%d node\n", nobj );
      e = topoiAllocEntry(topo);
      if (!e) return -1;
      e->dim = 1;
      e->maxcoords.coords[0] = nobj;
      e->topoType = TOPO_NODE;
      break;
    case HWLOC_OBJ_SOCKET:
      if (verbose) printf( "%d sockets\n", nobj );
      e = topoiAllocEntry(topo);
      if (!e) return -1;
      e->dim = 1;
      e->maxcoords.coords[0] = nobj;
      e->topoType = TOPO_SOCKET;
      break;
    case HWLOC_OBJ_CORE:
      if (verbose) printf( "%d cores\n", nobj );
      e = topoiAllocEntry(topo);
      if (!e) return -1;
      e->dim = 1;
      e->maxcoords.coords[0] = nobj;
      e->topoType = TOPO_CORE;
      cpuset = hwloc_bitmap_alloc();
      err = hwloc_get_last_cpu_location(topology,cpuset,0);
      if (err) {
	if (verbose)
	  perror( "get_location failed:" );
      }
      else {
	unsigned int index, nfound=0, foundindex;
	hwloc_bitmap_foreach_begin(index,cpuset);
	if (verbose) printf("cpuset contains index %u\n", index );
	nfound++;
	foundindex = index;
	hwloc_bitmap_foreach_end();
	if (nfound == 1) {
	  e->coords.coords[0] = foundindex;
	  e->coords.coords[1] = -1;
	}
      }
  if (isMultithreaded) {
    err = hwloc_get_thread_cpubind(topology,threadid,cpuset,0);
  }
  hwloc_bitmap_free(cpuset);
      break;
    default:
      break;
    }
    if (verbose) {
      for (i=0; i<nobj; i++) {
	hwloc_obj_snprintf(string, sizeof(string), topology,
			   hwloc_get_obj_by_depth(topology, depth, i),
			   "#", 0);
	printf("\t%d\t%s\n", i, string);
      }
    }
  }

  /* Cleanup */
  hwloc_topology_destroy(topology);
  return err;
}
