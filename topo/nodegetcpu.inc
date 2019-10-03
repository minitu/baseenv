/* -*- Mode: C; c-basic-offset:4 ; -*- */
#ifdef HAVE_UTMPX_H
#include <utmpx.h>
#endif
static int topoiGetNodeInfoGETCPU(int isMultithreaded, topoinfo_t *topo);

static int topoiGetNodeInfoGETCPU(int isMultithreaded, topoinfo_t *topo)
{
    int cpunum, chip=0;
    topoentry_t *e;

    cpunum = sched_getcpu();
    if (verbose) printf("cpu # = %d\n", cpunum);
    /* On some systems, the cores are numbered consecutively and the same
       number of cores per chip.  In this case, we can determine the chip
       by looking at the cpunum */
    if (cpunum >= 0 && coresPerChip > 0) {
	chip   = cpunum / coresPerChip;
	/* Normalize the cpunum as per chip */
	cpunum = cpunum - chip*coresPerChip;
    }
    if (cpunum >= 0) {
	e = topoiAllocEntry(topo);
	if (!e) return -1;
	e->dim = 1;
	if (coresPerChip > 0) {
	    e->maxcoords.coords[0] = coresPerChip;
	}
	else {
	    e->maxcoords.coords[0] = -1;  /* Unknown max */
	}
	e->topoType = TOPO_CORE;
	e->coords.coords[0] = cpunum;
	e->coords.coords[1] = -1;
    }
    if (coresPerChip > 0 && cpunum >= 0) {
	e = topoiAllocEntry(topo);
	if (!e) return -1;
	e->dim = 1;
	if (chipsPerNode > 0)
	    e->maxcoords.coords[0] = chipsPerNode;
	else
	    e->maxcoords.coords[0] = -1; /* Unknown max */
	e->topoType = TOPO_SOCKET;
	e->coords.coords[0] = chip;
	e->coords.coords[1] = -1;
    }
    return 0;
}

