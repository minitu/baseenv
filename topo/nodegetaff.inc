/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include <sched.h>

/*
 * Use sched_getaffinity to find the thread number.
 * No chip info, so apply to a node.
 */
static int topoiGetNodeInfoGETAFFINITY(int isMultithreaded, topoinfo_t *topo);

static int topoiGetNodeInfoGETAFFINITY(int isMultithreaded, topoinfo_t *topo)
{
    topoentry_t *e;
    int       i, j, corenum, chip;
    cpu_set_t coremask;

    (void)sched_getaffinity(0, sizeof(coremask), &coremask);

    corenum = -1;
    for (i = 0; i < CPU_SETSIZE; i++) {
	if (CPU_ISSET(i, &coremask)) {
	    /* Found a core on which this thread can run.  Now
	       we check to see if that is unique.  If not, then
	       we can't determine on which core we expect the
	       thread to run. */
	    int run = 0;
	    for (j = i + 1; j < CPU_SETSIZE; j++) {
		if (CPU_ISSET(j, &coremask)) run++;
	    }
	    if (!run)
		corenum = i;
	    break;
	}
    }
    if (corenum >= 0 && coresPerChip > 0) {
	chip   = corenum / coresPerChip;
	/* Normalize the cpunum as per chip */
	corenum = corenum - chip*coresPerChip;
    }
    if (corenum >= 0) {
	e = topoiAllocEntry(topo);
	if (!e) return 1;
	e->dim = 1;
	if (coresPerChip > 0) {
	    e->maxcoords.coords[0] = coresPerChip;
	}
	else {
	    e->maxcoords.coords[0] = -1;  /* Unknown max */
	}
	e->topoType            = TOPO_CORE;
	e->coords.coords[0]    = corenum;
	e->coords.coords[1]    = -1;
    }
    if (coresPerChip > 0 && corenum >= 0) {
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


