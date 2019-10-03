/*
 * Test different options for stencil (halo) communication
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>   /* For isdigit */

#include "mpi.h"
#include "stencil.h"

static void getRangeList(const char *range, int *nvalsP, int **vals);

void myAbort(MPI_Comm comm, int rc, const char *msg)
{
    int wrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    fprintf(stderr, "[%d]: %s\n", wrank, msg);
    fflush(stderr);
    MPI_Abort(comm, rc);
}

#define NARGS 10
void getArgs(int argc, char *argv[],
	     probDesc *pd, procInfo *pi, stencilInfo *si, FILE **fp,
             int *maxcommtypes)
{
    int rank, size, i;
    int nvals, *vals;
    MPI_Comm comm;
    int args[NARGS], dims[2];

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* New argument processing:
       stencil --size xxx --energy nn --iters nn --px x --py y -o filename
               --norma --normadtype --normacontig
    */

    /* Set defaults */
    si->nSizes        = -1;
    si->nSize         = 0;
    pd->energy        = 1;           /* energy to be injected per iteration */
    si->niters        = 10;          /* number of iterations */
    dims[0]           = 0;
    dims[1]           = 0;
    MPI_Dims_create(size, 2, dims);
    pi->px            = dims[0];     /* 1st dim processes */
    pi->py            = dims[1];     /* 2nd dim processes */
    pi->useRMA        = 1;           /* Run RMA tests with contig data */
    pi->useRMAdtype   = 1;           /* Run RMA tests with datatypes */
    pi->skipCongruent = 1;           /* Skip comm congruent with WORLD */
    pi->withCart      = 0;
    /* maxcommtypes is set by default on entry to this routine */
    if (rank == 0) {
	for (i=1; i<argc; i++) {
	    char *v = argv[i];
	    if (*v == '-') v++;
	    if (*v == '-') v++;
	    if (strcmp(v, "size") == 0) {
		/* size is either a single value, a comma separated list, or
		   a range start:end:incr */
		getRangeList(argv[++i], &nvals, (int **)&vals);
		si->nSizes = nvals;
		si->nSize  = vals;
	    }
	    else if (strcmp(v, "energy") == 0) {
		pd->energy = atoi(argv[++i]);
	    }
	    else if (strcmp(v, "iters") == 0) {
		si->niters = atoi(argv[++i]);
	    }
	    else if (strcmp(v, "px") == 0) {
		pi->px = atoi(argv[++i]);
	    }
	    else if (strcmp(v, "py") == 0) {
		pi->py = atoi(argv[++i]);
	    }
	    else if (strcmp(v, "o") == 0) {
		*fp = fopen(argv[++i], "w");
		if (!*fp) {
		    myAbort(MPI_COMM_WORLD, 1, "Unable to open output file");
		}
	    }
	    else if (strcmp(v, "norma") == 0) {
		pi->useRMA      = 0;
		pi->useRMAdtype = 0;
	    }
	    else if (strcmp(v, "normadtype") == 0) {
		pi->useRMAdtype = 0;
	    }
	    else if (strcmp(v, "normacontig") == 0) {
		pi->useRMA      = 0;
	    }
	    else if (strcmp(v, "noskip") == 0) {
		pi->skipCongruent = 0;
	    }
	    else if (strcmp(v, "withcart") == 0) {
		pi->withCart      = 1;
	    }
	    else if (strcmp(v, "worldonly") == 0) {
		*maxcommtypes = 1;
	    }
	    else {
		fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
		fprintf(stderr, "%s --size xxx --energy nn --iters nn --px x --py y -o filename --norma --normadtype --normacontig --noskip --worldonly\n", argv[0]);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}

	/* Sanity check */
	if (si->nSizes <= 0) {
	    /* Create a default */
	    getRangeList("128,256,512,1024", &nvals, (int **)&vals);
	    si->nSizes = nvals;
	    si->nSize  = vals;
	}

	if(pi->px * pi->py != size)
	    myAbort(comm, 1, "px*py != size of COMM_WORLD");
	/* All sizes must evenly divide process count */
	for (i=0; i<nvals; i++) {
	    if (vals[i] % pi->py != 0)
		myAbort(comm, 2, "n mod py != 0");
	    if(vals[i] % pi->px != 0)
		myAbort(comm, 3, "n mod px != 0");
	}

	/* distribute arguments */
	args[0] = nvals; args[1] = pd->energy; args[2] = si->niters;
	args[3] = pi->px; args[4] = pi->py; args[5] = pi->useRMA;
	args[6] = pi->useRMAdtype; args[7] = pi->skipCongruent;
	args[8] = pi->withCart; args[9] = *maxcommtypes;
	MPI_Bcast(args, NARGS, MPI_INT, 0, comm);
	MPI_Bcast(vals, nvals, MPI_INT, 0, comm);
    }
    else {
	MPI_Bcast(args, NARGS, MPI_INT, 0, comm);
	nvals=args[0];  pd->energy=args[1]; si->niters=args[2];
	pi->px=args[3]; pi->py=args[4]; pi->useRMA=args[5];
	pi->useRMAdtype=args[6];
	pi->skipCongruent = args[7];
	pi->withCart      = args[8];
	*maxcommtypes     = args[9];
	si->nSize = (int *)calloc(nvals, sizeof(int));
	MPI_Bcast(si->nSize, nvals, MPI_INT, 0, comm);
    }
    si->nSizes = nvals;
}

/* first size is either a single value, a comma separated list, or
   a range start:end:incr.  Does not support a combination (e.g.,
   you can't do 10,15,20:40:10 */
static void getRangeList(const char *range, int *nvalsP, int **valsP)
{
    int nvals = 0, *vals;
    const char *p = range;
    int   isList = 1;

    /* Check for list or range */
    while (p && *p) {
	if (*p == ',') {
	    isList = 1;
	    break;
	}
	else if (*p == ':') {
	    isList = 2;
	    break;
	}
	p++;
    }

    p = range;
    if (isList == 1) {
	/* Find the number of values; check input for correctness */
	while (p && *p) {
	    if (isdigit(*p)) {
		while (isdigit(*p)) p++;
		if (!*p) {
		    /* Single value. */
		    nvals++;
		    break;
		}
		if (*p == ',') {
		    nvals++;
		    p++;
		    continue;
		}
		if (*p == ':') {
		    /* start:end:incr */
		    myAbort(MPI_COMM_WORLD, 1,
			    "range start:end:incr not implemented");
		    return;  /* Keep compiler happy */
		}
	    }
	    else {
		fprintf(stderr, "Invalid range format %s", range);
		myAbort(MPI_COMM_WORLD, 1, "n[,n2...]");
		return; /* Keep compiler happy */
	    }
	}
	p = range;
	*valsP = vals = (int *)calloc(nvals,sizeof(int));
	while (p && *p) {
	    int v = 0;
	    while (isdigit(*p)) {
		v = 10*v + (*p++ - '0');  /* Assumes ASCII */
	    }
	    *vals++ = v;
	    if (!*p) break;
	    if (*p++ != ',') {
		/* PANIC.  We should have caught this already */
		myAbort(MPI_COMM_WORLD, 1, "Internal error reading range");
		return; /* Keep compiler happy */
	    }
	}
    }
    else if (isList == 2) {
	/* range of the form start:end:incr */
	int start, end, incr, i;
	start = 0;
	while (*p && isdigit(*p))
	    start = start * 10 + (*p++ - '0');
	if (*p != ':') {
	    myAbort(MPI_COMM_WORLD, 1, "Unexpected character in range");
	    return; /* Keep compiler happy */
	}
	p++;
	end = 0;
	while (*p && isdigit(*p))
	    end = end * 10 + (*p++ - '0');

	incr = 1;
	if (*p == ':') {
	    incr = 0;
	    p++;
	    while (*p && isdigit(*p))
		incr = incr * 10 + (*p++ - '0');
	}
	if (*p != 0) {
	    myAbort(MPI_COMM_WORLD, 1, "Unexpected character in range");
	    return; /* Keep compiler happy */
	}

	nvals = 1 + (end - start + incr-1) / incr;
	*valsP = vals = (int *)calloc(nvals,sizeof(int));
	vals[0] = start;
	for (i=1; i<nvals; i++) {
	    vals[i] = vals[i-1] + incr;
	}
    }

    *nvalsP = nvals;
}

/* Q: Have multiple versions of this, to compare Cart_create with default
   and with SMP-aware reordering? */
void getDecomp(MPI_Comm comm, procInfo *pi, stencilInfo *si)
{
    int rank, size;
    int rx, ry;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    pi->comm = comm;

    /* determine my coordinates (x,y) -- r=x*a+y in the 2d processor array */
    rx = si->rx = rank % pi->px;
    ry = si->ry = rank / pi->px;
    /* determine my four neighbors */
    pi->north = (ry-1)*pi->px+rx; if(ry-1 < 0)   pi->north = MPI_PROC_NULL;
    pi->south = (ry+1)*pi->px+rx; if(ry+1 >= pi->py) pi->south = MPI_PROC_NULL;
    pi->west  = ry*pi->px+rx-1;   if(rx-1 < 0)   pi->west = MPI_PROC_NULL;
    pi->east  = ry*pi->px+rx+1;   if(rx+1 >= pi->px) pi->east = MPI_PROC_NULL;
    /* decompose the domain */
    si->bx = si->n/pi->px; /* block size in x */
    si->by = si->n/pi->py; /* block size in y */
    si->offx = rx*si->bx;  /* offset in x */
    si->offy = ry*si->by;  /* offset in y */
}


void getDecompCart(MPI_Comm comm, procInfo *pi, stencilInfo *si)
{
    int    rank, crank, size;
    int    pdims[2] = {0,0};
    int    periods[2] = {0,0};
    int    coords[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Override the process decomposition using DIMS_CREATE */
    MPI_Dims_create(size, 2, pdims);
    pi->px = pdims[0];
    pi->py = pdims[1];

    /* Select reordering to let system order ranks */
    MPI_Cart_create(comm, 2, pdims, periods, 1, &pi->comm);
    MPI_Comm_set_name(pi->comm, "Cart comm");

    /* determine my coordinates (x,y) -- r=x*a+y in the 2d processor array */
    MPI_Comm_rank(pi->comm, &crank);
    MPI_Cart_coords(pi->comm, crank, 2, coords);
    si->rx = coords[0];
    si->ry = coords[1];

    /* determine my four neighbors */
    MPI_Cart_shift(pi->comm, 0, 1, &pi->west, &pi->east);
    MPI_Cart_shift(pi->comm, 1, 1, &pi->north, &pi->south);

    /* decompose the domain */
    si->bx   = si->n/pi->px;   /* block size in x */
    si->by   = si->n/pi->py;   /* block size in y */
    si->offx = si->rx*si->bx;  /* offset in x */
    si->offy = si->ry*si->by;  /* offset in y */
}

void getDecompNodecart(MPI_Comm comm, procInfo *pi, stencilInfo *si, int wsock)
{
    int    rank, nrank, size;
    int    pdims[2] = {0,0};
    int    periods[2] = {0,0};
    int    coords[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Override the process decomposition using DIMS_CREATE */
    MPI_Dims_create(size, 2, pdims);
    pi->px = pdims[0];
    pi->py = pdims[1];

    /* Select reordering to let system order ranks */

    if (wsock) {
        MPIX_SetUseSocket(1);
    }
    MPIX_Nodecart_create(comm, 2, pdims, periods, 1, &pi->comm);
    if (wsock) {
        MPIX_SetUseSocket(0);
        MPI_Comm_set_name(pi->comm, "Nodecart-h");
    }
    else
        MPI_Comm_set_name(pi->comm, "Nodecart");

    /* determine my coordinates (x,y) -- r=x*a+y in the 2d processor array */
    MPI_Comm_rank(pi->comm, &nrank);
    MPIX_Nodecart_coords(pi->comm, nrank, 2, coords);
    si->rx = coords[0];
    si->ry = coords[1];

    /* determine my four neighbors */
    MPIX_Nodecart_shift(pi->comm, 0, 1, &pi->west, &pi->east);
    MPIX_Nodecart_shift(pi->comm, 1, 1, &pi->north, &pi->south);

    /* Sanity check */
    { int err = 0;
	if (pi->west != MPI_PROC_NULL && (pi->west < 0 || pi->west >= size)) {
	    fprintf(stderr, "%d: west = %d, not in [0,%d]\n", rank, pi->west,
		    size-1);
	    fflush(stderr);
	    err++;
	}
	if (pi->east != MPI_PROC_NULL && (pi->east < 0 || pi->east >= size)) {
	    fprintf(stderr, "%d: east = %d, not in [0,%d]\n", rank, pi->east,
		    size-1);
	    fflush(stderr);
	    err++;
	}
	if (pi->north != MPI_PROC_NULL &&
	    (pi->north < 0 || pi->north >= size)) {
	    fprintf(stderr, "%d: north = %d, not in [0,%d]\n", rank,
		    pi->north, size-1);
	    fflush(stderr);
	    err++;
	}
	if (pi->south != MPI_PROC_NULL &&
	    (pi->south < 0 || pi->south >= size)) {
	    fprintf(stderr, "%d: south = %d, not in [0,%d]\n", rank, pi->south,
		    size-1);
	    fflush(stderr);
	    err++;
	}
	if (err > 0) MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* decompose the domain */
    si->bx   = si->n/pi->px;   /* block size in x */
    si->by   = si->n/pi->py;   /* block size in y */
    si->offx = si->rx*si->bx;  /* offset in x */
    si->offy = si->ry*si->by;  /* offset in y */
}


/* Initialize the structure used to record times for stencil operations */
void initST(stencilTrialTime *st, int ntimes)
{
    int i;

    st->navail = ntimes;
    st->ncur   = -1;
    st->st = (stencilTime *)calloc(ntimes, sizeof(stencilTime));
    if (!st->st) myAbort(MPI_COMM_WORLD, 1,
		      "Unable to allocate memory for timing data");
    for (i=0; i<ntimes; i++) {
	st->st[i].commInit     = 0;
	st->st[i].commStart    = 0;
	st->st[i].commComplete = 0;
	st->st[i].commPack     = 0;
	st->st[i].commUnpack   = 0;
	st->st[i].compInterior = 0;
	st->st[i].compBndy     = 0;
	st->st[i].compUpdate   = 0;
    }
}

/* Initialize three heat sources */
void initProbDesc(probDesc *pd, stencilInfo *si)
{
    int i, n = si->n;     /* n is the global mesh size. */
    int sources[NSOURCES][2]; /* = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*8/9}}; */

    sources[0][0] = n/2;
    sources[0][1] = n/2;
    sources[1][0] = n/3;
    sources[1][1] = n/3;
    sources[2][0] = n*4/5;
    sources[2][1] = n*8/9;

    pd->locnsources = 0; /* number of sources in my area */
    for (i=0; i<NSOURCES; ++i) {
	/* determine which sources are in my patch */
	int locx = sources[i][0] - si->offx;
	int locy = sources[i][1] - si->offy;
	if(locx >= 0 && locx < si->bx && locy >= 0 && locy < si->by) {
	    pd->locsources[pd->locnsources][0] = locx+1; /* offset by halo zone */
	    pd->locsources[pd->locnsources][1] = locy+1; /* offset by halo zone */
	    pd->locnsources++;
	}
    }
}

void allocMesh(stencilInfo *si)
{
    /* Allocate memory in a single block; this simplifies some of the
       versions (esp. RMA). */
    size_t asize = (si->bx+2)*(si->by+2);
    si->aold = (double*)calloc(2*asize,sizeof(double)); /* 1-wide halo zones! */
    si->anew = si->aold + asize;           /* 1-wide halo zones! */
}

void freeMesh(stencilInfo *si)
{
    free(si->aold);
}

void initMesh(probDesc *pd, stencilInfo *si)
{
    initMeshBase(pd, si->bx, si->by, si->aold, si->anew);
}

void initMeshBase(probDesc *pd, int bx, int by, double *restrict aold,
		   double *restrict anew)
{
    int i, j;

    /* We must initialize the halo because that includes the boundary
       conditions, which are a==0 (fixed Temperature edges)
       for this example */
    for(j=0; j<=by+1; ++j) {
	for(i=0; i<=bx+1; ++i) {
	    anew[ind(i,j)] = 0;
	    aold[ind(i,j)] = 0;
	}
    }

    /* Initialize heat sources */
    for(i=0; i<pd->locnsources; ++i) {
	aold[ind(pd->locsources[i][0],pd->locsources[i][1])] +=
	    pd->energy; /* heat source */
    }
}

# if 0
/* These are not yet correct because they don't work correctly if there is
   more than one process */
#include <math.h>
void printarr(FILE *fp, double *a, int n, int bx) {
    /* does nothing right now, should record each "frame" as image */
    int size = 5;
    int i, j;

    size = 1024 / n;
    if (size < 1) size = 1;
    else if (size > 5) size = 5;
    /* Without the height and width here, browsers seem to use a default that
       is too small */
    fprintf(fp, "<html>\n<body>\n<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"%d\" height=\"%d\">",size*n, size*n);

    fprintf(fp, "\n<rect x=\"0\" y=\"0\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(0,0,0);stroke:rgb(0,0,0)\"/>", size*n, size*n);
    for(i=1; i<n+1; ++i)
	for(j=1; j<n+1; ++j) {
	    int rgb = (a[ind(i,j)] > 0) ? rgb = (int)round(255.0*a[ind(i,j)]) : 0.0;
	    if(rgb>255) rgb=255;
	    if(rgb) fprintf(fp, "\n<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(%i,0,0);stroke:rgb(%i,0,0)\"/>", size*(i-1), size*(j-1), size, size, rgb, rgb);
	}
    fprintf(fp, "</svg>\n</body>\n</html>");
}
#endif
void printSoln(stencilInfo *si, double *v, const char *name)
{
    char fname[256];
    FILE *fp;
    int  bx = si->bx, by = si->by, n = si->n;

    /* Form the output file name from the label, size
     */
    snprintf(fname, sizeof(fname), "heat-%s-%d.svg", name, n);

    fp = fopen(fname, "a");
    if (!fp) {
	myAbort(MPI_COMM_WORLD, 1, "Unable to open solution output file");
    }
    /* To make this correct for the parallel case, v should be a
       global version of the array, with each process contributing its
       piece */
#if 0
    printarr(fp, v, n, bx);
#else
    /* This just outputs each grid to the designated file */
    {
	int wrank, wsize, r, i, j;
	MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	for (r=0; r<wsize; r++) {
	    if (r == wrank) {
		fprintf(fp,"Rank %d proc[%d,%d], grid offset[%d,%d]:\n",
			r, si->rx, si->ry, si->offx, si->offy);
		/* Include boundaries and any ghost cells */
		for(j=0; j<by+2; ++j) {
		    for(i=0; i<bx+2; ++i) {
			fprintf(fp, "%.2e ", v[ind(i,j)]);
		    }
		    fprintf(fp, "\n");
		}
		fflush(fp);
	    }
	    /* A faster version could use the Seq routines, but this should be
	       adequate for debugging and requires no external routines */
	    MPI_Barrier(MPI_COMM_WORLD);
	}
    }
#endif
    fclose(fp);
}

#ifdef HAVE_TOPOINIT
#include "topoinfo.h"
static topoinfo_t *tinfo=0;

void initTopology(int v)
{
    topoInit(v, &tinfo);
}

/* This is a temporary that prints data - should be replaced by a
   routine to return data, and have a separate print routine */
void getNodeCommunicationStats(procInfo *pi, FILE *fp)
{
    topodist_t *dt = 0;
    int        sends[5], ns=0;
    MPI_Comm   nodecomm;
    int        nOnNode, nOffNode, totalOffNode;
    int        nrank, nsize, nnodes, wrank;
    int nodeOnNode, nodeOffNode, nn[2], navg[2], tn[2], tavg;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (pi->south != MPI_PROC_NULL) sends[ns++] = pi->south;
    if (pi->north != MPI_PROC_NULL) sends[ns++] = pi->north;
    if (pi->east  != MPI_PROC_NULL) sends[ns++] = pi->east;
    if (pi->west  != MPI_PROC_NULL) sends[ns++] = pi->west;

    topodistInit(pi->comm, ns, sends, ns, sends, tinfo, &dt);
    topodistNodeCommInfo(dt, MPI_COMM_WORLD, &nOnNode, &nOffNode,
			 &nodecomm, &totalOffNode);
    MPI_Comm_rank(nodecomm, &nrank);
    MPI_Comm_size(nodecomm, &nsize);
    nnodes = (nrank == 0);
    MPI_Allreduce(MPI_IN_PLACE, &nnodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    /*
     * Also compute: min, max, average:
     * Number of sends per process out of node, in node
     * Number of sends per node out of node, in node
     */
    /* Per node numbers */
    MPI_Allreduce(&nOnNode, &nodeOnNode, 1, MPI_INT, MPI_SUM, nodecomm);
    MPI_Allreduce(&nOffNode, &nodeOffNode, 1, MPI_INT, MPI_SUM, nodecomm);

    /* min and max over nodes */
    nn[0] = nodeOffNode;
    nn[1] = -nodeOffNode;
    MPI_Allreduce(MPI_IN_PLACE, nn, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (nrank == 0) { navg[0] = nodeOffNode; navg[1] = nodeOnNode;
	tn[0] = totalOffNode; tn[1] = - totalOffNode; tavg = totalOffNode; }
    else { navg[0] = 0; navg[1] = 0;
	tn[0] = 0; tn[1] = -1000000000; tavg = 0.0; }
    MPI_Allreduce(MPI_IN_PLACE, navg, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    navg[0] /= nnodes;
    navg[1] /= nnodes;
    MPI_Allreduce(MPI_IN_PLACE, tn, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tavg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    tavg /= nnodes;

    MPI_Comm_free(&nodecomm);
    topodistFree(dt);

    /* Eventually make this use verbose */
    if (0 && nrank == 0) {
	printf("[%d]: (node size = %d) totalOffNode = %d\n", wrank, nsize,
	       totalOffNode);
    }
    /* Over all nodes */
    if (wrank == 0) {
	fprintf(fp, "Processes communicating off node (max,min,avg): %d,%d,%d\n",
		tn[0], -tn[1], tavg);
	fprintf(fp, "Sends/node Offnode(max,min,avg): %d,%d,%d\n",
		nn[0],-nn[1], navg[0]);
	fflush(fp);
    }
}

void endTopology(void)
{
    topoFinalize(&tinfo);
}
#endif
