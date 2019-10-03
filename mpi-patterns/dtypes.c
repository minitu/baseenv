#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Needed for restrict definition if using an old C compiler */
#include "baseenv.h"

/* See
/Users/gropp/projects/Papers/2011/performance_requirements/progs/mpi_dttest
 for more tests
*/
/*
 * Explore the performance of MPI derived datatypes, in different uses:
 *
 *   Pack/Unpack
 *   Point-to-point communication
 *   RMA
 *   Collective communication
 *
 * In each case, the performance is compared to simple "user" code that
 * uses a temporary buffer.  Experiments are run with different patterns and
 * sizes to explore the impact of different communication regiemes (e.g.,
 * eager and rendezvous) and memory hierarchy.
 *
 * Also explores greater range of datatype parameters, including block
 * sizes.
 */

/* Record times for a single communication, using user-pack and datatypes.
   In the special case of no communication (pack or unpack), the "comm"
   field is used (not pack or unpack) */
typedef struct {
    double *pack,     /* Time for user pack code */
	   *comm,     /* Time for user communication */
           *unpack,   /* Time for user unpack code */
           *dtcomm;   /* Time for communication with datatype */
} commTimes2;
/*
   This structure is for each datatype, and contains arrays of times for
   each component of the test
*/
typedef struct {
    enum { DT_NONE, DT_CONTIG, DT_VECTOR, DT_INDEX, DT_BLOCKINDEX, DT_STRUCT,
	   DT_RESIZED } dttype;
    int  dtcount;        /* Number of instances of the datatype to use */
    int  parms[4];       /* Parameters to describe datatype */
    int  uparms[4];      /* Parmaeters used for user-version of pack/unpack */
    int  *aparms;
    MPI_Datatype basetype; /* Base type to use */
    const char *label;   /* Short string for case; see below */
    MPI_Aint nbytes;     /* Amount of data moved in each operation */
    int    navail, ncur; /* sizeof timing arrays and next free index */
    double dtcreate;     /* Time to create and commit datatype */
    commTimes2 pack,      /* Pack */
	unpack,          /* Unpack */
	pt2pt,           /* Point to point communication */
	rmaPut,          /* RMA Put with datatype */
	rmaGet;          /* RMA Get with datatype */
} dtypeTrials;

/* Label values:
   [C/V/I/B/S] [-/C/V/S] [l] [m] [-s]
   where

   C/V/I/B/S/R - Datatype is Contiguous, Vector, Indexed, Block-indexed,
                 struct, or resized.  Second set is the same,
                 with - meaning base type was predefined

   l         - size of base type, e.g. 8 for double
   m         - count (? use s/m/l for small, medium,large?)
   s         - stride (? 1/s/m/l, allow 1 for stride, as well as - (N/A))

 */

void initTestTypes(dtypeTrials dtimes[], int dtimesSize, MPI_Datatype basetype);
void initDT(dtypeTrials *dt, int ntimes);
void myAbort(MPI_Comm comm, int rc, const char *msg);
void syncPartners(MPI_Comm comm, int rank, int master, int partner);
void initWin(MPI_Aint winsize, int basesize, MPI_Win *win_p, void *winmem);
void printResults(dtypeTrials dtimes[]);
void getMinAvgTimes(double *t, int n, double *tmin_p, double *tavg_p);
void createDatatype(dtypeTrials *dtimes, MPI_Datatype basetype,
		    MPI_Datatype *dt);
int allocBuffers(MPI_Datatype dt, int dtcount, MPI_Comm comm,
		 void *inbuf_p, void *pkoutbuf_p, void *outbuf_p,
		 void *upkbuf_p,
		 MPI_Aint *insize_p, MPI_Aint *pksize_p, MPI_Aint *outsize_p);
void freeBuffers(void *inbuf, void *outbuf, void *pkbuf, void *upkbuf);
void UserVecPack(const double *restrict inbuf, double *restrict pkbuf,
		 int count, int stride, int blklen, double *tpack, int *ucount);
void UserVecUnpack(const double *restrict pkbuf, double *restrict outbuf,
		   int count, int stride, int blklen,
		   double *tunpack);

static int verbose = 0;
static int withRma = 0;

#define MAX_TEST_DTYPES 100

int main(int argc, char *argv[])
{
    MPI_Datatype dt, basetype;
    MPI_Comm     comm;
    MPI_Win      win;
    double      *winmem;
    MPI_Aint     insize, outsize, pksize;
    double      *inbuf, *outbuf, *pkoutbuf, *upkbuf;
    int          position, typesize, ucount;
    double       t;
    dtypeTrials  dtimes[MAX_TEST_DTYPES];
    int          tidx=0, dtcount;
    int          nt, ntimes = 10, r, rcount=100;
    int          wrank, wsize, rank, master, partner;

    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    master  = 0;
    partner = wsize - 1;

    basetype = MPI_DOUBLE;
    MPI_Type_size(basetype, &typesize);
    initTestTypes(dtimes, MAX_TEST_DTYPES, basetype);

    for (tidx=0; dtimes[tidx].dttype != DT_NONE; tidx++) {
	int dtsize;

	if (verbose && wrank == 0)
	    printf("Type test %d: %s\n", tidx, dtimes[tidx].label);
	initDT(&dtimes[tidx], ntimes);

	createDatatype(&dtimes[tidx], dtimes[tidx].basetype, &dt);

	dtcount = dtimes[tidx].dtcount;
	allocBuffers(dt, dtcount, comm, &inbuf, &pkoutbuf, &outbuf, &upkbuf,
		     &insize, &pksize, &outsize);
	MPI_Type_size(dt, &dtsize);
	dtimes[tidx].nbytes = dtcount * (MPI_Aint)dtsize;

	if (withRma) {
	    MPI_Aint extent, lb;
	    MPI_Type_get_extent(dt, &lb, &extent);
	    extent   *= dtcount;
	    initWin(extent, typesize, &win, &winmem);
	}

	for (nt=0; nt<ntimes; nt++) {
	    int packedSize;
	    dtimes[tidx].ncur = nt;

	    if (verbose && wrank == 0)
		printf("\tTest %d\n", nt);

	    /* Perform operations and time them */
	    if (verbose && wrank == 0) printf("\tPack: %ld %d\n",
					      (long)pksize, dtcount);
	    t = MPI_Wtime();
	    for (r=0; r<rcount; r++) {
		position = 0;
		MPI_Pack(inbuf, dtcount, dt, pkoutbuf, pksize, &position, comm);
	    }
	    dtimes[tidx].pack.dtcomm[nt]  = (MPI_Wtime() - t)/rcount;

	    if (verbose && wrank == 0) printf("\tUnpack: %ld %ld\n",
					      (long)pksize, (long)position);
	    packedSize = position;
	    t = MPI_Wtime();
	    for (r=0; r<rcount; r++) {
		position   = 0;
		MPI_Unpack(pkoutbuf, packedSize, &position,
			   outbuf, dtcount, dt, comm);
	    }
	    dtimes[tidx].unpack.dtcomm[nt] = (MPI_Wtime() - t)/rcount;

	    if (verbose && wrank == 0) printf("\tPt2pt:\n");
	    syncPartners(comm, rank, master, partner);
	    t = MPI_Wtime();
	    for (r=0; r<rcount; r++) {
		if (rank == master) {
		    MPI_Send(inbuf, dtcount, dt, partner, 2, comm);
		    MPI_Recv(outbuf, dtcount, dt, partner, 3, comm,
			     MPI_STATUS_IGNORE);
		}
		else if (rank == partner) {
		    MPI_Recv(outbuf, dtcount, dt, master, 2, comm,
			     MPI_STATUS_IGNORE);
		    MPI_Send(inbuf, dtcount, dt, master, 3, comm);
		}
	    }
	    dtimes[tidx].pt2pt.dtcomm[nt] = (MPI_Wtime() - t)/rcount;

	    if (withRma) {
		if (verbose && wrank == 0) printf("\tRMA:\n");
		syncPartners(comm, rank, master, partner);
		t = MPI_Wtime();
		MPI_Win_fence(0,win);
		if (rank == master)
		    MPI_Put(inbuf, dtcount, dt, partner, 0, dtcount, dt, win);
		MPI_Win_fence(0,win);
		dtimes[tidx].rmaPut.dtcomm[nt] = MPI_Wtime() - t;

		syncPartners(comm, rank, master, partner);
		t = MPI_Wtime();
		MPI_Win_fence(0,win);
		if (rank == master)
		    MPI_Get(outbuf, dtcount, dt, partner, 0, dtcount, dt, win);
		MPI_Win_fence(0,win);
		dtimes[tidx].rmaGet.dtcomm[nt] = MPI_Wtime() - t;
	    }

	    /* Perform the operations without datatypes */
	    switch (dtimes[tidx].dttype) {
  	    case DT_CONTIG:
		break;
	    case DT_VECTOR:
		UserVecPack(inbuf, upkbuf,
			    dtimes[tidx].uparms[0], /* count */
			    dtimes[tidx].uparms[2], /* stride */
			    dtimes[tidx].uparms[1], /* blklen */
			    &dtimes[tidx].pack.comm[nt], &ucount);
		dtimes[tidx].pt2pt.pack[nt] = dtimes[tidx].pack.comm[nt];
		break;
	    case DT_RESIZED:
		UserVecPack(inbuf, upkbuf,
			    dtimes[tidx].dtcount,   /* count */
			    dtimes[tidx].uparms[2], /* stride */
			    dtimes[tidx].uparms[1], /* blklen */
			    &dtimes[tidx].pack.comm[nt], &ucount);
		dtimes[tidx].pt2pt.pack[nt] = dtimes[tidx].pack.comm[nt];
		break;
	    default:
		break;
	    }
	    /* Sanity check on packing: ucount*sizeof(element type) must
	       equal the nbytes for the datatype */
	    if (ucount * sizeof(double) != dtimes[tidx].nbytes) {
		fprintf(stderr, "Ucount = %d, size = %d, but nbytes = %ld\n",
			ucount, ucount * 8, (long)dtimes[tidx].nbytes);
		myAbort(MPI_COMM_WORLD, 1, "Incorrect user packed size");
	    }

	    syncPartners(comm, rank, master, partner);
	    t = MPI_Wtime();
	    if (rank == master) {
		MPI_Send(upkbuf, ucount, basetype, partner, 2, comm);
		MPI_Recv(upkbuf, ucount, basetype, partner, 3, comm,
			 MPI_STATUS_IGNORE);
	    }
	    else if (rank == partner) {
		MPI_Recv(upkbuf, ucount, basetype, master, 2, comm,
			 MPI_STATUS_IGNORE);
		MPI_Send(upkbuf, ucount, basetype, master, 3, comm);
	    }
	    dtimes[tidx].pt2pt.comm[nt] = MPI_Wtime() - t;

	    switch (dtimes[tidx].dttype) {
  	    case DT_CONTIG:
		break;
	    case DT_VECTOR:
		UserVecUnpack(upkbuf, outbuf,
			      dtimes[tidx].parms[0], /* count */
			      dtimes[tidx].parms[2], /* stride */
			      dtimes[tidx].parms[1], /* blklen */
			      &dtimes[tidx].unpack.comm[nt]);
		dtimes[tidx].pt2pt.unpack[nt] = dtimes[tidx].unpack.comm[nt];
		break;
	    case DT_RESIZED:
		UserVecUnpack(upkbuf, outbuf,
			      dtimes[tidx].dtcount,                 /* count */
			      dtimes[tidx].parms[0]/sizeof(double), /* stride */
			      typesize/sizeof(double),              /* blklen */
			      &dtimes[tidx].unpack.comm[nt]);
		dtimes[tidx].pt2pt.unpack[nt] = dtimes[tidx].unpack.comm[nt];
		break;
	    default:
		break;
	    }
	} /* for nt */
	/* Free buffers and windows */
	freeBuffers(inbuf, outbuf, pkoutbuf, upkbuf);
	if (withRma) {}
	MPI_Type_free(&dt);
    } /* for tidx */

    /* Make the report */
    printResults(dtimes);

    MPI_Finalize();
    return 0;
}

/* Initialize the structure used to record times for stencil operations */
void initCommTimes(commTimes2 *op, int ntimes);

void initDT(dtypeTrials *dt, int ntimes)
{
    dt->navail = ntimes;
    dt->ncur   = -1;
    initCommTimes(&dt->pack, ntimes);
    initCommTimes(&dt->unpack, ntimes);
    initCommTimes(&dt->pt2pt, ntimes);
    initCommTimes(&dt->rmaPut, ntimes);
    initCommTimes(&dt->rmaGet, ntimes);
}

void initCommTimes(commTimes2 *op, int ntimes)
{
    op->pack   = (double *)malloc(ntimes*sizeof(double));
    op->comm   = (double *)malloc(ntimes*sizeof(double));
    op->unpack = (double *)malloc(ntimes*sizeof(double));
    op->dtcomm = (double *)malloc(ntimes*sizeof(double));
    if (!op->pack || !op->comm || !op->unpack || !op->dtcomm) {
	myAbort(MPI_COMM_WORLD, 1, "Unable to allocate commTimes arrays");
    }
}

void getMinAvgTimes(double *t, int n, double *tmin_p, double *tavg_p)
{
    int i;
    double tmin=1e10, tavg=0;
    for (i=0; i<n; i++) {
	tavg += t[i];
	if (t[i] < tmin) tmin = t[i];
    }
    *tmin_p = tmin;
    *tavg_p = tavg / n;
}

void myAbort(MPI_Comm comm, int rc, const char *msg)
{
    int wrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    fprintf(stderr, "[%d]: %s\n", wrank, msg);
    fflush(stderr);
    MPI_Abort(comm, rc);
}

void syncPartners(MPI_Comm comm, int rank, int master, int partner)
{
    if (verbose && (rank == master || rank == partner))
	printf("\tStarting syncPartners (%d,%d)\n", master, partner);
    if (rank == master) {
        MPI_Send(NULL, 0, MPI_INT, partner, 0, comm);
	MPI_Recv(NULL, 0, MPI_INT, partner, 1, comm, MPI_STATUS_IGNORE);
    }
    else if (rank == partner) {
        MPI_Recv(NULL, 0, MPI_INT, master, 0, comm, MPI_STATUS_IGNORE);
	MPI_Send(NULL, 0, MPI_INT, master, 1, comm);
    }
}

/* winmen is really a void** pointer to a pointer to hold the window memory */
void initWin(MPI_Aint winsize, int basesize, MPI_Win *win_p, void *winmem)
{
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_lock", "true");
    MPI_Info_set(info, "same_disp_unit", "true");
    MPI_Info_set(info, "same_size", "true");
    MPI_Info_set(info, "accumulate_ops", "same_op");
    MPI_Win_allocate(winsize, basesize, info, MPI_COMM_WORLD, winmem, win_p);
    MPI_Info_free(&info);
}

/*
 * Common Datatype Test Routines
a */

int allocBuffers(MPI_Datatype dt, int dtcount, MPI_Comm comm,
		 void *inbuf_p, void *pkoutbuf_p, void *outbuf_p,
		 void *upkbuf_p,
		 MPI_Aint *insize_p, MPI_Aint *pksize_p, MPI_Aint *outsize_p)
{
    MPI_Aint lb, extent;
    int      i, pksize, insize, typesize;
    double   *inbuf, *outbuf, *pkoutbuf, *upkbuf;

    /* Allocate buffers and windows for type operations */
    /* The idea here is that space is allocated to fit the datatype.
     */
    MPI_Type_get_extent(dt, &lb, &extent);
    extent   *= dtcount;
    inbuf    = (double *)malloc(extent);
    MPI_Pack_size(dtcount, dt, comm, &pksize);
    pkoutbuf = (double *)malloc(pksize);
    outbuf   = (double *)malloc(extent);
    MPI_Type_size(dt, &typesize);
    typesize *= dtcount;
    upkbuf   = (double *)malloc(typesize);

    if (verbose) {
	printf("\tAlloc: %ld %ld %ld\n", (long)extent, (long)pksize,
	       (long)typesize);
    }

    /* Initilize buffers */
    insize = extent / sizeof(double);
    for (i=0; i<insize; i++) {
	inbuf[i] = i;
	outbuf[i] = -i;
    }
    insize = pksize / sizeof(double);
    for (i=0; i<insize; i++)
	pkoutbuf[i] = -i;
    insize = typesize / sizeof(double);
    for (i=0; i<insize; i++)
	upkbuf[i] = -i;

    *(void **)inbuf_p     = (void *)inbuf;
    *(void **)outbuf_p    = (void *)outbuf;
    *(void **)pkoutbuf_p  = (void *)pkoutbuf;
    *(void **)upkbuf_p    = (void *)upkbuf;
    *insize_p  = extent;
    *outsize_p = extent;
    *pksize_p  = pksize;

    return 0;
}

void freeBuffers(void *inbuf, void *outbuf, void *pkbuf, void *upkbuf)
{
    if (verbose)
	printf("\tFree buffers\n");
    free(inbuf); free(outbuf); free(pkbuf); free(upkbuf);
}

/*
 * Common Vector Test Routines
 */

void UserVecPack(const double *restrict inbuf, double *restrict pkbuf,
		 int count, int stride, int blklen, double *tpack, int *ucount)
{
    double t;
    int    i, j;

    if (verbose) {
	printf("\tUvecpack %d %d %d\n", count, stride, blklen);
    }

    if (blklen == 1) {
	t = MPI_Wtime();
	for (i=0; i<count; i++) {
	    pkbuf[i] = inbuf[i*stride];
	}
	*tpack = MPI_Wtime() - t;
	*ucount = count;
    }
    else if (blklen == 4) {
	t = MPI_Wtime();
	j = 0;
	for (i=0; i<count; i++) {
	    pkbuf[j]   = inbuf[i*stride];
	    pkbuf[j+1] = inbuf[i*stride+1];
	    pkbuf[j+2] = inbuf[i*stride+2];
	    pkbuf[j+3] = inbuf[i*stride+3];
	    j+=4;
	}
	*tpack = MPI_Wtime() - t;
	*ucount = count*4;
    }
    else {
	t = MPI_Wtime();
	for (i=0; i<count; i++) {
	    for (j=0; j<blklen; j++)
		pkbuf[i*blklen+j] = inbuf[i*stride+j];
	}
	*tpack = MPI_Wtime() - t;
	*ucount = count * blklen;
    }
}

void UserVecUnpack(const double *restrict pkbuf, double *restrict outbuf,
		   int count, int stride, int blklen,
		   double *tunpack)
{
    double t;
    int    i, j;

    if (verbose) {
	printf("\tUvecunpack %d %d %d\n", count, stride, blklen);
    }

    if (blklen == 1) {
	t = MPI_Wtime();
	for (i=0; i<count; i++) {
	    outbuf[i*stride] = pkbuf[i];
	}
	*tunpack = MPI_Wtime() - t;
    }
    else if (blklen == 4) {
	t = MPI_Wtime();
	j=0;
	for (i=0; i<count; i++) {
	    outbuf[i*stride]   = pkbuf[j];
	    outbuf[i*stride+1] = pkbuf[j+1];
	    outbuf[i*stride+2] = pkbuf[j+2];
	    outbuf[i*stride+3] = pkbuf[j+3];
	    j+=4;
	}
	*tunpack = MPI_Wtime() - t;
    }
    else {
	t = MPI_Wtime();
	for (i=0; i<count; i++) {
	    for (j=0; j<blklen; j++)
		outbuf[i*stride+j] = pkbuf[i*blklen+j];
	}
	*tunpack = MPI_Wtime() - t;
    }
}

/* Output is
   type p1 p2 userpk/comm/unpack dtcreate dtcomm
 */
void printResults(dtypeTrials dtimes[])
{
    int wrank, tidx;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    if (wrank == 0) {
	double tucommmin, tucommavg, tcommmin, tcommavg;
	double tupkmin, tupkavg, tuupkmin, tuupkavg;
	double urate, dtrate;      /* In GB/sec */
	printf("Datatype\top type\tbytes\tUser    \tDatatype\tU Packmin\tU PackAvg\tU Op min\tU Op avg\tU unpk min\tU unpk avg\tDT min  \tDT avg\n");
	for (tidx=0; dtimes[tidx].dttype != DT_NONE; tidx++) {
	    const char *label = dtimes[tidx].label;

	    /* For each op, get data and print results */
	    /* Pack: user in comm, dt in dtcomm */
	    getMinAvgTimes(dtimes[tidx].pack.comm, dtimes[tidx].ncur+1,
			   &tucommmin, &tucommavg);
	    getMinAvgTimes(dtimes[tidx].pack.dtcomm, dtimes[tidx].ncur+1,
			   &tcommmin, &tcommavg);
	    urate  = dtimes[tidx].nbytes / (tucommmin) / 1e9;
	    dtrate = dtimes[tidx].nbytes / (tcommmin) / 1e9;
	    printf("%10s\t%.6s\t%ld\t%.2e\t%.2e\t        \t        \t%.2e\t%.2e\t        \t        \t%.2e\t%.2e\n",
		   label, "pack", (long)dtimes[tidx].nbytes,
		   urate, dtrate,
		   tucommmin, tucommavg, tcommmin, tcommavg);

	    getMinAvgTimes(dtimes[tidx].unpack.comm, dtimes[tidx].ncur+1,
			   &tucommmin, &tucommavg);
	    getMinAvgTimes(dtimes[tidx].unpack.dtcomm, dtimes[tidx].ncur+1,
			   &tcommmin, &tcommavg);
	    urate  = dtimes[tidx].nbytes / (tucommmin) / 1e9;
	    dtrate = dtimes[tidx].nbytes / (tcommmin) / 1e9;
	    printf("%10s\t%.6s\t%ld\t%.2e\t%.2e\t        \t        \t%.2e\t%.2e\t        \t        \t%.2e\t%.2e\n",
		   label, "unpack", (long)dtimes[tidx].nbytes,
		   urate, dtrate,
		   tucommmin, tucommavg, tcommmin, tcommavg);

	    getMinAvgTimes(dtimes[tidx].pt2pt.comm, dtimes[tidx].ncur+1,
			   &tucommmin, &tucommavg);
	    getMinAvgTimes(dtimes[tidx].pt2pt.pack, dtimes[tidx].ncur+1,
			   &tupkmin, &tupkavg);
	    getMinAvgTimes(dtimes[tidx].pt2pt.unpack, dtimes[tidx].ncur+1,
			   &tuupkmin, &tuupkavg);
	    getMinAvgTimes(dtimes[tidx].pt2pt.dtcomm, dtimes[tidx].ncur+1,
			   &tcommmin, &tcommavg);
	    urate  = dtimes[tidx].nbytes / (tupkmin + tucommmin + tuupkmin) / 1e9;
	    dtrate = dtimes[tidx].nbytes / (tcommmin) / 1e9;
	    printf("%10s\t%.6s\t%ld\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
		   label, "pt2pt", (long)dtimes[tidx].nbytes,
		   urate, dtrate,
		   tupkmin, tupkavg, tucommmin, tucommavg, tuupkmin, tuupkavg,
		   tcommmin, tcommavg);
	}
    }
}

void initTestTypes(dtypeTrials dtimes[], int dtimesSize, MPI_Datatype basetype)
{
    int typesize;
    int idx = 0;
    int count, blklen;
    /* Create the description of the tests to run */
    /* Eventually, use command line options and/or read a file.
       For example, options might be v:1000:1:1000 for a vector with
       count = 1000, blklen = 1, stride = 1000 */

    /* Create vectors with different lengths and strides */
    for (blklen = 1; blklen<5; blklen*=2) {
	for (count=1000; count<10000; count *=2) {
	    if (idx >= dtimesSize)
		myAbort(MPI_COMM_WORLD, 1, "Too many tests defined");
	    dtimes[idx].dttype   = DT_VECTOR;
	    dtimes[idx].dtcount  = 1;
	    dtimes[idx].parms[0] = count;     /* count */
	    dtimes[idx].parms[1] = blklen;    /* blklen */
	    dtimes[idx].parms[2] = count;     /* stride */
	    dtimes[idx].uparms[0] = count;
	    dtimes[idx].uparms[1] = blklen;
	    dtimes[idx].uparms[2] = count;
	    dtimes[idx].aparms   = 0;
	    dtimes[idx].basetype = basetype;
	    dtimes[idx].label    = "V-mm";
	    idx++;

	    if (verbose && blklen > 1)
		printf("count = %d, blklen = %d, count %% blklen = %d\n",
		       count, blklen, count % blklen);
	    if (blklen > 1 && (count % blklen) == 0) {
		if (idx >= dtimesSize)
		    myAbort(MPI_COMM_WORLD, 1, "Too many tests defined");
		/* Also build a version using a contiguous type */
		dtimes[idx].dttype    = DT_VECTOR;
		dtimes[idx].dtcount   = 1;
		dtimes[idx].parms[0]  = count;     /* count */
		dtimes[idx].parms[1]  = 1;         /* blklen */
		dtimes[idx].parms[2]  = count;     /* stride */
		dtimes[idx].uparms[0] = count;
		dtimes[idx].uparms[1] = blklen;
		dtimes[idx].uparms[2] = count;
		dtimes[idx].aparms   = 0;
		MPI_Type_contiguous(blklen, basetype, &dtimes[idx].basetype);
		dtimes[idx].label    = "Vb-mm";
		idx++;
	    }

	    /* Resized version */
	    if (blklen == 1) {
		if (idx >= dtimesSize)
		    myAbort(MPI_COMM_WORLD, 1, "Too many tests defined");
		MPI_Type_size(basetype, &typesize);
		dtimes[idx].dttype    = DT_RESIZED;
		dtimes[idx].dtcount   = count;
		dtimes[idx].parms[0]  = count*typesize;
		dtimes[idx].uparms[0] = count;
		dtimes[idx].uparms[1] = blklen;
		dtimes[idx].uparms[2] = count;
		dtimes[idx].aparms    = 0;
		dtimes[idx].basetype  = basetype;
		dtimes[idx].label     = "R-mm";
		idx++;
	    }
	    else {
		if (idx >= dtimesSize)
		    myAbort(MPI_COMM_WORLD, 1, "Too many tests defined");
		/* Create a contiguous type with blklen of the basetype */
		/* The stride is still in terms of the base type */
		MPI_Type_contiguous(blklen, basetype, &dtimes[idx].basetype);
		MPI_Type_size(basetype, &typesize);
		dtimes[idx].dttype    = DT_RESIZED;
		dtimes[idx].dtcount   = count;
		dtimes[idx].parms[0]  = count*typesize;
		dtimes[idx].uparms[0] = count;
		dtimes[idx].uparms[1] = blklen;
		dtimes[idx].uparms[2] = count;
		dtimes[idx].aparms    = 0;
		dtimes[idx].label     = "Rb-mm";
		idx++;
	    }
	}
    }

    if (idx >= dtimesSize)
	myAbort(MPI_COMM_WORLD, 1, "Too many tests defined");
    dtimes[idx].dttype   = DT_NONE;
}

void createDatatype(dtypeTrials *dtimes, MPI_Datatype basetype,
		    MPI_Datatype *dt)
{
    double t;
    t = MPI_Wtime();
    switch (dtimes->dttype) {
    case DT_CONTIG:
	t = MPI_Wtime();
	MPI_Type_contiguous(dtimes->parms[0], basetype, dt);
	MPI_Type_commit(dt);
	break;
    case DT_VECTOR:
	MPI_Type_vector(dtimes->parms[0], dtimes->parms[1],
			dtimes->parms[2], basetype, dt);
	MPI_Type_commit(dt);
	break;
    case DT_RESIZED:
	MPI_Type_create_resized(basetype, 0,
				(MPI_Aint)dtimes->parms[0], dt);
	MPI_Type_commit(dt);
	break;
    case DT_INDEX:
    case DT_BLOCKINDEX:
    case DT_STRUCT:
    default:
	myAbort(MPI_COMM_WORLD, 1, "Unsupported datatype for test");
	break;
    }
    dtimes->dtcreate = MPI_Wtime() - t;
}

/*
 * Add a 3D vector face type, with xy, yz, and xz options
 * Only really interesting if using a halo.
 */

/*
 * File input option:
 * type,label,parms [\]
 *
 * type:
 *
 */
