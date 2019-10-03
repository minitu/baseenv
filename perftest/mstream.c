#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>

/* mpptestconf.h is needed for restrict, for poor or ancient C compilers
   that don't implement modern C */
#include "mpptestconf.h"

/*
 * The purpose of this program is to quickly measure the STREAM performance
 * for a node, using separate processes for each core.  This will need to
 * be run as an MPI program, and with one process per core.
 *
 * There are several differences with the basic STREAM tests.  These are:
 *    Runs for different sizes (upto the STREAM_ARRAY_SIZE)
 *    Runs using routines as well as main program
 *    Runs using dynamically allocated memory
 */

typedef struct {
    int    ntest;
    double tmin, tmax, tsum;
} Time_t;
typedef struct {
    int    n;   /* Array size */
    int    nc;  /* Number of processes (cores) */
    Time_t tcopy, tscale, tadd, ttriad, rtcopy, rtscale, rtadd, rttriad;
} Result_t;
typedef struct {
    int      nc, rsize;  /* Number of cores; size of results array */
    int      nr;         /* Number of used entries in results */
    MPI_Comm comm;       /* Comm for these cores/processes */
    Result_t *results;   /* Results for these operations */
} NProcResult_t;

#ifndef STREAM_ARRAY_SIZE
#define STREAM_ARRAY_SIZE 10000000
#endif
#ifndef STREAM_TYPE
#define STREAM_TYPE double
#endif
#ifndef NTIMES
#define NTIMES 10
#endif
#ifndef MAX_PROCS
#define MAX_PROCS 32
#endif

static STREAM_TYPE a[STREAM_ARRAY_SIZE],
    b[STREAM_ARRAY_SIZE],
    c[STREAM_ARRAY_SIZE];

static int debug = 0;
enum { DEBUG_BASE = 1, DEBUG_ALL=10 };

void clearResults(int n, int nc, Result_t *results);

void rcopy(int n, const STREAM_TYPE *restrict a, STREAM_TYPE *restrict b);
void rscale(int n, const STREAM_TYPE *restrict a, STREAM_TYPE *restrict b,
	    STREAM_TYPE scalar);
void radd(int n, const STREAM_TYPE *restrict a, const STREAM_TYPE *restrict b,
	  STREAM_TYPE *restrict c);
void rtriad(int n, const STREAM_TYPE *restrict a, const STREAM_TYPE *restrict b,
	    STREAM_TYPE scalar, STREAM_TYPE *restrict c);
void dummy(int n, STREAM_TYPE *a, STREAM_TYPE *b, STREAM_TYPE *c);

int main(int argc, char *argv[])
{
    int            provided, wsize, wrank, j, k, nr;
    STREAM_TYPE    *ap, *bp, *cp;
    STREAM_TYPE    scalar;
    int            n, nc;
    NProcResult_t Timings[MAX_PROCS];
    Result_t      *results=0;
    MPI_Comm      newcomm;
    double        t;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Each process performs stream tests */
    scalar = 1.00001;

    /* For each number of processes (cores), starting at 1, to wsize */
    for (nc=1; nc<=wsize; nc++) {
	if (debug && wrank == 0) printf("nc = %d\n", nc);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_split(MPI_COMM_WORLD, wrank < nc, wrank, &Timings[nc-1].comm);
	if (debug && wrank == 0) printf("split complete\n");
	if (wrank >= nc) {
	    continue;
	}
	Timings[nc-1].nc    = nc;
	Timings[nc-1].rsize = 16;
	Timings[nc-1].results =
	    (Result_t *)malloc(Timings[nc-1].rsize*sizeof(Result_t));
	/* Temps to make code clearer */
	newcomm = Timings[nc-1].comm;
	results = Timings[nc-1].results;
	nr = 0;

	if (debug) {
	    int ssize;
	    MPI_Comm_size(newcomm, &ssize);
	    printf( "nc = %d, ssize = %d\n", nc, ssize);
	    fflush(stdout);
	    if (ssize != nc) {
		fprintf(stderr,"Inconsistent size for newcomm: %d\n", ssize);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}

	/* Small values of n need a different test setup because the
	   amount of time is so small */
	n = 2048;
	while (n <= STREAM_ARRAY_SIZE) {
	    /* Initialize the next test */
	    ap = (STREAM_TYPE *)malloc(n*sizeof(STREAM_TYPE));
	    bp = (STREAM_TYPE *)malloc(n*sizeof(STREAM_TYPE));
	    cp = (STREAM_TYPE *)malloc(n*sizeof(STREAM_TYPE));
	    if (!ap || !bp || !cp) {
		fprintf(stderr,"Unable to allocated %d words for ap,bp,cp\n",
			n);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	    for (j=0; j<n; j++) {
		a[j] = 1.0;
		b[j] = 2.0;
		c[j] = 0.0;
		ap[j] = 1.0;
		bp[j] = 2.0;
		cp[j] = 0.0;
	    }
	    clearResults(n,nc,results+nr);

	    /* Run the tests */
	    for (k=0; k<NTIMES; k++) {
		if (debug >= DEBUG_ALL && wrank == 0) printf("tcopy:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		for (j=0; j<n; j++) {
		    c[j] = a[j];
		}
		t = MPI_Wtime() - t;
		results[nr].tcopy.tsum += t;
		if (t > 0 && t < results[nr].tcopy.tmin) results[nr].tcopy.tmin = t;
		if (t > results[nr].tcopy.tmax) results[nr].tcopy.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("tscale:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		for (j=0; j<n; j++) {
		    b[j] = scalar * c[j];
		}
		t = MPI_Wtime() - t;
		results[nr].tscale.tsum += t;
		if (t > 0 && t < results[nr].tscale.tmin) results[nr].tscale.tmin = t;
		if (t > results[nr].tscale.tmax) results[nr].tscale.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("tadd:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		for (j=0; j<n; j++) {
		    c[j] = a[j] + b[j];
		}
		t = MPI_Wtime() - t;
		results[nr].tadd.tsum += t;
		if (t > 0 && t < results[nr].tadd.tmin) results[nr].tadd.tmin = t;
		if (t > results[nr].tadd.tmax) results[nr].tadd.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("ttriad:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		for (j=0; j<n; j++) {
		    a[j] = b[j]+scalar*c[j];
		}
		t = MPI_Wtime() - t;
		results[nr].ttriad.tsum += t;
		if (t > 0 && t < results[nr].ttriad.tmin) results[nr].ttriad.tmin = t;
		if (t > results[nr].ttriad.tmax) results[nr].ttriad.tmax = t;

		dummy(n, a, b, c);

		if (debug >= DEBUG_ALL && wrank == 0) printf("rcopy:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		rcopy(n, c, a);
		t = MPI_Wtime() - t;
		results[nr].rtcopy.tsum += t;
		if (t > 0 && t < results[nr].rtcopy.tmin) results[nr].rtcopy.tmin = t;
		if (t > results[nr].rtcopy.tmax) results[nr].rtcopy.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("rscale:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		rscale(n, b, c, scalar);
		t = MPI_Wtime() - t;
		results[nr].rtscale.tsum += t;
		if (t > 0 && t < results[nr].rtscale.tmin) results[nr].rtscale.tmin = t;
		if (t > results[nr].rtscale.tmax) results[nr].rtscale.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("radd:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		radd(n, a, b, c);
		t = MPI_Wtime() - t;
		results[nr].rtadd.tsum += t;
		if (t > 0 && t < results[nr].rtadd.tmin) results[nr].rtadd.tmin = t;
		if (t > results[nr].rtadd.tmax) results[nr].rtadd.tmax = t;

		if (debug >= DEBUG_ALL && wrank == 0) printf("rtriad:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		rtriad(n, b, c, scalar, a);
		t = MPI_Wtime() - t;
		results[nr].rttriad.tsum += t;
		if (t > 0 && t < results[nr].rttriad.tmin) results[nr].rttriad.tmin = t;
		if (t > results[nr].rttriad.tmax) results[nr].rttriad.tmax = t;
	    }
	    free(ap);
	    free(bp);
	    free(cp);
	    n *= 2;
	    nr++;
	}
	Timings[nc-1].nr = nr;
    }

    /* Output results */
    if (debug) printf("Output results\n");

    for (nc=1; nc<=wsize; nc++) {
	/* Times */
	if (wrank == 0) {
	    printf("Times For %d processes\n", nc);
	    printf("n\tcopy    \tscale   \tadd     \ttriad   \trcopy   \trscale  \tradd    \trtriad\n");
	}
	if (wrank < nc) {
	    nr      = Timings[nc-1].nr;
	    newcomm = Timings[nc-1].comm;
	    results = Timings[nc-1].results;
	    for (k=0; k<nr; k++) {
		double times[8];
		times[0] = results[k].tcopy.tmin;
		times[1] = results[k].tscale.tmin;
		times[2] = results[k].tadd.tmin;
		times[3] = results[k].ttriad.tmin;
		times[4] = results[k].rtcopy.tmin;
		times[5] = results[k].rtscale.tmin;
		times[6] = results[k].rtadd.tmin;
		times[7] = results[k].rttriad.tmin;
		MPI_Allreduce(MPI_IN_PLACE, times, 8, MPI_DOUBLE,
			      MPI_MAX, newcomm);
		if (wrank == 0) {
		    printf("%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			   results[k].n, times[0], times[1], times[2], times[3],
			   times[4], times[5], times[6], times[7]);
		}
	    }
	    /* Rates */
	    if (wrank == 0) {
		printf("Rates For %d processes\n", nc);
		printf("n\tcopy    \tscale   \tadd     \ttriad   \trcopy   \trscale  \tradd    \trtriad\n");
	    }
	    for (k=0; k<nr; k++) {
		double rates[8];
		double sz;
		/* Rates in MB/sec (not MiB/sec) */
		sz       = results[k].n * sizeof(STREAM_TYPE) * 1.0e-6;
		rates[0] = sz / results[k].tcopy.tmin;
		rates[1] = sz / results[k].tscale.tmin;
		rates[2] = sz / results[k].tadd.tmin;
		rates[3] = sz / results[k].ttriad.tmin;
		rates[4] = sz / results[k].rtcopy.tmin;
		rates[5] = sz / results[k].rtscale.tmin;
		rates[6] = sz / results[k].rtadd.tmin;
		rates[7] = sz / results[k].rttriad.tmin;
		MPI_Allreduce(MPI_IN_PLACE, rates, 8, MPI_DOUBLE,
			      MPI_MIN, newcomm);
		if (wrank == 0) {
		    printf("%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			   results[k].n, rates[0], rates[1], rates[2], rates[3],
			   rates[4], rates[5], rates[6], rates[7]);
		}
	    }
	}
    }
    if (wrank == 0) fflush(stdout);

    /* Clean up */
    if (debug && wrank == 0) printf("Cleaning up...\n");
    for (nc=0; nc<wsize; nc++) {
	if (wrank < nc) {
	    free(Timings[nc].results);
	}
	MPI_Comm_free(&Timings[nc].comm);
    }

    MPI_Finalize();
    return 0;
}

void clearResults(int n, int nc, Result_t *results)
{
    results->n             = n;
    results->nc            = nc;
    results->tcopy.tmin    = 1.0e12;
    results->tcopy.tmax    = 0;
    results->tcopy.tsum    = 0;
    results->tcopy.ntest   = NTIMES;
    results->tscale.tmin   = 1.0e12;
    results->tscale.tmax   = 0;
    results->tscale.tsum   = 0;
    results->tscale.ntest  = NTIMES;
    results->tadd.tmin     = 1.0e12;
    results->tadd.tmax     = 0;
    results->tadd.tsum     = 0;
    results->tadd.ntest    = NTIMES;
    results->ttriad.tmin   = 1.0e12;
    results->ttriad.tmax   = 0;
    results->ttriad.tsum   = 0;
    results->ttriad.ntest  = NTIMES;
    results->rtcopy.tmin   = 1.0e12;
    results->rtcopy.tmax   = 0;
    results->rtcopy.tsum   = 0;
    results->rtcopy.ntest  = NTIMES;
    results->rtscale.tmin  = 1.0e12;
    results->rtscale.tmax  = 0;
    results->rtscale.tsum  = 0;
    results->rtscale.ntest = NTIMES;
    results->rtadd.tmin    = 1.0e12;
    results->rtadd.tmax    = 0;
    results->rtadd.tsum    = 0;
    results->rtadd.ntest   = NTIMES;
    results->rttriad.tmin  = 1.0e12;
    results->rttriad.tmax  = 0;
    results->rttriad.tsum  = 0;
    results->rttriad.ntest = NTIMES;
}
