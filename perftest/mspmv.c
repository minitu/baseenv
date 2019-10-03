#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>

#include "mpptestconf.h"

/*
 * The purpose of this program is to quickly measure the spmv performance
 * for a node, using separate processes for each core.  This will need to
 * be run as an MPI program, and with one process per core.
 */

typedef struct {
    int    ntest;
    double tmin, tmax, tsum;
} Time_t;
typedef struct {
    int    n;   /* Array size */
    int    nflop; /* Number of floating point operations */
    int    nc;  /* Number of processes (cores) */
    Time_t tspmv, rtspmv;
} Result_t;
typedef struct {
    int      nc, rsize;  /* Number of cores; size of results array */
    int      nr;         /* Number of used entries in results */
    MPI_Comm comm;       /* Comm for these cores/processes */
    Result_t *results;   /* Results for these operations */
} NProcResult_t;

#ifndef MAX_ROW
#define MAX_ROW 10000000
#endif
#ifndef MAX_NNZ
#define MAX_NNZ (5*MAX_ROW)
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

static STREAM_TYPE a[MAX_NNZ],
    x[MAX_ROW],
    y[MAX_ROW];
static int ia[MAX_ROW+1], ja[MAX_NNZ];

static int debug = 0;
enum { DEBUG_BASE = 1, DEBUG_ALL=10 };

void clearResults(int n, int nc, Result_t *results);
void initMatrix(int n, STREAM_TYPE *a, int *ia, int *ja, STREAM_TYPE *x,
		STREAM_TYPE *y);

void rspmv(int nr, const STREAM_TYPE *restrict a, const int *restrict ia,
	   const int *restrict ja, const STREAM_TYPE *restrict x,
	   STREAM_TYPE *restrict y);
void dummy(int n, STREAM_TYPE *a, STREAM_TYPE *x, STREAM_TYPE *y);

int main(int argc, char *argv[])
{
    int            provided, wsize, wrank, j, k, nr;
    STREAM_TYPE    scalar;
    int            n, nc;
    NProcResult_t Timings[MAX_PROCS];
    Result_t      *results=0;
    MPI_Comm      newcomm;
    double        t;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Each process performs spmv tests */
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
	n = 1024;
	while (n <= MAX_ROW) {
	    /* Initialize the next test */
	    initMatrix(n, a, ia, ja, x, y);
	    clearResults(n,nc,results+nr);

	    /* Run the tests */
	    for (k=0; k<NTIMES; k++) {
		STREAM_TYPE sum;
		int         idx;
		if (debug >= DEBUG_ALL && wrank == 0) printf("tspmv:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		for (j=0; j<n; j++) {
		    sum = 0;
		    for (idx=ia[j];idx<ia[j+1];idx++)
			sum += a[idx]*x[ja[idx]];
		    y[j] = sum;
		}
		t = MPI_Wtime() - t;
		results[nr].tspmv.tsum += t;
		results[nr].nflop = 2*ia[n];
		if (t > 0 && t < results[nr].tspmv.tmin) results[nr].tspmv.tmin = t;
		if (t > results[nr].tspmv.tmax) results[nr].tspmv.tmax = t;

		dummy(n, a, x, y);

		if (debug >= DEBUG_ALL && wrank == 0) printf("rspmv:%d\n", n);
		MPI_Barrier(newcomm);
		t = MPI_Wtime();
		rspmv(n, a, ia, ja, x, y);
		t = MPI_Wtime() - t;
		results[nr].rtspmv.tsum += t;
		results[nr].nflop = 2*ia[n];  /* Yes, same for both */
		if (t > 0 && t < results[nr].rtspmv.tmin) results[nr].rtspmv.tmin = t;
		if (t > results[nr].rtspmv.tmax) results[nr].rtspmv.tmax = t;

	    }
	    n *= 4; /* Keep the 2d matrix square */
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
	    printf("n\tspmv    \trspmv\n");
	}
	if (wrank < nc) {
	    nr      = Timings[nc-1].nr;
	    newcomm = Timings[nc-1].comm;
	    results = Timings[nc-1].results;
	    for (k=0; k<nr; k++) {
		double times[2];
		times[0] = results[k].tspmv.tmin;
		times[1] = results[k].rtspmv.tmin;
		MPI_Allreduce(MPI_IN_PLACE, times, 2, MPI_DOUBLE,
			      MPI_MAX, newcomm);
		if (wrank == 0) {
		    printf("%d\t%.2e\t%.2e\n",
			   results[k].n, times[0], times[1]);
		}
	    }
	    /* Rates in MF/sec */
	    if (wrank == 0) {
		printf("Rates For %d processes\n", nc);
		printf("n\tspmv    \trspmv\n");
	    }
	    for (k=0; k<nr; k++) {
		double rates[2];
		double sz;
		/* Rates in MFLOP/sec (not MiFLOP/sec) */
		sz       = results[k].nflop * 1.0e-6;
		rates[0] = sz / results[k].tspmv.tmin;
		rates[1] = sz / results[k].rtspmv.tmin;
		MPI_Allreduce(MPI_IN_PLACE, rates, 2, MPI_DOUBLE,
			      MPI_MIN, newcomm);
		if (wrank == 0) {
		    printf("%d\t%.2e\t%.2e\n",
			   results[k].n, rates[0], rates[1]);
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
    results->tspmv.tmin    = 1.0e12;
    results->tspmv.tmax    = 0;
    results->tspmv.tsum    = 0;
    results->tspmv.ntest   = NTIMES;
    results->rtspmv.tmin   = 1.0e12;
    results->rtspmv.tmax   = 0;
    results->rtspmv.tsum   = 0;
    results->rtspmv.ntest  = NTIMES;
}

void initMatrix(int n, STREAM_TYPE *a, int *ia, int *ja, STREAM_TYPE *x,
		STREAM_TYPE *y)
{
    /* Create a pentadiagonal matrix, representing very roughly a finite
       difference approximation to the Laplacian on a square n x n mesh */
    int meshn = sqrt((double)n);
    int i, j, row, nnz;

    row = 0;
    nnz = 0;
    for (i=0; i<meshn; i++) {
	for (j=0; j<meshn; j++) {
	    ia[row] = nnz;
	    if (i>0) { ja[nnz] = row - meshn; a[nnz] = -1.0; nnz++; }
	    if (j>0) { ja[nnz] = row - 1; a[nnz] = -1.0; nnz++; }
	    ja[nnz] = row; a[nnz] = 4.0; nnz++;
	    if (j<meshn-1) { ja[nnz] = row + 1; a[nnz] = -1.0; nnz++; }
	    if (i<meshn-1) { ja[nnz] = row + meshn; a[nnz] = -1.0; nnz++; }
	    row++;
	}
    }
    ia[row] = nnz;

    /* Create the source (x) vector and clear the dest vector */
    for (i=0; i<n; i++) {
	x[i] = 1.0;
	y[i] = 0.0;
    }

}
