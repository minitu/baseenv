#include "baseenv.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* mputil allows us to time on any number of cores */
#include "mputil.h"

/* Define DO_HEAT_COMPUTATION to compute the sum of new values within the
   loop.  This typically affects the code generated */

/* */

#define ind(i,j,k) (k)*((by+2)*(bx+2))+(j)*(bx+2)+(i)

void sweep(int bx, int by, int bz, double *restrict anew,
	   const double *restrict aold, double *heatPtr);
void sweepb(int bx, int by, int bz, double *restrict anew,
	    const double *restrict aold, double *heatPtr);
void sweepi(int bx, int by, int bz, double *restrict anew,
	    const double *restrict aold, double *heatPtr);
void sweepbi(int bx, int by, int bz, double *restrict anew,
	     const double *restrict aold, double *heatPtr);


void sweep(int bx, int by, int bz, double *restrict anew,
	   const double *restrict aold, double *heatPtr)
{
    double heat;

    // update grid points
    heat = 0.0;
//#pragma @ICE loop=sweep
    for (int k=1; k<bz+1; ++k) {
	for(int j=1; j<by+1; ++j) {
	    for(int i=1; i<bx+1; ++i) {
		anew[ind(i,j,k)] = aold[ind(i,j,k)]/2.0 + (aold[ind(i-1,j,k)] + aold[ind(i+1,j,k)] + aold[ind(i,j-1,k)] + aold[ind(i,j+1,k)] + aold[ind(i,j,k-1)] + aold[ind(i,j,k+1)])/6.0/2.0;
#ifdef DO_HEAT_COMPUTATION
		heat += anew[ind(i,j,k)];
#endif
	    }
	}
    }

    *heatPtr = heat;
}

void sweepi(int bx, int by, int bz, double *restrict anew,
	    const double *restrict aold, double *heatPtr)
{
    double heat;
    const double *restrict ac;
    double *restrict an;
    int   yoffset = bx+2;
    int   zoffset = (bx+2)*(by+2);
    double s1 = 0.5, s2 = 1.0/(6.0*2.0);

    // update grid points
    heat = 0.0;
    for (int k=1; k<bz+1; ++k) {
	for(int j=1; j<by+1; ++j) {
	    ac = aold + ind(0,j,k);
	    an = anew + ind(0,j,k);
	    for(int i=1; i<bx+1; ++i) {
		an[i] = ac[i]*s1 + (ac[i-1] + ac[i+1] + ac[i-yoffset] + ac[i+yoffset] + ac[i-zoffset] + ac[i+zoffset])*s2;
#ifdef DO_HEAT_COMPUTATION
		heat += an[i];
#endif
	    }
	}
    }

    *heatPtr = heat;
}

#define BLOCKSIZE 8
void sweepb(int bx, int by, int bz, double *restrict anew,
	     const double *restrict aold, double *heatPtr)
{
    double heat;
    int jt, jmax;

    // update grid points
    heat = 0.0;
    for (jt=1; jt<by+1; jt+= BLOCKSIZE) {
	jmax = jt + BLOCKSIZE;
	if (jmax > by+1) jmax=by+1;
	for (int k=1; k<bz+1; ++k) {
	    for(int j=jt; j<jmax; ++j) {
		for(int i=1; i<bx+1; ++i) {
		    anew[ind(i,j,k)] = aold[ind(i,j,k)]/2.0 + (aold[ind(i-1,j,k)] + aold[ind(i+1,j,k)] + aold[ind(i,j-1,k)] + aold[ind(i,j+1,k)] + aold[ind(i,j,k-1)] + aold[ind(i,j,k+1)])/6.0/2.0;
#ifdef DO_HEAT_COMPUTATION
		heat += anew[ind(i,j,k)];
#endif
		}
	    }
	}
    }

    *heatPtr = heat;
}

void sweepbi(int bx, int by, int bz, double *restrict anew,
	     const double *restrict aold, double *heatPtr)
{
    double heat;
    const double *restrict ac;
    double *restrict an;
    int   yoffset = bx+2;
    int   zoffset = (bx+2)*(by+2);
    double s1 = 0.5, s2 = 1.0/(6.0*2.0);
    int jt, jmax;

    // update grid points
    heat = 0.0;
    for (jt=1; jt<by+1; jt+= BLOCKSIZE) {
	jmax = jt + BLOCKSIZE;
	if (jmax > by+1) jmax=by+1;
	for (int k=1; k<bz+1; ++k) {
	    for(int j=jt; j<jmax; ++j) {
		ac = aold + ind(0,j,k);
		an = anew + ind(0,j,k);
		for(int i=1; i<bx+1; ++i) {
		    an[i] = ac[i]*s1 + (ac[i-1] + ac[i+1] + ac[i-yoffset] + ac[i+yoffset] + ac[i-zoffset] + ac[i+zoffset])*s2;
#ifdef DO_HEAT_COMPUTATION
		    heat += an[i];
#endif
		}
	    }
	}
    }

    *heatPtr = heat;
}

int main(int argc, char **argv)
{
    int n, i;
    int k, ntest=16;
    double *aold, *anew, temp;
    clock_t t0, t1;
    double  t, rate, tmin;

    MPUTIL_INIT(0);

    //MPUTIL_SETDEBUG(1);

    if (argc != 2) {
	fprintf(stderr, "usage: %s n\n", argv[0]);
	MPUTIL_ABORT(1);
    }
    n  = atoi(argv[1]);

    if (n > 1024) {
	fprintf(stderr, "Maximum value of n is 1024\n");
	MPUTIL_ABORT(1);
    }
    aold = (double *)malloc(n*n*n*sizeof(double));
    anew = (double *)malloc(n*n*n*sizeof(double));

    if (!aold || !anew) {
	fprintf(stderr, "Unable to allocate %d words of storage\n", n*n*n);
	MPUTIL_ABORT(1);
    }

    for (i=0; i<n*n*n; i++) { aold[i] = i; anew[i] = -1.0; }

    MPUTIL_BEGIN;

    /* avoid cold start */
    sweep(n-2, n-2, n-2, anew, aold, &temp);
    temp = 0;

    MPUTIL_LABEL("Method");
    MPUTIL_OUTAPP("\tTime\tRate\n");

    MPUTIL_LABEL("Basic:  ");
    MPUTIL_SYNC;
    tmin = 1.e30;
    for (k=0; k<ntest; k++) {
	t0   = clock();
	sweep(n-2, n-2, n-2, anew, aold, &temp);
	sweep(n-2, n-2, n-2, aold, anew, &temp);
	t1   = clock();
	t    = (t1 - t0) * 1.e-6;
	MPUTIL_GETDMAX(t,t);
	if (t < tmin) tmin = t;
    }

    rate = 2*8*(n-2)*(n-2)*(n-2) / tmin;
    /* time then rate */
    MPUTIL_OUTAPP("\t%.2e\t%.2e\n", tmin, rate*1.e-6);

    /* avoid cold start */
    sweepi(n-2, n-2, n-2, anew, aold, &temp);
    temp = 0;

    MPUTIL_LABEL("BasicI:  ");
    MPUTIL_SYNC;
    tmin = 1.e30;
    for (k=0; k<ntest; k++) {
	t0   = clock();
	sweepi(n-2, n-2, n-2, anew, aold, &temp);
	sweepi(n-2, n-2, n-2, aold, anew, &temp);
	t1   = clock();
	t    = (t1 - t0) * 1.e-6;
	MPUTIL_GETDMAX(t,t);
	if (t < tmin) tmin = t;
    }

    rate = 2*8*(n-2)*(n-2)*(n-2) / tmin;
    /* time then rate */
    MPUTIL_OUTAPP("\t%.2e\t%.2e\n", tmin, rate*1.e-6);

    /* avoid cold start */
    sweepb(n-2, n-2, n-2, anew, aold, &temp);
    temp = 0;

    MPUTIL_LABEL("Blocked:");
    MPUTIL_SYNC;
    tmin = 1.e30;
    for (k=0; k<ntest; k++) {
	t0   = clock();
	sweepb(n-2, n-2, n-2, anew, aold, &temp);
	sweepb(n-2, n-2, n-2, aold, anew, &temp);
	t1   = clock();
	t    = (t1 - t0) * 1.e-6;
	MPUTIL_GETDMAX(t,t);
	if (t < tmin) tmin = t;
    }

    rate = 2*8*(n-2)*(n-2)*(n-2) / tmin;
    /* time then rate */
    MPUTIL_OUTAPP("\t%.2e\t%.2e\n", tmin, rate*1.e-6);

    /* avoid cold start */
    sweepbi(n-2, n-2, n-2, anew, aold, &temp);
    temp = 0;

    MPUTIL_LABEL("BlockedI:");
    MPUTIL_SYNC;
    tmin = 1.e30;
    for (k=0; k<ntest; k++) {
	t0   = clock();
	sweepbi(n-2, n-2, n-2, anew, aold, &temp);
	sweepbi(n-2, n-2, n-2, aold, anew, &temp);
	t1   = clock();
	t    = (t1 - t0) * 1.e-6;
	MPUTIL_GETDMAX(t,t);
	if (t < tmin) tmin = t;
    }

    rate = 2*8*(n-2)*(n-2)*(n-2) / tmin;
    /* time then rate */
    MPUTIL_OUTAPP("\t%.2e\t%.2e\n", tmin, rate*1.e-6);

    MPUTIL_END;

    free(aold);
    free(anew);

    MPUTIL_FINALIZE;

    return 0;
}
