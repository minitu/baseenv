#include "mputil.h"
#include <stdio.h>

/* Vector example code */
#include <time.h>

#define LENV 32000
#define LENV2 1000000
// OFFSET roughly sqrt(LENV2)
#define OFFSET 1000

__attribute__ ((aligned(16))) double X[LENV2], Y[LENV2];
void dummy(double *X, double *Y);

int main(int argc, char *argv[])
{
    int     ntimes=2000, ntimes2=200;
    clock_t start_t, end_t;
    double  t, c;

    MPUTIL_INIT(0);

    MPUTIL_BEGIN;

    for (int i=0; i<LENV2; i++) {
	X[i] = 1.0;
	Y[i] = 1.0;
    }
    c = argc + 1.3;

    MPUTIL_MASTER_BEGIN;
    printf("Next iteration\n");fflush(stdout);
    MPUTIL_MASTER_END;

    MPUTIL_LABEL("X=X+cY for\t%d", LENV);
    MPUTIL_SYNC;
    start_t = clock();
    for (int nl=0; nl<ntimes; nl++) {
	for (int i=0; i<LENV; i++) {
	    X[i] = X[i] + c * Y[i];
	}
	dummy(X, Y);
    }
    end_t = clock();
    t     = (end_t - start_t) / 1000000.0;
    // need to take min/max over participating processes
    MPUTIL_OUTAPP("\t%.2e\n", t/(ntimes));

    MPUTIL_LABEL("X=X+cY for %d", LENV2);
    MPUTIL_SYNC;
    start_t = clock();
    for (int nl=0; nl<ntimes2; nl++) {
	for (int i=0; i<LENV2; i++) {
	    X[i] = X[i] + c * Y[i];
	}
	dummy(X, Y);
    }
    end_t = clock();
    t     = (end_t - start_t) / 1000000.0;
    // need to take min/max over participating processes
    MPUTIL_OUTAPP("\t%.2e\n", t/(ntimes2));

    // Only the "interior" processed (see merged loop)
    MPUTIL_LABEL("Y=X--2X+X+ for %d", LENV2);
    MPUTIL_SYNC;
    start_t = clock();
    for (int nl=0; nl<ntimes2; nl++) {
	for (int i=OFFSET; i<LENV2-OFFSET; i++) {
	    Y[i] = X[i-1] - 2*X[i] + X[i+1];
	}
	dummy(X, Y);
    }
    end_t = clock();
    t     = (end_t - start_t) / 1000000.0;
    // need to take min/max over participating processes
    MPUTIL_OUTAPP("\t%.2e\n", t/(ntimes2));

    MPUTIL_LABEL("Y=X-n-2X+X+n for %d", LENV2);
    MPUTIL_SYNC;
    start_t = clock();
    for (int nl=0; nl<ntimes2; nl++) {
	for (int i=OFFSET; i<LENV2-OFFSET; i++) {
	    Y[i] = X[i-OFFSET] - 2*X[i] + X[i+OFFSET];
	}
	dummy(X, Y);
    }
    end_t = clock();
    t     = (end_t - start_t) / 1000000.0;
    // need to take min/max over participating processes
    MPUTIL_OUTAPP("\t%.2e\n", t/(ntimes2));

    MPUTIL_LABEL("Y=X-n-X--4X+X++X+n for %d", LENV2);
    MPUTIL_SYNC;
    start_t = clock();
    for (int nl=0; nl<ntimes2; nl++) {
	for (int i=OFFSET; i<LENV2-OFFSET; i++) {
	    Y[i] = X[i-OFFSET] + X[i-1] - 4*X[i] + X[i+1] + X[i+OFFSET];
	}
	dummy(X, Y);
    }
    end_t = clock();
    t     = (end_t - start_t) / 1000000.0;
    // need to take min/max over participating processes
    MPUTIL_OUTAPP("\t%.2e\n", t/(ntimes2));

    MPUTIL_END;

    MPUTIL_FINALIZE;

    return 0;
}
