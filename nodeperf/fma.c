#include "baseenv.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <stdlib.h>

/*
 * A simple code to compute some floating point performance values for
 * several different operations.  These complement the STREAM values.
 *
 * While this code does time the sparse matrix-vector product, as noted below,
 * these timings are too simple to be used for serious analysis.
 */

void initVars(int n, double *a, double *b, double *c);
double getmintime(struct timeval *tStart, struct timeval *tEnd, double told);

int main(int argc, char *argv[])
{
    double *a, *x, *y;
    int    i, k, n;
    struct timeval  tStart, tEnd;
    double tadd, tsadd, tfma, tsfma, tmult;

    /* Some arrays will be used; the sizes should be small enough to be
       within cache */
    n = 100;
    if (argc > 1) n = atoi(argv[1]);

    a = (double *)malloc(n*sizeof(double));
    x = (double *)malloc(n*sizeof(double));
    y = (double *)malloc(n*sizeof(double));

    initVars(n, a, x, y);

    /* Perform a matrix-vector multiply: y = A*x */
    /* Warning: To use this for timing, you need to (a) handle cold start
       (b) perform enough tests to make timing quantum relatively small */

    tadd  = 1.0e10;
    tsfma = 1.0e10;
    tfma  = 1.0e10;
    tsadd = 1.0e10;
    tmult = 1.0e10;
    for (k=0; k<10; k++) {
	double sum, sa, sb;
	int    kk;
	gettimeofday( &tStart, 0 );
	for (kk=0; kk<100; kk++) {
	    for (i=0; i<n; i++) {
		a[i] = a[i] + x[i];
	    }
	}
	gettimeofday( &tEnd, 0 );
	tadd = getmintime( &tStart, &tEnd, tadd);

	gettimeofday( &tStart, 0 );
	for (kk=0; kk<100; kk++) {
	    for (i=0; i<n; i++) {
		a[i] = a[i] + 1.0;
	    }
	}
	gettimeofday( &tEnd, 0 );
	tsadd = getmintime( &tStart, &tEnd, tsadd);

	sum = k;
	sa  = 1.0;
	sb  = 1.001;
	if (argc > 10) sa = -sa;
	gettimeofday( &tStart, 0 );
	for (kk=0; kk<100; kk++) {
	    for (i=0; i<n; i++) {
		sum = sum + sa*sb;
	    }
	}
	gettimeofday( &tEnd, 0 );
	tsfma = getmintime( &tStart, &tEnd, tsfma);

	gettimeofday( &tStart, 0 );
	for (kk=0; kk<100; kk++) {
	    for (i=0; i<n; i++) {
		a[i] = a[i] + x[i]*y[i];
	    }
	}
	gettimeofday( &tEnd, 0 );
	tfma = getmintime( &tStart, &tEnd, tfma);

	initVars(n, a, x, y);
	gettimeofday( &tStart, 0 );
	for (kk=0; kk<100; kk++) {
	    for (i=0; i<n; i++) {
		a[i] = a[i] * x[i];
	    }
	}
	gettimeofday( &tEnd, 0 );
	tmult = getmintime( &tStart, &tEnd, tmult);

    }
    /* Compute with the result to keep the compiler for marking the
       code as dead */
    for (i=0; i<n; i++) {
	if (a[i] < 0) {
	    fprintf(stderr,"a[%d]=%f, fails consistency test\n", i, a[i]);
	}
    }
    tadd  /= 100;
    tsadd /= 100;
    tfma  /= 100;
    tmult /= 100;
    printf("Op\tn\tTime     \tRate\n");
    printf("%s\t%d\t%.2e\t%.2e\n", "Add", n, tadd, n / tadd);
    printf("%s\t%d\t%.2e\t%.2e\n", "SAdd", n, tsadd, n / tsadd);
    printf("%s\t%d\t%.2e\t%.2e\n", "SFMA", n, tsfma, n / tsfma);
    printf("%s\t%d\t%.2e\t%.2e\n", "FMA", n, tfma, n / tfma);
    printf("%s\t%d\t%.2e\t%.2e\n", "Mult", n, tmult, n / tmult);

    free(a); free(x); free(y);
    return 0;
}

void initVars(int n, double *a, double *b, double *c)
{
    int i;
    for (i=0; i<n; i++) {
	a[i] = 1.0e-10;
	b[i] = 1.0e-10;
	c[i] = i*b[i];
    }
}

double getmintime(struct timeval *tStart, struct timeval *tEnd, double told)
{
    double t = (tEnd->tv_sec - tStart->tv_sec) +
	    1.0e-6 * (tEnd->tv_usec - tStart->tv_usec);
    if (t < told) return t;
    return told;
}
