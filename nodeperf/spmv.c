#include "baseenv.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <stdlib.h>

void rspmv( int nrows, const int *restrict ia, const int * restrict ja,
	    const double * restrict a, const double *restrict x,
	    double *restrict y);
void rspmv2( int nrows, const int *restrict ia, const int * restrict ja,
	    const double * restrict a, const double *restrict x,
	    double *restrict y);

/*D
  spmv-c - Report the performance of an example sparse matrix-vector multiply

  Input parameter:
.  n - The mesh size for the matrix.  The number of rows in the matrix is
  n*n

Notes:
 This routine is an example of sparse matrix-vector multiply, using
 CSR (compressed sparse row format).  Because this code is in C, the
 index values in the 'ia' and 'ja' arrays are zero origin rather than one
 origin as they would be in Fortran.

 While this code does time the sparse matrix-vector product,
 these timings are too simple to be used for serious analysis.

  D*/
int main(int argc, char *argv[])
{
    int    *ia, *ja;
    double *a, *x, *y;
    int    row, i, j, k, idx, n, nnzMax, nnz, nrows;
    struct timeval  tStart, tEnd;
    double t, tmin, trmin, tr2min;

    n = 10;
    if (argc > 1) n = atoi(argv[1]);

    /* Warning: if n > sqrt(2^31), you may get integer overflow */

    /* Allocate enough storage for the matrix.  We allocate more than
       is needed in order to simplify the code */
    nrows  = n * n;
    nnzMax = nrows * 5;
    ia = (int*)malloc(nrows*sizeof(int));
    ja = (int*)malloc(nnzMax*sizeof(int));
    a  = (double*)malloc(nnzMax*sizeof(double));
    /* Allocate the source and result vectors */
    x = (double*)malloc(nrows*sizeof(double));
    y = (double*)malloc(nrows*sizeof(double));

    /* Create a pentadiagonal matrix, representing very roughly a finite
       difference approximation to the Laplacian on a square n x n mesh */
    row = 0;
    nnz = 0;
    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    ia[row] = nnz;
	    /* Put the diagonal first.  Some versions will  take advantage
	       of this */
	    ja[nnz] = row; a[nnz] = 4.0; nnz++;
	    if (i>0) { ja[nnz] = row - n; a[nnz] = -1.0; nnz++; }
	    if (j>0) { ja[nnz] = row - 1; a[nnz] = -1.0; nnz++; }
	    if (j<n-1) { ja[nnz] = row + 1; a[nnz] = -1.0; nnz++; }
	    if (i<n-1) { ja[nnz] = row + n; a[nnz] = -1.0; nnz++; }
	    row++;
	}
    }
    ia[row] = nnz;

    /* Create the source (x) vector */
    for (i=0; i<nrows; i++) x[i] = 1.0;

    /* Perform a matrix-vector multiply: y = A*x */
    /* Warning: To use this for timing, you need to (a) handle cold start
       (b) perform enough tests to make timing quantum relatively small */

    tmin   = 1.0e10;
    trmin  = 1.0e10;
    tr2min = 1.0e10;
    for (k=0; k<10; k++) {
	gettimeofday( &tStart, 0 );
	for (row=0; row<nrows; row++) {
	    double sum = 0.0;
	    for (idx=ia[row]; idx<ia[row+1]; idx++) {
		sum += a[idx] * x[ja[idx]];
	    }
	    y[row] = sum;
	}
	gettimeofday( &tEnd, 0 );
	t = (tEnd.tv_sec - tStart.tv_sec) +
	    1.0e-6 * (tEnd.tv_usec - tStart.tv_usec);
	if (t < tmin) tmin = t;

	gettimeofday( &tStart, 0 );
	rspmv(nrows,ia,ja,a,x,y);
	gettimeofday( &tEnd, 0 );
	t = (tEnd.tv_sec - tStart.tv_sec) +
	    1.0e-6 * (tEnd.tv_usec - tStart.tv_usec);
	if (t < trmin) trmin = t;

	gettimeofday( &tStart, 0 );
	rspmv2(nrows,ia,ja,a,x,y);
	gettimeofday( &tEnd, 0 );
	t = (tEnd.tv_sec - tStart.tv_sec) +
	    1.0e-6 * (tEnd.tv_usec - tStart.tv_usec);
	if (t < tr2min) tr2min = t;

	/* Compute with the result to keep the compiler for marking the
	   code as dead */
	for (row=0; row<nrows; row++) {
	    if (y[row] < 0) {
		fprintf(stderr,"y[%d]=%f, fails consistency test\n", row, y[row]);
	    }
	}
    }

    printf("Time for Sparse Ax\n");
    printf("Op                  \tNrows\tnnz\tTime\n");
    printf("Inline              \t%d\t%d\t%.2e\n", nrows, nnz, tmin);
    printf("Routine             \t%d\t%d\t%.2e\n", nrows, nnz, trmin);
    printf("Routine, unroll by 4\t%d\t%d\t%.2e\n", nrows, nnz, tr2min);

    free(ia); free(ja); free(a); free(x); free(y);
    return 0;
}

void rspmv( int nrows, const int *restrict ia, const int * restrict ja,
	    const double * restrict a, const double *restrict x,
	    double *restrict y)
{
    int row, idx;

    for (row=0; row<nrows; row++) {
	double sum = 0.0;
	for (idx=ia[row]; idx<ia[row+1]; idx++) {
	    sum += a[idx] * x[ja[idx]];
	}
	y[row] = sum;
    }
}

void rspmv2( int nrows, const int *restrict ia, const int * restrict ja,
	     const double * restrict a, const double *restrict x,
	     double *restrict y)
{
    int row, idx, in, in2;

    idx = ia[0];
    in  = ia[1];
    for (row=0; row<nrows; row++) {
	double sum = a[idx++]*x[row];
	in2 = ia[row+2];
	if (in - idx == 4) {
	    int i0 = ja[idx], i1 = ja[idx+1], i2 = ja[idx+2], i3=ja[idx+3];
//	    double x0=x[i0], x1=x[i1],x2=x[i2], x3=x[i3];
//	    sum += (a[idx]*x0+a[idx+1]*x1)+(a[idx+2]*x2+a[idx+3]*x3);
	    double xx[4]={x[i0],x[i1],x[i2],x[i3]};
	    sum += a[idx]*xx[0]+a[idx+1]*xx[1]+a[idx+2]*xx[2]+a[idx+3]*xx[3];
	    idx=in;
	}
	else {
	    while (idx<in) {
		sum += a[idx] * x[ja[idx]]; idx++;
	    }
	}
	y[row] = sum;
	in     = in2;
    }
}
