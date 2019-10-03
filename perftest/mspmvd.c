#include "mpptestconf.h"

#ifndef STREAM_TYPE
#define STREAM_TYPE double
#endif
void dummy(int n, STREAM_TYPE *a, STREAM_TYPE *b, STREAM_TYPE *c);
void rspmv(int nr, const STREAM_TYPE *restrict a, const int *restrict ia,
	   const int *restrict ja, const STREAM_TYPE *restrict x,
	   STREAM_TYPE *restrict y);

void rspmv(int nr, const STREAM_TYPE *restrict a, const int *restrict ia,
	   const int *restrict ja, const STREAM_TYPE *restrict x,
	   STREAM_TYPE *restrict y)
{
    int j, idx;
    STREAM_TYPE sum;
    for (j=0; j<nr; j++) {
	sum = 0;
	for (idx=ia[j];idx<ia[j+1];idx++)
	    sum += a[idx]*x[ja[idx]];
	y[j] = sum;
    }
}
void dummy(int n, STREAM_TYPE *a, STREAM_TYPE *b, STREAM_TYPE *c){}
