#include "mpptestconf.h"

#ifndef STREAM_TYPE
#define STREAM_TYPE double
#endif

void rcopy(int n, const STREAM_TYPE *restrict a, STREAM_TYPE *restrict b) {
    int i;
    for (i=0; i<n; i++)
	b[i] = a[i];
}
void rscale(int n, const STREAM_TYPE *restrict a, STREAM_TYPE *restrict b,
	    STREAM_TYPE scalar)
{
    int i;
    for (i=0; i<n; i++)
	b[i] = scalar * a[i];
}
void radd(int n, const STREAM_TYPE *restrict a, const STREAM_TYPE *restrict b,
	  STREAM_TYPE *restrict c)
{
    int i;
    for (i=0; i<n; i++)
	c[i] = a[i] + b[i];
}
void rtriad(int n, const STREAM_TYPE *restrict a, const STREAM_TYPE *restrict b,
	    STREAM_TYPE scalar, STREAM_TYPE *restrict c)
{
    int i;
    for (i=0; i<n; i++)
	c[i] = a[i] + scalar*b[i];
}
void dummy(int n, STREAM_TYPE *a, STREAM_TYPE *b, STREAM_TYPE *c){}
