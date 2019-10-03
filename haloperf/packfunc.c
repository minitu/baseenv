/* Routines to perform manual back/unpack */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "haloperf.h"
#include "packfunc.h"

#define NTEST 100
/* All routines return the time taken by the pack operation */

double pksimple( double *inbuf, double *outbuf, int count, int stride )
{
    int i, k;
    double t1, *p;

    t1 = MPI_Wtime();
    for (k=0; k<NTEST; k++) {
	p = inbuf;
	for (i=0; i<count; i++) {
	    outbuf[i] = *p;
	    p += stride;
	}
    }
    t1 = MPI_Wtime() - t1;
    return t1 / NTEST;
}

double pkconst( const double *inbuf, double *outbuf, int count, int stride )
{
    int i, k;
    double t1;
    const double *p;

    t1 = MPI_Wtime();
    for (k=0; k<NTEST; k++) {
	p = inbuf;
	for (i=0; i<count; i++) {
	    outbuf[i] = *p;
	    p += stride;
	}
    }
    t1 = MPI_Wtime() - t1;
    return t1 / NTEST;
}

double pkrestrict( const double *restrict inbuf, double *restrict outbuf, 
		   int count, int stride )
{
    int i, k;
    double t1;

    t1 = MPI_Wtime();
    for (k=0; k<NTEST; k++) {
	const double *restrict p = inbuf;
	for (i=0; i<count; i++) {
	    outbuf[i] = *p;
	    p += stride;
	}
    }
    t1 = MPI_Wtime() - t1;
    return t1 / NTEST;
}

double pksplit( const double *restrict inbuf, double *restrict outbuf, 
		int count, int stride )
{
    int i, k;
    double t1;
    register double tmp1, tmp2;
    int      cby2 = count/2;

    t1 = MPI_Wtime();
    for (k=0; k<NTEST; k++) {
	const double *restrict p1;
	const double *restrict p2;
	double *restrict o2;
	p1 = inbuf;
	p2 = inbuf + cby2 * stride;
	o2 = outbuf + cby2;
	for (i=0; i<cby2; i++) {
	    tmp1 = *p1;
	    tmp2 = *p2;
	    outbuf[i] = tmp1;
	    *o2++  = tmp2;
	    p1 += stride;
	    p2 += stride;
	}
    }
    t1 = MPI_Wtime() - t1;
    return t1 / NTEST;
}



