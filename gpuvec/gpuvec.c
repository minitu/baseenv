#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>   /* For pow */

#include "gpuvec.h"
FILE *outfp=0;

#ifdef USE_DYNAMIC_ALLOC
__attribute__((aligned(16))) BASETYPE *restrict a, *restrict b, *restrict c,
    *restrict d, *restrict e;
/* 2D arrays need to be statically allocated */
__attribute__((aligned(16))) BASETYPE
    aa[LEN2][LEN2],bb[LEN2][LEN2],cc[LEN2][LEN2],tt[LEN2][LEN2];
#if 0
__attribute__((aligned(16))) BASETYPE *restrict*restrict aa,
    *restrict*restrict bb,
    *restrict*restrict cc,
    *restrict*restrict tt;
#endif
int memalignsize=16;
/* Note that OpenACC doesn't permit static or global variables in accellerator
   code sections. This global will be communicated to the individual kernels */
int veclen=LEN;
#undef LEN
#define SETVECLEN int LEN=veclen;
#else
__attribute__((aligned(16))) BASETYPE a[LEN],b[LEN],c[LEN],d[LEN],e[LEN],
    aa[LEN2][LEN2],bb[LEN2][LEN2],cc[LEN2][LEN2],tt[LEN2][LEN2];
#define SETVECLEN
#endif

#ifdef RUNTIME_NTIMES
int ntimes = 200000;
#else
#ifndef ntimes
/* The default value is 200000 */
#define ntimes 200000
#endif
#endif

/* External routines */
/* PGI seems to want the acc routine declarations together, and before other
   prototypes */
#pragma acc routine (dummy) seq
#pragma acc routine (dummyab) seq
#pragma acc routine (dummyabc) seq
/* nohost can't be used it the target is the host processor - this is probably
   a design error in OpenACC */
/* I've also had trouble with the PGI compiler even when targeting the GPU.
   For compilers that cannot handle nohost properly, define IGNORE_NOHOST */
#if defined(DEVICE_IS_HOST) || defined(IGNORE_NOHOST)
#pragma acc routine (set1d) vector
#else
#pragma acc routine (set1d) vector nohost
#endif

/* Query: could we make these BASETYPE restrict * ? */
/* Query: The 2D arrays are a problem, particularly if not 2D allocation.
   Consider a subset of tests that only uses 1D allocations, but may
   use 2D operations with explicit index computation */
int dummy(BASETYPE *restrict a, BASETYPE *restrict b, BASETYPE *restrict c,
	  BASETYPE *restrict d, BASETYPE *restrict e,
	  BASETYPE aa[LEN2][LEN2], BASETYPE bb[LEN2][LEN2],
	  BASETYPE cc[LEN2][LEN2], BASETYPE s);

int dummyab(BASETYPE *restrict a, BASETYPE *restrict b);
int dummyabc(BASETYPE *restrict a, BASETYPE *restrict b, BASETYPE *restrict c);

void printSummary(const char *name, double t, long ops, BASETYPE checkval,
		  const char *shortdesc);
void getCheck(BASETYPE *a, int len, BASETYPE *checkval);
int set1d(BASETYPE *restrict a, int len, BASETYPE val);

int set1d(BASETYPE *restrict a, int len, BASETYPE val)
{
#pragma acc loop
    for (int i=0; i<len; i++) a[i] = val;
    return 0;
}

int vpvts(BASETYPE);
int vsumr(void);

int paxpy(BASETYPE);
int paxpby(BASETYPE, BASETYPE);
int paxpbypcz(BASETYPE, BASETYPE, BASETYPE);
int pmdot(void);

int stencil3d(int,int,int);


int main(int argc, char **argv)
{
    double  s = 3.14159, s1 = 0.291735, s2 = 0.3456;
    int i;
    int bx, by, bz;
    double l3;
    SETVECLEN;

    outfp = stdout;
    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-outfile") == 0) {
	    outfp = fopen(argv[i+1],"w");
	    if (!outfp) {
		fprintf(stderr, "Unable to open %s for writing\n", argv[i+1]);
		exit(1);
	    }
	    i++;
	}
#ifdef USE_DYNAMIC_ALLOC
	else if (strcmp(argv[i], "-len") == 0) {
	    veclen = atoi(argv[i+1]);
	    LEN    = veclen;
	    i++;
	}
#endif
	else {
	    fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
	    exit(1);
	}
    }

    /* Determine bx, by, and bz so that (bx+2)(by+2)(bz+2) <= LEN.
       Let bx=by=bz for now */
    l3 = pow((double)LEN,1./3.);
    bz = by = bx = l3-2;
    /* Check */
    if ((bx+2)*(by+2)*(bz+2) > LEN) {
	fprintf(stderr, "Error in computing bx,by,bz = (%d,%d,%d)\n",
		bx, by, bz);
	exit(1);
    }

#ifdef USE_DYNAMIC_ALLOC
    posix_memalign((void **)&a, memalignsize, LEN*sizeof(BASETYPE));
    posix_memalign((void **)&b, memalignsize, LEN*sizeof(BASETYPE));
    posix_memalign((void **)&c, memalignsize, LEN*sizeof(BASETYPE));
    posix_memalign((void **)&d, memalignsize, LEN*sizeof(BASETYPE));
    posix_memalign((void **)&e, memalignsize, LEN*sizeof(BASETYPE));
#endif
    fprintf(outfp, "Basetype is %s, data length %d\n", BASETYPE_NAME, LEN);
    fprintf(outfp, "Stencil size (%d,%d,%d)\n", bx, by, bz);
#ifdef COMPILE_OPTS
    fprintf(outfp, "Compile options\t%s\n", COMPILE_OPTS);
#endif
    fprintf(outfp, "Loop\tTime(sec)\tRate    \tChecksum\tDesc\n");

    /* For running on the host instead of the accelerator, do we
       need to call acc_set_num_cores(NUM_CORES);
       where NUM_CORES is defined to the number of cores to use? */

#pragma acc enter data create(a[LEN],b[LEN],c[LEN],d[LEN],e[LEN])
    vpvts(s);
    vsumr();
    paxpy(s);
    paxpby(s, s1);
    paxpbypcz(s, s1, s2);
    pmdot();
    stencil3d(bx, by, bz);
#pragma acc exit data copyout(a[LEN],b[LEN],c[LEN],d[LEN],e[LEN])

    /* Ensure that the copyout happens */
    dummy(a, b, c, d, e, aa, bb, cc, 0);

    if (outfp != stdout) fclose(outfp);

#ifdef USE_DYNAMIC_ALLOC
free(a);
free(b);
free(c);
free(d);
free(e);
#endif

    return 0;
}

void printSummary(const char *name, double t, long ops, BASETYPE checkval,
		  const char *shortdesc)
{
    double rate = 0.0;
    if (t > 0) rate = (double)ops/t;
    fprintf(outfp, "%s\t%.2e\t%.2e\t%.6e\t%s\n", name, t, rate, checkval,
	    shortdesc ? shortdesc : "");
}

/* This routine computes a checksum for a result of length len. It attempts
   to force the compiler to perform the operation in scalar mode. */
void getCheck(BASETYPE *a, int len, BASETYPE *checkval)
{
    BASETYPE sum = 0.0;

#pragma acc update self(a[0:len])

#if defined(_CRAYC)
#pragma _CRI novector
#elif defined(__xlc__) || defined(PGI_COMPILER) || defined(__INTEL_COMPILER)
#pragma novector
#endif
    for (int i=0; i<len; i++) sum += a[i];
    *checkval = sum;
}

/* For scalar inputs, consider storing those in a device array as part of the
   init step, then moving them into a local variable inside the structured
   block in order to avoid the potential cost of a data transfer */
/* Here start the vector test routines */
/* ---------------------------------------------------------------------- */
int vpvts(BASETYPE s)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(a,b)
#pragma acc parallel
    {
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 1.0);
    }

    start_t = clock();
/* Note that this requires the value s to be moved to the accelerator.
   We could instead use a short array for scalar motion and use a scalar
   solely on the accelerator to avoid including the time to move the
   scalar to/from the accelerator */
#pragma acc data present(a,b)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		a[i] += b[i] * s;
	    }
	    dummyab(a, b);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    getCheck(a, LEN, &checkval);
    printSummary("vpvts", clock_dif_sec, 2*LEN, checkval,
		 "vector plus vector times scalar");

    return 0;
}

int vsumr(void)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval, sum;
    SETVECLEN;

#pragma acc data present(a)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 1.0);
    }

    start_t = clock();
#pragma acc data present(a,b)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel reduction (+:sum)
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		sum += a[i];
	    }
	    dummyab(a, b);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    /* Query: This forces the system to move the sum scalar at the end of
       the accelerator section - consider making sum a local variable,
       and assigning it to a[0] */
    checkval = sum;
    printSummary("vsumr", clock_dif_sec, LEN, checkval, "vector sum reduction");

    return 0;
}


/* PETSc - These are drawn from the vector operations that PETSc provides
   in the PETSc Vector */

int paxpy(BASETYPE s)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(a,b)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 1.0);
    }

    start_t = clock();
#pragma acc data present(a,b)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		b[i] = s * a[i] + b[i];
	    }
	    dummyab(a,b);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    getCheck(b, LEN, &checkval);
    printSummary("paxpy", clock_dif_sec, LEN*2, checkval, "Y=aX+Y");

    return 0;
}

int paxpby(BASETYPE s1, BASETYPE s2)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(a,b)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 0.3);
    }

    start_t = clock();
#pragma acc data present(a,b)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		b[i] = s1 * a[i] + s2 * b[i];
	    }
	    dummyab(a, b);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    getCheck(b, LEN, &checkval);
    printSummary("paxpby", clock_dif_sec, LEN*3, checkval, "Y=aX+bY");

    return 0;
}

int paxpbypcz(BASETYPE s1, BASETYPE s2, BASETYPE s3)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(a, b, c)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.3);
	set1d(b, LEN, 1.0);
	set1d(c, LEN, 0.5);
    }

    start_t = clock();
#pragma acc data present(a, b, c)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		c[i] = s1 * a[i] + s2 * b[i] + s3 * c[i];
	    }
	    dummyabc(a, b, c);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    getCheck(c, LEN, &checkval);
    printSummary("paxpbypcz", clock_dif_sec, LEN*5, checkval, "Z=aX+bY+cZ");

    return 0;
}

int pmdot(void)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE s1, s2, s3, s4;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(a,b,c,d,e)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 1.0);
	set1d(c, LEN, 0.3);
	set1d(d, LEN, -0.2);
	set1d(e, LEN, -0.5);
    }

    start_t = clock();
#pragma acc data present(a,b,c,d,e)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
	    s1 = s2 = s3 = s4 = 0;
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		s1 += a[i] * b[i];
		s2 += a[i] * c[i];
		s3 += a[i] * d[i];
		s4 += a[i] * e[i];
	    }
	    dummy(a, b, c, d, e, aa, bb, cc, s1+s2+s3+s4);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    checkval = s1+s2+s3+s4;
    printSummary("pmdot", clock_dif_sec, LEN*4*2, checkval, "Multiple dot product");

    return 0;
}

/* Other PETSc routines include max and min by element, and element-wise
   product and division; these are not yet tested here. */

/* ---------------------------------------------------------------------- */
#define ind(i,j,k) (k)*((by+2)*(bx+2))+(j)*(bx+2)+(i)
int stencil3d(int bx, int by, int bz)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(LIST-OF-VARS-TO-INITIALIZE)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 1.0);
    }

    start_t = clock();
#pragma acc data present(LIST-OF-VARS-INCLUDES-ALL-INCLUDING-CALL-TO-DUMMY)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
/* NOTE: THIS IS UNDER DEVELOPMENT.  IT PROBABLY NEEDS A LOOP COLAPSE */
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int k=1; k<bz+1; ++k) {
		for(int j=1; j<by+1; ++j) {
		    for(int i=1; i<bx+1; ++i) {
			b[ind(i,j,k)] = a[ind(i,j,k)]/2.0 + (a[ind(i-1,j,k)] + a[ind(i+1,j,k)] + a[ind(i,j-1,k)] + a[ind(i,j+1,k)] + a[ind(i,j,k-1)] + a[ind(i,j,k+1)])/6.0/2.0;
		    }
		}
	    }
	    dummyab(a, b);
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    /* Perform a computation on a result of the operation to serve as
       a check on the calculation. For vector output, use getCheck.
       If there are multiple outputs, use getCheck on each one and add
       the results together */
    getCheck(b, LEN, &checkval);
    printSummary("stencil3d", clock_dif_sec, 8*bx*by*bz, checkval,
		 "3d stencil");

    return 0;
}

/* Routine Template */
#if 0
int vpvts(BASETYPE s)
{
    clock_t  start_t, end_t, clock_dif;
    double   clock_dif_sec;
    BASETYPE checkval;
    SETVECLEN;

#pragma acc data present(LIST-OF-VARS-TO-INITIALIZE)
#pragma acc parallel
    {
	/* Initialize all variables used in test on GPU */
	set1d(a, LEN, 0.5);
	set1d(b, LEN, 1.0);
    }

    start_t = clock();
#pragma acc data present(LIST-OF-VARS-INCLUDES-ALL-INCLUDING-CALL-TO-DUMMY)
#ifdef ACC_USE_KERNELS
#pragma acc kernels
#else
#pragma acc parallel
#endif
    {
	for (int nl = 0; nl < ntimes; nl++) {
#ifndef ACC_USE_KERNELS
#pragma acc loop
#endif
	    for (int i = 0; i < LEN; i++) {
		VECTOR-LOOP-TO-COMPUTE;
	    }
	    dummy(a, b, c, d, e, aa, bb, cc, 0.);
	    // or dummyab(a,b) or dummyabc(a,b,c)
	}
    }
    end_t = clock(); clock_dif = end_t - start_t;
    clock_dif_sec = (double) (clock_dif/1000000.0);
    clock_dif_sec = clock_dif_sec/ntimes;
    /* Perform a computation on a result of the operation to serve as
       a check on the calculation. For vector output, use getCheck.
       If there are multiple outputs, use getCheck on each one and add
       the results together */
    getCheck(a, LEN, &checkval);
    printSummary("NAME", clock_dif_sec, NUMBER-OF-OPS, checkval, "DESC");

    return 0;
}
#endif
