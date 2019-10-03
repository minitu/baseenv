/* -*- Mode: C; c-basic-offset:4 ; -*- */

/* haloperf.h has the results from configure */
#include "haloperf.h"

#define MAX_NEIGHBOR 8
#define MAX_DIMS     3

/* Temporary to hold topology information - should be provided by the
   Topology API */
typedef struct TopoDesc { int coords[3]; } TopoDesc;

/* Describe the "halo" */
typedef struct {
    int nNeighbors;
    int nTest, isMaster;
    double *sbuf[MAX_NEIGHBOR], *rbuf[MAX_NEIGHBOR]; 
    int partners[MAX_NEIGHBOR];
    int partnerClique[MAX_NEIGHBOR];   /* Used to support options for
					  SMP process topology */
    TopoDesc partnerCoords[MAX_NEIGHBOR];
    int dims[MAX_DIMS], coords[MAX_DIMS], ndims;
    MPI_Comm comm;
} HaloElement;

/* Describe a test */
typedef struct HaloTest {
    char testName[40];
    int  (*initTest)( HaloElement *, int, struct HaloTest * );
    int  (*runTest)( HaloElement *, int, struct HaloTest *, double * );
    int  (*freeTest)( HaloElement *, int, struct HaloTest * );
    void *data;
} HaloTest;  

/* Description of the tests and their results.  This combines a particular 
   halo exchange implementation (HaloTest) with a description of the 
   halo data (HaloElement) */
typedef struct TestResult {
    HaloTest    *test;
    HaloElement *halo;
    int          size;
    double       maxTime, minTime, avgTime;
} TestResult;
