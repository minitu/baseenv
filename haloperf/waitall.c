/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include "halocompare.h"
#include "waitall.h"

/* Forward refs */

static int runWaitallTest( HaloElement *halo, int size, HaloTest *self, 
			   double *time )
{
    MPI_Request r[2*MAX_NEIGHBOR];
    int nNbr  = halo->nNeighbors;
    int nTest = halo->nTest;
    double t;
    int i, k;

    MPI_Barrier( MPI_COMM_WORLD );
    
    t = MPI_Wtime();
    for (i=0; i<nTest; i++) {
	for (k=0; k<nNbr; k++) {
	    MPI_Isend( halo->sbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], i, halo->comm, &r[k] );
	    MPI_Irecv( halo->rbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], i, halo->comm, &r[nNbr+k] );
	}
	MPI_Waitall( 2*nNbr, r, MPI_STATUSES_IGNORE );
    }
    t = MPI_Wtime() - t;

    if (t < 100* MPI_Wtick()) {
	return 1;
    }
    *time = t / nTest;
    return 0;
}

static int freeWaitallTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

static int initWaitallTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

int setupWaitallTest( HaloTest *self )
{
    strncpy( self->testName, "Waitall", sizeof(self->testName) );
    self->initTest = initWaitallTest;
    self->runTest  = runWaitallTest;
    self->freeTest = freeWaitallTest;

    self->data     = 0;
    return 0;
}

/* This is a special test for *either* send or *receive*, not both */
static int runWaitallTest1( HaloElement *halo, int size, HaloTest *self, 
			   double *time )
{
    MPI_Request r[MAX_NEIGHBOR];
    int nNbr  = halo->nNeighbors;
    int nTest = halo->nTest;
    int master = halo->isMaster;
    double t;
    int i, k;

    MPI_Barrier( MPI_COMM_WORLD );
    
    t = MPI_Wtime();
    for (i=0; i<nTest; i++) {
	if (master) {
	    for (k=0; k<nNbr; k++) {
		/* printf( "sending to %d\n", halo->partners[k] ); */
		MPI_Isend( halo->sbuf[k], size, MPI_DOUBLE, 
			   halo->partners[k], i, halo->comm, &r[k] );
	    }
	}
	else {
	    for (k=0; k<nNbr; k++) {
		/* printf( "receiving from %d\n", halo->partners[k] ); */
		MPI_Irecv( halo->rbuf[k], size, MPI_DOUBLE, 
			   halo->partners[k], i, halo->comm, &r[k] );
	    }
	}
	/* fflush(stdout); */
	MPI_Waitall( nNbr, r, MPI_STATUSES_IGNORE );
    }
    t = MPI_Wtime() - t;

    if (t < 100* MPI_Wtick() && nNbr > 0) {
	return 1;
    }
    *time = t / nTest;
    return 0;
}

int updateWaitall1Test( HaloTest *self )
{
    strncpy( self->testName, "Waitall (1 sender)", sizeof(self->testName) );
    self->runTest  = runWaitallTest1;
    
    return 0;
}
