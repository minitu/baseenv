/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include "mpi.h"
#include <string.h>
#include "halocompare.h"
#include "irecv.h"

static int runIrecvTest( HaloElement *halo, int size, HaloTest *self, 
			 double *time )
{
    MPI_Request r[MAX_NEIGHBOR];
    int nNbr  = halo->nNeighbors;
    int nTest = halo->nTest;
    double t;
    int i, k;

    MPI_Barrier( MPI_COMM_WORLD );
    
    t = MPI_Wtime();
    for (i=0; i<nTest; i++) {
	for (k=0; k<nNbr; k++) {
	    MPI_Irecv( halo->rbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], i, halo->comm, &r[k] );
	}
	for (k=0; k<nNbr; k++) {
	    MPI_Send( halo->sbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], i, halo->comm );
	}
	MPI_Waitall( nNbr, r, MPI_STATUSES_IGNORE );
    }
    t = MPI_Wtime() - t;

    if (t < 100* MPI_Wtick()) {
	return 1;
    }
    *time = t / nTest;
    return 0;
}

static int freeIrecvTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

static int initIrecvTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

int setupIrecvTest( HaloTest *self )
{
    strncpy( self->testName, "Send/Irecv", sizeof(self->testName) );
    self->initTest = initIrecvTest;
    self->runTest  = runIrecvTest;
    self->freeTest = freeIrecvTest;

    self->data     = 0;
    return 0;
}
