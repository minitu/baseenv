/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include "mpi.h"
#include <string.h>
#include "halocompare.h"
#include "persist.h"

static int runPersistTest( HaloElement *halo, int size, HaloTest *self, 
			   double *time )
{
    MPI_Request r[2*MAX_NEIGHBOR];
    int nNbr  = halo->nNeighbors;
    int nTest = halo->nTest;
    double t;
    int i, k;

    for (k=0; k<nNbr; k++) {
	MPI_Send_init( halo->sbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], 0, halo->comm, &r[k] );
	MPI_Recv_init( halo->rbuf[k], size, MPI_DOUBLE, 
		       halo->partners[k], 0, halo->comm, &r[nNbr+k] );
    }

    MPI_Barrier( MPI_COMM_WORLD );
    
    t = MPI_Wtime();
    for (i=0; i<nTest; i++) {
	MPI_Startall( 2*nNbr, r );
	MPI_Waitall( 2*nNbr, r, MPI_STATUSES_IGNORE );
    }
    t = MPI_Wtime() - t;

    for (k=0; k<2*nNbr; k++) {
	MPI_Request_free( &r[k] );
    }
    
    if (t < 100* MPI_Wtick()) {
	return 1;
    }
    *time = t / nTest;
    return 0;
}

static int freePersistTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

static int initPersistTest( HaloElement *halo, int size, HaloTest *self )
{
    return 0;
}

int setupPersistTest( HaloTest *self )
{
    strncpy( self->testName, "Persist", sizeof(self->testName) );
    self->initTest = initPersistTest;
    self->runTest  = runPersistTest;
    self->freeTest = freePersistTest;

    self->data     = 0;
    return 0;
}
