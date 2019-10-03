/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "halocompare.h"
#include "phased.h"

#define MAX_COMM_PHASES 4

/* This structure contains information on the messages to send and receive */
typedef struct {
    int nSendPartners, SendPartners[MAX_NEIGHBOR];
    int nRecvPartners, RecvPartners[MAX_NEIGHBOR];
    void *sBuf[MAX_NEIGHBOR], *rBuf[MAX_NEIGHBOR];
} ExchangePhase;

/* An exchange is made up of one or more communication phases */
typedef struct {
    int           nPhases;         /* Number of phases of communication */
    ExchangePhase phase[MAX_COMM_PHASES];
} Exchange;

static void printPhases( HaloTest *self );


static int runPhasedTest( HaloElement *halo, int size, HaloTest *self, 
			  double *time )
{
    MPI_Request r[2*MAX_NEIGHBOR];
    int nTest = halo->nTest;
    int phase;
    MPI_Comm comm = halo->comm;
    double t;
    int i, j;
    Exchange *desc = (Exchange *)self->data;

    MPI_Barrier( MPI_COMM_WORLD );
    
    t = MPI_Wtime();
    for (j=0; j<nTest; j++) {
	for (phase=0; phase<desc->nPhases; phase++) {
	    ExchangePhase * const ephase = &(desc->phase[phase]);
	    int rcnt = 0;
	    /* rationale: Get flits out the door, then start the
	       receives so that they'll match as the isend envelopes 
	       arrive.  An alternative is to post the receives first
	       to ensure that the sends match as they arrive */
	    for (i=0; i<ephase->nSendPartners; i++) {
		MPI_Isend( ephase->sBuf[i], size, MPI_DOUBLE,
			   ephase->SendPartners[i], j, comm, &r[rcnt++] );
	    }
	    for (i=0; i<ephase->nRecvPartners; i++) {
		MPI_Irecv( ephase->rBuf[i], size, MPI_DOUBLE,
			   ephase->RecvPartners[i], j, comm, &r[rcnt++] );
	    }
	    MPI_Waitall( rcnt, r, MPI_STATUSES_IGNORE );
	}
    }
    t = MPI_Wtime() - t;

    
    if (t < 100* MPI_Wtick()) {
	return 1;
    }
    *time = t / nTest;
    return 0;
}

static int freePhasedTest( HaloElement *halo, int size, HaloTest *self )
{
    free( self->data );
    return 0;
}

static int initPhasedTest( HaloElement *halo, int size, HaloTest *self )
{
    Exchange *desc = (Exchange *)malloc( sizeof(Exchange) );
    int i, isOdd, sendPhase, recvPhase;
	
    self->data = (void *)desc;
    
    /* Compute the phases based on the stencil type, inferred from the
       number of partners (2, 4 = star, 8 = box).
       The rule is:
       star:  2 phases: odd, then even senders (based on sum of x,y
              coordinates).  This is easy, because each process
	      does all of its sending or receiving in a single phase,
	      thus we need only compute which phase a point belongs to
       box:   4 phases:  ( x % 2 ) + 2 ( y % 2),
              where (x,y) are the coordinates of the process and
              % is the mod operation

      One alternative is to have 2 groups of 2 phases:
      the star stencil (2 phases, based on the same odd/even
      computation) and the complement (a "times" stecil), 
      also broken into 2 phases).  
      The other alternative is to 

      star: 2 or 4 neighbors
      box : 8 neighbors

      The Box version is not implemented yet.  A more general approach
      can be used by picking a color for the rank=0 process, and have that
      process communicate to its partners, who then pick another color.
      The communication continues (in parallel) until a coloring is determined
      (the devil is in the details, of course).  
	*/
    desc->nPhases = 2;
    isOdd = 0;
    for (i=0; i<halo->ndims; i++) 
	isOdd += halo->coords[i];
    isOdd = isOdd & 0x1;
    sendPhase = isOdd;
    recvPhase = 1 - isOdd;
    desc->phase[sendPhase].nSendPartners = halo->nNeighbors;
    desc->phase[sendPhase].nRecvPartners = 0;
    desc->phase[recvPhase].nRecvPartners = halo->nNeighbors;
    desc->phase[recvPhase].nSendPartners = 0;
    for (i=0; i<halo->nNeighbors; i++) {
	desc->phase[sendPhase].SendPartners[i] = halo->partners[i];
	desc->phase[sendPhase].sBuf[i]         = halo->sbuf[i];
	desc->phase[recvPhase].RecvPartners[i] = halo->partners[i];
	desc->phase[recvPhase].rBuf[i]         = halo->rbuf[i];
    }

    if (0)
	printPhases( self );
    
    return 0;
}

int setupPhasedTest( HaloTest *self )
{
    strncpy( self->testName, "Phased", sizeof(self->testName) );
    self->initTest = initPhasedTest;
    self->runTest  = runPhasedTest;
    self->freeTest = freePhasedTest;

    self->data     = 0;
    return 0;
}

static void printPhases( HaloTest *self )
{
    Exchange *desc = (Exchange *)self->data;
    int r, size, rank;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    for (r=0; r<size; r++) {
	int phase, i;
	MPI_Barrier( MPI_COMM_WORLD );
	if (r == rank) {
	    for (phase=0; phase < desc->nPhases; phase ++ ) {
		if (desc->phase[phase].nSendPartners) {
		    printf( "%d: %d senders (%d)\n%d:", 
			    rank, phase, desc->phase[phase].nSendPartners, 
			    rank );
		    for (i=0; i<desc->phase[phase].nSendPartners; i++) {
			printf( "%d ", desc->phase[phase].SendPartners[i] );
		    }
		    printf( "\n" );
		}
		if (desc->phase[phase].nRecvPartners) {
		    printf( "%d: %d receivers (%d)\n%d:", 
			    rank, phase, desc->phase[phase].nRecvPartners,
			    rank );
		    for (i=0; i<desc->phase[phase].nRecvPartners; i++) {
			printf( "%d ", desc->phase[phase].RecvPartners[i] );
		    }
		    printf( "\n" );
		}
	    }
	}
	fflush( stdout );
    }
}

