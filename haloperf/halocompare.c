/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 
   This program probes the behavior of an MPI implementation on 
   simple halo exchanges.  It is designed to be easily (though not trivially) 
   extended, and to test several "dimensions" of options.

   Process Layout
   There are two broad categories of process layout:
   Using MPI_COMM_WORLD and an application-computed layout
   Using MPI_Cart_create to create a "good" layout.
 */
#include "halocompare.h"
/*#include "meshtopo.h" */
#include "topoinfo.h"
/* Add the prototypes for each methods setup function */
#include "irecv.h"
#include "waitall.h"
#include "persist.h"
#include "phased.h"
/* Add prototypes for topology discovery */
#if 0
#ifdef HAVE_FINDCLIQUES_H 
#include "findcliques.h"
#endif
#endif

/* Update this value when more communicators are defined in initHaloComm */
#define MAX_HALOCOMM_TYPE 5

/* Verbose can be 0 (no output), 1 (output from rank==0), or >1 (output from 
   all ranks) */
static int verbose = 0;
static int wrank = -1;       /* Rank of process in MPI_COMM_WORLD */
char *title = 0;
#define MAX_RESULTS 128
int nResults = 0;
TestResult results[MAX_RESULTS];
int hasTopo = 0;
int hasTopoMesh = 0;
topoinfo_t *tinfo;

/* Routines to help with tests */
int initHaloComm( HaloElement *, int, int );
int initHaloElement( HaloElement *, int, int );
int setNeighbors2( HaloElement *, int );
int setNeighbors2y( HaloElement *, int );
int setNeighbors4( HaloElement *, int );
int setNeighbors8( HaloElement *, int );
int MakeOnePoint( HaloElement *, HaloElement * );
int coordToRank2( int, int, int [], int );
int FindTimes( double, double *, double *, double * );
int ConfirmPattern( HaloElement * );
void IntFindMinMax( int v, int *vmin, int *vmax, MPI_Comm comm, int commrank );
void printMeshNodes( FILE * );

/* Make it easier to track down problems */
#define TRACE
#ifdef TRACE
#define FUNCENTER(a) {if (verbose > 1 || (verbose == 1 && wrank == 0)){printf("Entering " a "\n" ); fflush(stdout);}}
#define FUNCEXIT(a)  {if (verbose > 1 || (verbose == 1 && wrank == 0)){printf("Exiting " a "\n" ); fflush(stdout);}}
#else
#define FUNCENTER(a)
#define FUNCEXIT(a)
#endif

#define WAITALLTEST 0x1
#define IRECVTEST   0x2
#define PERSISTTEST 0x4
#define PHASEDTEST  0x8
/* Turn on all tests */
static int enabledTests = -1;
/* Save information about possible node hierarchy */
#define MAX_NODE_SIZE 256
static int nodeSize = -1;
static int numnodes = -1;
static int myCliqueNum  = -1;
static int myCliqueSize = -1;
static int myCliqueRanks[MAX_NODE_SIZE];
static int minClique = -1, maxClique = 1;

int main( int argc, char *argv[] )
{
    int size, err, type;
    int desiredNbrs = 2;
    int isPeriodic  = 0;
    int i, k, nTestComm = 0, maxCount = 2048;
    HaloElement  halo[16], halo1world, halo1cart, halo1smp;
    HaloTest     haloTests[16], halotest1world, halotest1cart, halotest1smp;
    char myName[MPI_MAX_PROCESSOR_NAME+1];

    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );

    /* See if topology information is available.  Depending on how
       the topology package is configured, this will determine which
       processes are on the same node (clique) and the properties of the
       network. */
    if (topoInit(0, &tinfo) ==0) hasTopo = 1;

    /*if (TopoInitMesh() == 0) hasTopo = 1; */

    /* Get modifiers on the tests (maximum partners, message sizes) */
    for (i=1; i<argc; i++) {
	if (strcmp(argv[i],"-nbrs") == 0) {
	    i++;
	    desiredNbrs = strtol( argv[i], 0, 10 );
	}
	else if (strcmp( argv[i], "-periodic" ) == 0) {
	    isPeriodic = 1;
	}
	else if (strcmp( argv[i], "-count" ) == 0) {
	    i++;
	    maxCount = strtol( argv[i], 0, 10 );
	}
	else if (strcmp( argv[i], "-title" ) == 0) {
	    i++;
	    title = strdup( argv[i] );
	}
	else if (strcmp( argv[i], "-verbose" ) == 0 ||
		 strcmp( argv[i], "-v" ) == 0) {
	    verbose++;
	}
	else if (strcmp( argv[i], "-nowaitall" ) == 0) {
	    enabledTests = enabledTests & ~WAITALLTEST;
	}
	else if (strcmp( argv[i], "-noirecvtest" ) == 0) {
	    enabledTests = enabledTests & ~IRECVTEST;
	}
	else if (strcmp( argv[i], "-nopersist" ) == 0) {
	    enabledTests = enabledTests & ~PERSISTTEST;
	}
	else if (strcmp( argv[i], "-nophased" ) == 0) {
	    enabledTests = enabledTests & ~PHASEDTEST;
	}
	else if (strcmp(argv[i],"-smpsize") == 0) {
	    i++;
	    nodeSize = strtol( argv[i], 0, 10 );
	}
	else if (strcmp( argv[i], "-smp" ) == 0) {
	    /* Try to extract this information the topo data */
	    myCliqueSize = MAX_NODE_SIZE;
	    err = topoNodeEnumeration(tinfo, &numnodes,
				      &myCliqueNum,
				      &myCliqueSize,
				      myCliqueRanks);
/*
	    err = MPE_FindCliqueFromNodename( MPI_COMM_WORLD, 
					      16, MAX_NODE_SIZE, 
					      &myCliqueNum, &myCliqueSize,
					      myCliqueRanks );
*/
	    if (err) {
		fprintf( stderr, "Unable to find cliques when -smp selected; error return %d\n", err );
		MPI_Abort( MPI_COMM_WORLD, 1 );
	    }
	}
	else {
	    if (wrank == 0) {
		fprintf( stderr, "Unrecognized arg %s\n", argv[i] );
		fflush( stderr );
	    }
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
    }

    /* All processes will either have <= 0 or > 0 for this value */
    if (myCliqueSize > 0) {
	IntFindMinMax( myCliqueSize, &minClique, &maxClique,
		       MPI_COMM_WORLD, 0 );
    }

    if (wrank == 0) {
	int  nameLen;
	char *s = getenv( "HALO_TITLE" );
	MPI_Get_processor_name( myName, &nameLen );
	printf( "Test with %d processes on %s %s\n", size, myName,
		isPeriodic ? "(periodic)" : "" );
	if (nodeSize > 0)
	    printf( "Using a node size of %d\n", nodeSize );
	if (myCliqueSize > 0) {
	    printf( "Process 0 is in a clique of %d processes\n", myCliqueSize );
	    printf( "Cliques range from %d to %d in size\n", minClique, maxClique );
	}
	if (hasTopo) {
	    int mindims[5], maxdims[5], qtorus[5], ndims=5;
	    err = topoMeshContainer(tinfo, &ndims, mindims, maxdims, qtorus);
            if (!err) {
                hasTopoMesh = 1;
	        printf("Processes within ");
	        for (i=0; i<ndims; i++) {
	            printf("[%d,%d]%c", mindims[i], maxdims[i],
		           (i<ndims-1)?'x':'\n');
                }
	    }
	}
	/* Add any extra title information to the output page. */
	if (title && title[0]) {
	    printf( "%s\n", title );
	    /* Note allocated above, free here */
	    free( title );
	}
	if (s && s[0]) {
	    printf( "%s\n", s );
	}
	
    }
    /* Make sure all processes agree on whether the topology is a mesh */
    MPI_Bcast(&hasTopoMesh, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Initialize the tests */
    if (enabledTests & WAITALLTEST)
	setupWaitallTest( &haloTests[nTestComm++] );
    if (enabledTests & IRECVTEST)
	setupIrecvTest( &haloTests[nTestComm++] );
    if (enabledTests & PERSISTTEST)
	setupPersistTest( &haloTests[nTestComm++] );
    /* The phased test only works for star stencils */
    if ((enabledTests & PHASEDTEST) && desiredNbrs < 8)
	setupPhasedTest( &haloTests[nTestComm++] );

    /* Initialize the experiements */
    /* This particular ordering of the loops will give us, for each communicator
       type, the behavior of each of the implemented methods */
    nResults = 0;
    /* type is for the type of Communicator */
    for (type=0; type<=MAX_HALOCOMM_TYPE; type++) {
	initHaloElement( &halo[type], MAX_NEIGHBOR, maxCount );
	if (halo[type].comm != MPI_COMM_WORLD && 
	    halo[type].comm != MPI_COMM_NULL) { 
	    MPI_Comm_free( &halo[type].comm ); 
	}
	initHaloComm( &halo[type], type, isPeriodic );
	/* Skip unavailable communicators */
	if (halo[type].comm == MPI_COMM_NULL) {
	    continue;
	}
	switch (desiredNbrs) {
	case 2:
	    setNeighbors2( &halo[type], isPeriodic );
	    break;
	case 4:
	    setNeighbors4( &halo[type], isPeriodic );
	    break;
	case 8:
	    setNeighbors8( &halo[type], isPeriodic );
	    break;
	default:
	    fprintf( stderr, "%d neighbors not supported\n", desiredNbrs );
	    fflush(stdout);
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	if (ConfirmPattern( &halo[type] )) {
	    fprintf( stderr, "PANIC: Bad pattern for type %d\n", type );
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	/* Get basic information about the node heirarchy.  This can be
	   used to determine how many "off-node" communication links 
	   are used in different process to processor mappings. */
	if (myCliqueSize > 0) {
	    TopoFindPartnerCliquenum( halo[type].nNeighbors, 
				      halo[type].partners, halo[type].comm, 
				      myCliqueNum, halo[type].partnerClique );
	}

	if (verbose) {
	    char commName[MPI_MAX_OBJECT_NAME+1];
	    int  nameLen;
	    if (wrank == 0) {
		printf( "Mesh dims = [%d,%d]\n", halo[type].dims[0],
			halo[type].dims[1] );
		fflush(stdout);
	    }
	    MPI_Comm_get_name( halo[type].comm, commName, &nameLen );
	    if (verbose > 1) {
		/* FIXME: Needs a synchronized print */
		printf( "Partners for %s (%d,%d):\t", commName, 
			halo[type].coords[0], halo[type].coords[1] );
		for (i=0; i<halo[type].nNeighbors; i++) {
		    printf( "%d ", halo[type].partners[i] );
		}
		printf( "\n" );
	    }
	    else {
		printf( "Communicator %s (use -v -v for detailed partners)\n", 
			commName );
	    }
	    fflush(stdout);
	}

	if (nResults >= MAX_RESULTS) {
	    fprintf( stderr, "Too many tests created (%d)\n", nResults );
	    fflush( stderr );
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}

	for (k=0; k<nTestComm; k++) {
	    results[nResults].test = &haloTests[k];
	    results[nResults].halo = &halo[type];
	    results[nResults].size = maxCount;
	    results[nResults].maxTime = 0;
	    nResults++;
	}
    }

    /* Add a test for the single sender */
    /* 
       Still to do:
       If SMP aware, for a 2-d mesh, need these tests:
       1) All neighbors in clique
       2) All but one neighbor in clique
       3) Two without (vertex in 2-d mesh)
       4) All processes in the same clique send and receive as normal;
          others receive as necessary.
     */
    if (MakeOnePoint( &halo[0], &halo1world )) {
	fprintf( stderr, "Aborting because MakeOnePoint failed\n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    setupWaitallTest( &halotest1world );
    updateWaitall1Test( &halotest1world );
    results[nResults].test    = &halotest1world;
    results[nResults].halo    = &halo1world;
    results[nResults].size    = maxCount;
    results[nResults].maxTime = 0;
    nResults++;

    /* Use the cart communicator */
    if (MakeOnePoint( &halo[3], &halo1cart )) {
	fprintf( stderr, "Aborting because MakeOnePoint failed\n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    setupWaitallTest( &halotest1cart );
    updateWaitall1Test( &halotest1cart );
    results[nResults].test    = &halotest1cart;
    results[nResults].halo    = &halo1cart;
    results[nResults].size    = maxCount;
    results[nResults].maxTime = 0;
    nResults++;

    /* Use the smp-aware communicator, if available */
    if (halo[5].comm != MPI_COMM_NULL) {
	if (MakeOnePoint( &halo[5], &halo1smp )) {
	    fprintf( stderr, "Aborting because MakeOnePoint failed\n" );
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	setupWaitallTest( &halotest1smp );
	updateWaitall1Test( &halotest1smp );
	results[nResults].test    = &halotest1smp;
	results[nResults].halo    = &halo1smp;
	results[nResults].size    = maxCount;
	results[nResults].maxTime = 0;
	nResults++;
    }


    /* For each test, run the timings */
    for (i=0; i<nResults; i++) {
	/* For each timing, run it several times to find a stable minimum 
	   time */
	double minTime, testTime, testMinTime, testAvgTime, t, tstats[3];
	int    nTrial=10, nSame=0;

	if (verbose && wrank == 0) {
	    printf( "Starting test %d (%s)\n", i, results[i].test->testName );
	    fflush(stdout);
	}
	minTime = 1.0e6;
	results[i].test->initTest( results[i].halo, results[i].size, 
				   results[i].test );
	while (nTrial-- && nSame < 3) {
	    err = results[i].test->runTest( results[i].halo, results[i].size, 
					    results[i].test, &testTime );
	    if (err) {
		if (wrank == 0) {
		    printf( "Error running test %d\n", i );
		}
		break;
	    }
	    t = testTime;
	    FindTimes( t, &testTime, &testMinTime, &testAvgTime );
	    /* We select on the maximum time measured, but store the
	       other statistics (max,min,average) */
	    if (testTime < minTime) {
		minTime = testTime;
		tstats[0] = testTime;
		tstats[1] = testMinTime;
		tstats[2] = testAvgTime;
		nSame = 0;
	    }
	    else nSame++;
	}
	results[i].test->freeTest( results[i].halo, results[i].size, 
				   results[i].test );
	results[i].maxTime = tstats[0];
	results[i].minTime = tstats[1];
	results[i].avgTime = tstats[2];
    }

    /* Report the results.  Eventually, we could sort these in various ways */
    if (wrank == 0) {
	/* Print the heading - FIXME: not lined up yet */
	printf( "method             \tcomm                   \tnbrs\t#doubles\tmaxT\tminT\trate\n" );
	/* To do: add (min)max # of off node partners if topology support 
	   available */
    }
    for (i=0; i<nResults ; i++) {
	char commName[MPI_MAX_OBJECT_NAME+1];
	int  nameLen, maxNbrs, d;
	int  myOffnode, minOffnode, maxOffnode, totalOffnode;
	int  match, commsize;
	double rate, avgOffnode;
	int  mindist, maxdist, avgdist;
	double dnbrs;
	d = results[i].halo->nNeighbors;
	MPI_Reduce( &d, &maxNbrs, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
	/* If cliqueNum information available, compute statistics
	   on the number of on and off-node (clique) communication */
	if (myCliqueSize > 0) {
	    int j;
	    match = 0;
	    for (j=0; j<results[i].halo->nNeighbors; j++) 
		if (results[i].halo->partnerClique[j] == myCliqueNum) match++;
	    myOffnode = results[i].halo->nNeighbors - match;
	    IntFindMinMax( myOffnode, &minOffnode, &maxOffnode, 
			   MPI_COMM_WORLD, 0 );
	    MPI_Reduce( &myOffnode, &totalOffnode, 1, MPI_INT, MPI_SUM, 
			0, MPI_COMM_WORLD );
	    MPI_Comm_size( results[i].halo->comm, &commsize );
	    avgOffnode = (double)totalOffnode / (double)commsize ;
	}
	if (hasTopoMesh) {
	    int dist[MAX_NEIGHBOR], ii, nbrs;
	    topodist_t *tdinfo;

	    topodistInit(results[i].halo->comm,
			 results[i].halo->nNeighbors,
			 results[i].halo->partners,
			 results[i].halo->nNeighbors,
			 results[i].halo->partners,
			 tinfo, &tdinfo);
	    for (ii=0; ii<results[i].halo->nNeighbors; ii++)
		topoMeshHopDistance(tdinfo, results[i].halo->partners[ii],
				    &dist[ii]);
#if 0
	    TopoMeshDistances( results[i].halo->nNeighbors,
			       results[i].halo->partners,
			       results[i].halo->comm, dist );
#endif
	    mindist = 100000000;
	    maxdist = -1;
	    avgdist = 0;
	    for (ii=0; ii<results[i].halo->nNeighbors; ii++) {
		if (dist[ii] < mindist) mindist = dist[ii];
		if (dist[ii] > maxdist) maxdist = dist[ii];
		avgdist += dist[ii];
	    }
	    if (results[i].halo->nNeighbors == 0) {
		mindist = 0; maxdist = 0;
	    }
	    nbrs = 1.0;
	    if (results[i].halo->nNeighbors > 0) 
		nbrs = results[i].halo->nNeighbors;
	    ii = maxdist;
	    MPI_Reduce( &ii, &maxdist, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
	    ii = mindist;
	    MPI_Reduce( &ii, &mindist, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD );
	    ii = avgdist;
	    MPI_Reduce( &ii, &avgdist, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
	    ii = nbrs;
	    MPI_Reduce( &ii, &nbrs, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
	    dnbrs = nbrs;

	    topodistFree(tdinfo);
		
	}
	if (wrank == 0) {
	    MPI_Comm_get_name( results[i].halo->comm, commName, &nameLen );
	    /* Get overall rate in B/sec .  Use maximum time (minimum rate).
	       This is very important for the 1-sender case, where the min
	       time will be near zero for the processes that have no 
	       communication */
	    rate = results[i].size * sizeof(double) * maxNbrs / 
		results[i].maxTime;
	    /* Covert to MB (not MiB) per sec */
	    rate = rate / 1.e6;
	    printf( "%-20s\t%-23s\t%2d\t%7d\t%.2e\t%.2e\t%4.0f", 
		    results[i].test->testName, commName, maxNbrs, 
		    results[i].size, results[i].maxTime, results[i].minTime, 
		    rate );
	    /* Provide basic information about mapping of processes to 
	       nodes */
	    if (myCliqueSize > 0) 
		printf( "\t(%1d)%1d\t%5.2f", minOffnode, maxOffnode,
			avgOffnode);
	    /* Access mesh topology information if available */
	    if (hasTopoMesh) {
		printf( "\t(%d)%d %5.2f", mindist, maxdist, (double)avgdist/dnbrs );
	    }
	    printf( "\n" );
	    fflush(stdout);
	}
    }

    /* Temporary - to see what nodes were assigned.  This prints once for
       each SMP */
    if (myCliqueSize > 0) {
        printMeshNodes( stdout );
    }

    /* Clean up */
    for (type=0; type<=MAX_HALOCOMM_TYPE; type++) {
	if (halo[type].comm != MPI_COMM_WORLD && 
	    halo[type].comm != MPI_COMM_NULL) { 
	    MPI_Comm_free( &halo[type].comm ); 
	}
    }

    topoFinalize(&tinfo);
    MPI_Finalize();
    
    return 0;
}

/* 
 * Service routines
 */
int initHaloElement( HaloElement *halo, int maxNeighbors, int maxCount )
{
    int i, size;

    FUNCENTER("initHaloElement");

    halo->nNeighbors = maxNeighbors;
    halo->nTest      = 100;
    halo->isMaster   = 0;
    for (i=0; i<maxNeighbors; i++) {
	halo->sbuf[i] = (double *)malloc( maxCount * sizeof(double) );
	halo->rbuf[i] = (double *)malloc( maxCount * sizeof(double) );
	if (!halo->sbuf[i] || !halo->rbuf[i]) {
	    printf( "Unable to allocate halo buffers of size %ld\n",
		    (long)(maxCount * sizeof(double)) );
	    fflush(stdout);
	    return 1;
	}
    }

    /* Defaults for the partners, communicator */
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    halo->comm = MPI_COMM_WORLD;
    for (i=0; i<maxNeighbors; i++) {
	halo->partners[i] = (1+i) % size;
    }

    FUNCEXIT("initHaloElement");
    return 0;
}

/* Based on type, create one of several communicators for a halo test.
   The comm will be named.
   
   Note: The use of this routine references types only up to 5
*/
int initHaloComm( HaloElement *halo, int type, int isPeriodic )
{
    int dims[MAX_DIMS], periods[MAX_DIMS];
    int key, rank, size;
    int i;
#ifdef HAVE_FINDCLIQUES_H
    int err;
#endif

    FUNCENTER("initHaloComm");

    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    for (i=0; i<MAX_DIMS; i++) dims[i] = 0;
    MPI_Dims_create( size, 2, dims );
    for (i=0; i<MAX_DIMS; i++) halo->dims[i] = dims[i];
    for (i=0; i<MAX_DIMS; i++) periods[i] = isPeriodic;
    halo->ndims   = 2;
    /* Mark no coordinates available (yet) */
    for (i=0; i<MAX_DIMS; i++) halo->coords[i] = -1;

    switch (type) {
    case 0: halo->comm = MPI_COMM_WORLD;
	/* Name is already set to comm world */
	break;
    case 1:
	FUNCENTER("MPI_Comm_dup");
	MPI_Comm_dup( MPI_COMM_WORLD, &halo->comm );
	FUNCEXIT("MPI_Comm_dup");
	MPI_Comm_set_name( halo->comm, "dup of world" );
	break;
    case 2:  
	FUNCENTER("MPI even/odd comm");
	if (!(rank & 0x1)) {
	    key = rank/2;
	}
	else {
	    key = (rank - 1) / 2 + (size + 1)/2;
	}
	MPI_Comm_split( MPI_COMM_WORLD, 0, key, &halo->comm );
	FUNCEXIT("MPI even/odd comm");
	MPI_Comm_set_name( halo->comm, "even then odd order" );
	break;
    case 3: 
	FUNCENTER("MPI Cartcomm");
	MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, 1, &halo->comm );
	FUNCEXIT("MPI Cartcomm");
	MPI_Comm_set_name( halo->comm, "2-d Cart" );
	/* Add a check to see if this is any different from COMM_WORLD.
	   To start with, we'll just generate a warning message.
	   Eventually, we can give the program the option of 
	   skipping the tests that use Cart_create */
	{
	    int crank, isSame, allSame;
	    MPI_Comm_rank( halo->comm, &crank );
	    isSame = crank == rank;
	    MPI_Allreduce( &isSame, &allSame, 1, MPI_INT, MPI_LAND, 
			   MPI_COMM_WORLD );
	    if (rank == 0 && allSame) {
		printf( "Cart create is no different than COMM_WORLD\n" );
	    }
	}
	break;
    case 4:
	FUNCENTER("MPI interleave comm");
	if (rank < size/2) {
	    key = 2 * rank;
	}
	else {
	    key = 2*(size - (rank+1)) + 1;
	}
	MPI_Comm_split( MPI_COMM_WORLD, 0, key, &halo->comm );
	FUNCEXIT("MPI interleave comm");
	MPI_Comm_set_name( halo->comm, "Interleave/reverse" );
	break;
    case 5:
#ifdef HAVE_FINDCLIQUES_H
	FUNCENTER("MPI SMPcomm");
	err = MPE_Create2dMeshComm( MPI_COMM_WORLD, 0, 
				    &halo->comm, halo->dims, halo->coords );
	FUNCEXIT("MPI SMPcomm");
	if (err) 
	    halo->comm = MPI_COMM_NULL;
	if (halo->comm != MPI_COMM_NULL) 
	    MPI_Comm_set_name( halo->comm, "SMP-aware 2d mesh" );
	/* FIXME: if verbose, provide some information on the mapping of this
	   communicator. */
#else
	halo->comm = MPI_COMM_NULL;
#endif	
	break;
	/* The remaining types are not used yet.  See the definition 
	   of MAX_HALOCOMM_TYPE */
    case 6: 
	FUNCENTER("MPI 3d cartcomm");
	MPI_Cart_create( MPI_COMM_WORLD, 3, dims, periods, 1, &halo->comm );
	FUNCEXIT("MPI 3d cartcomm");
	MPI_Comm_set_name( halo->comm, "3-d Cart" );
	halo->ndims   = 3;
	break;
#if 0
    case 7:
	/* Skip this case until it is ready */
	if (nodeSize > 0 || myCliqueSize > 0) {
	}
	halo->comm = MPI_COMM_NULL;
	break;
#endif
    default:
	printf( "Unrecognized comm type %d\n", type );
	fflush(stdout);
	FUNCEXIT("initHaloComm");
	return 1;
    }

    FUNCEXIT("initHaloComm");
    return 0;
}

/* Set the partners, given a communicator */
int setNeighbors2( HaloElement *halo, int isPeriodic )
{
    int topoType, rank, drank;
    int i,j;

    FUNCENTER("setNeighbors2");

    MPI_Comm_rank( halo->comm, &rank );
    MPI_Topo_test( halo->comm, &topoType );

    if (topoType == MPI_CART) {
	MPI_Cart_coords( halo->comm, rank, 2, halo->coords );
	/* By using the 1st, rather than 0th, dimension, we get the
	   same partners as the code that doesn't use a cart topology
	   in the case where the comm is just a dup of MPI_COMM_WORLD.
	   This makes comparisions easier to make with the other choices */
	MPI_Cart_shift( halo->comm, 0, 1, 
			&halo->partners[0], &halo->partners[1] );
	if (isPeriodic) 
	    halo->nNeighbors = 2;
	else {
	    /* Compress the number of neighbors */
	    int iIn, iOut = 0;
	    for (iIn = 0; iIn < 2; iIn++) {
		if (halo->partners[iIn] != MPI_PROC_NULL) {
		    halo->partners[iOut++] = halo->partners[iIn];
		}
	    }
	    halo->nNeighbors = iOut;
	}
    }
    else {
	int nbrs=0;
	/* Compute the coords if necessary in the canonical 2-d mesh, 
	   which is row-major */
	if (halo->coords[0] < 0 || halo->coords[1] < 0) {
	    j = rank % halo->dims[1];
	    i = rank / halo->dims[1];
	    halo->coords[0] = i;
	    halo->coords[1] = j;
	}
	else {
	    i = halo->coords[0];
	    j = halo->coords[1];
	    if (verbose) {
		printf( "coords = (%d,%d) in [%d,%d]\n", i, j, 
			halo->dims[0], halo->dims[1] );
	    }
	}

	drank = coordToRank2( i-1, j, halo->dims, isPeriodic );
	/*printf( "%d: rank to left is %d\n", rank, drank );*/
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i+1, j, halo->dims, isPeriodic );
	/*printf( "%d: rank to right is %d\n", rank, drank );*/
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	halo->nNeighbors = nbrs;
    }

    FUNCEXIT("setNeighbors2");
    return 0;
}

int setNeighbors2y( HaloElement *halo, int isPeriodic )
{
    int topoType, rank, drank;
    int i,j;

    FUNCENTER("setNeighbors2");

    MPI_Comm_rank( halo->comm, &rank );
    MPI_Topo_test( halo->comm, &topoType );

    if (topoType == MPI_CART) {
	MPI_Cart_coords( halo->comm, rank, 2, halo->coords );
	/* By using the 1st, rather than 0th, dimension, we get the
	   same partners as the code that doesn't use a cart topology
	   in the case where the comm is just a dup of MPI_COMM_WORLD.
	   This makes comparisions easier to make with the other choices */
	MPI_Cart_shift( halo->comm, 1, 1, 
			&halo->partners[0], &halo->partners[1] );
	if (isPeriodic) 
	    halo->nNeighbors = 2;
	else {
	    /* Compress the number of neighbors */
	    int iIn, iOut = 0;
	    for (iIn = 0; iIn < 2; iIn++) {
		if (halo->partners[iIn] != MPI_PROC_NULL) {
		    halo->partners[iOut++] = halo->partners[iIn];
		}
	    }
	    halo->nNeighbors = iOut;
	}
    }
    else {
	int nbrs=0;
	/* Compute the coords in the canonical 2-d mesh, which is row-major */
	if (halo->coords[0] < 0 || halo->coords[1] < 0) {
	    j = rank % halo->dims[1];
	    i = rank / halo->dims[1];
	    halo->coords[0] = i;
	    halo->coords[1] = j;
	}
	else {
	    i = halo->coords[0];
	    j = halo->coords[1];
	    if (verbose) {
		printf( "coords = (%d,%d) in [%d,%d]\n", i, j, 
			halo->dims[0], halo->dims[1] );
	    }
	}

	drank = coordToRank2( i, j-1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i, j+1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	halo->nNeighbors = nbrs;
    }

    FUNCEXIT("setNeighbors2");
    return 0;
}

int setNeighbors4( HaloElement *halo, int isPeriodic )
{
    int topoType, size, rank, i, j;

    FUNCENTER("setNeighbors4");

    MPI_Comm_rank( halo->comm, &rank );
    MPI_Topo_test( halo->comm, &topoType );
    if (topoType == MPI_CART) {
	MPI_Cart_coords( halo->comm, rank, 2, halo->coords );
	MPI_Cart_shift( halo->comm, 0, 1, 
			&halo->partners[0], &halo->partners[1] );
	MPI_Cart_shift( halo->comm, 1, 1, 
			&halo->partners[2], &halo->partners[3] );
	if (isPeriodic) 
	    halo->nNeighbors = 4;
	else {
	    /* Compress the number of neighbors */
	    int iIn = 0, iOut = 0;
	    for (iIn = 0; iIn < 4; iIn++) {
		if (halo->partners[iIn] != MPI_PROC_NULL) {
		    halo->partners[iOut++] = halo->partners[iIn];
		}
	    }
	    halo->nNeighbors = iOut;
	}
    }
    else {
	int drank, nbrs=0;
	MPI_Comm_size( halo->comm, &size );
	/* Compute the coords in the canonical 2-d mesh, which is row-major */
	/* FIXME: For the SMP-defined mesh, ensure that this method of 
	   calculating neighbors matches the assumptions of the SMP 
	   communicator */
	if (halo->coords[0] < 0 || halo->coords[1] < 0) {
	    j = rank % halo->dims[1];
	    i = rank / halo->dims[1];
	    halo->coords[0] = i;
	    halo->coords[1] = j;
	}
	else {
	    i = halo->coords[0];
	    j = halo->coords[1];
	    if (verbose) {
		printf( "coords = (%d,%d) in [%d,%d]\n", i, j, 
			halo->dims[0], halo->dims[1] );
	    }
	}

	drank = coordToRank2( i-1, j, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i+1, j, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i, j-1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i, j+1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;

	halo->nNeighbors = nbrs;
    }

    FUNCEXIT("setNeighbors4");
    return 0;
}

/* This sets the 8 neighbors in a 2-d mesh (including the diagonal elements); 
   note that the same amount of data is sent to each neighbor in these tests */
int setNeighbors8( HaloElement *halo, int isPeriodic )
{
    int topoType, size, rank, i, j;
    int coords[MAX_DIMS];

    FUNCENTER("setNeighbors8");

    MPI_Comm_rank( halo->comm, &rank );
    MPI_Topo_test( halo->comm, &topoType );
    if (topoType == MPI_CART) {
	MPI_Cart_coords( halo->comm, rank, 2, halo->coords );
	MPI_Cart_shift( halo->comm, 0, 1, 
			&halo->partners[0], &halo->partners[1] );
	MPI_Cart_shift( halo->comm, 1, 1, 
			&halo->partners[2], &halo->partners[3] );
	/* + + */
	for (i=0; i<2; i++) coords[i] = halo->coords[i];
	if (! isPeriodic && 
	    (coords[0] + 1 >= halo->dims[0]  || coords[1] + 1 >= halo->dims[1]))
	    halo->partners[4] = MPI_PROC_NULL;
	else {
	    coords[0] = (coords[0] + 1) % halo->dims[0];
	    coords[1] = (coords[1] + 1) % halo->dims[1];
	    MPI_Cart_rank( halo->comm, coords, &halo->partners[4] );
	}
	/* + - */
	for (i=0; i<2; i++) coords[i] = halo->coords[i];
	if (! isPeriodic &&
	    (coords[0] + 1 >= halo->dims[0] || coords[1] <= 0)) 
	    halo->partners[5] = MPI_PROC_NULL;
	else {
	    coords[0] = (coords[0] + 1) % halo->dims[0];
	    coords[1] = (coords[1] - 1 + halo->dims[1]) % halo->dims[1];
	    MPI_Cart_rank( halo->comm, coords, &halo->partners[5] );
	}
	/* - + */
	for (i=0; i<2; i++) coords[i] = halo->coords[i];
	if (! isPeriodic &&
	    (coords[0] <= 0 || coords[1] + 1 >= halo->dims[1])) 
	    halo->partners[6] = MPI_PROC_NULL;
	else {
	    coords[0] = (coords[0] - 1 + halo->dims[0]) % halo->dims[0];
	    coords[1] = (coords[1] + 1) % halo->dims[1];
	    MPI_Cart_rank( halo->comm, coords, &halo->partners[6] );
	}
	/* - - */
	for (i=0; i<2; i++) coords[i] = halo->coords[i];
	if (! isPeriodic &&
	    (coords[0] <= 0 || coords[1] <= 0)) 
	    halo->partners[7] = MPI_PROC_NULL;
	else {
	    coords[0] = (coords[0] - 1 + halo->dims[0]) % halo->dims[0];
	    coords[1] = (coords[1] - 1 + halo->dims[1]) % halo->dims[1];
	    MPI_Cart_rank( halo->comm, coords, &halo->partners[7] );
	}

	if (isPeriodic) 
	    halo->nNeighbors = 8;
	else {
	    /* Compress the number of neighbors */
	    int iIn = 0, iOut = 0;
	    for (iIn = 0; iIn < 8; iIn++) {
		if (halo->partners[iIn] != MPI_PROC_NULL) {
		    halo->partners[iOut++] = halo->partners[iIn];
		}
	    }
	    halo->nNeighbors = iOut;
	}
    }
    else {
	int drank, nbrs = 0;
	MPI_Comm_size( halo->comm, &size );
	/* Compute the coords in the canonical 2-d mesh, which is row-major */
	if (halo->coords[0] < 0 || halo->coords[1] < 0) {
	    j = rank % halo->dims[1];
	    i = rank / halo->dims[1];
	    halo->coords[0] = i;
	    halo->coords[1] = j;
	}
	else {
	    i = halo->coords[0];
	    j = halo->coords[1];
	    if (verbose) {
		printf( "coords = (%d,%d) in [%d,%d]\n", i, j, 
			halo->dims[0], halo->dims[1] );
	    }
	}

	drank = coordToRank2( i-1, j, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i+1, j, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i, j-1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i, j+1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i-1, j-1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i-1, j+1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i+1, j-1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;
	drank = coordToRank2( i+1, j+1, halo->dims, isPeriodic );
	if (drank != MPI_PROC_NULL) 
	    halo->partners[nbrs++] = drank;

	halo->nNeighbors = nbrs;
    }

    FUNCEXIT("setNeighbors8");
    return 0;
}

/* Create a new halo description from an existing one that has a single sender
   near the "center" of the mesh.
   
   Note that haloOut contains copies, not dups, of the buffers and the
   communicator.
*/
int MakeOnePoint( HaloElement *haloIn, HaloElement *haloOut )
{
    int i, j, k, masterRank, d;
    int iIn, iOut = 0;

    FUNCENTER("MakeOnePoint");

    /* Set this default */
    haloOut->isMaster = 0;

    /* First, determine the rank of the sender */
    i = haloIn->dims[0] / 2;
    j = haloIn->dims[1] / 2;
    d = -1;
    if (haloIn->coords[0] == i && haloIn->coords[1] == j) {
	MPI_Comm_rank( haloIn->comm, &d );
	haloOut->isMaster = 1;
    }
    MPI_Allreduce( &d, &masterRank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

    if (masterRank < 0) {
	printf( "Failed to find the origin!\n" );
	fflush(stdout);
	return 1;
    }

    /* Now, for any process that has that sender as a neighbor, leave them
       in the partner list.  Otherwise, remove and update nNeighbors */
    haloOut->nTest = haloIn->nTest;
    haloOut->comm  = haloIn->comm;
    haloOut->ndims = haloIn->ndims;
    for (k=0; k<haloIn->ndims; k++) {
	haloOut->dims[k]   = haloIn->dims[k];
	haloOut->coords[k] = haloIn->coords[k];
    }
    
    /* Compress the number of neighbors */
    if (haloOut->isMaster) {
	for (iIn = 0; iIn < haloIn->nNeighbors; iIn++) {
	    haloOut->sbuf[iIn]       = haloIn->sbuf[iIn];
	    haloOut->rbuf[iIn]       = haloIn->rbuf[iIn];
	    haloOut->partners[iIn] = haloIn->partners[iIn];
	}
	haloOut->nNeighbors = haloIn->nNeighbors;
    }
    else {
	for (iIn = 0; iIn < haloIn->nNeighbors; iIn++) {
	    if (haloIn->partners[iIn] == masterRank) {
		haloOut->sbuf[iOut]       = haloIn->sbuf[iIn];
		haloOut->rbuf[iOut]       = haloIn->rbuf[iIn];
		haloOut->partners[iOut++] = haloIn->partners[iIn];
	    }
	}
	haloOut->nNeighbors = iOut;
    }

    /* Note that one process sends and the others receive; this is 
       indicated by the isMaster flag */

    FUNCEXIT("MakeOnePoint");
    return 0;
}

/* ----------------------------------------------------------------------- */
/* Helper function: given an i, j coordinate, and a dims[] array, 
   compute the corresponding rank.
   If the coords are out of range in the non-periodic case, return 
   MPI_PROC_NULL.  Otherwise, return the cooresponding periodic value.
*/
int coordToRank2( int i, int j, int dims[], int isPeriodic )
{
    int rank = MPI_PROC_NULL;
    if (isPeriodic) {
	/* Make the coords fit */
	if (dims[0] > 0) {
	    while (i < 0) i += dims[0];
	    while (i >= dims[0]) i -= dims[0];
	}
	if (dims[1] > 0) {
	    while (j < 0) j += dims[1];
	    while (j >= dims[1]) j -= dims[1];
	}
	rank = i * dims[1] + j;
    }
    else {
	if (i >= 0 && i < dims[0] && j >= 0 && j < dims[1]) {
	    rank = i * dims[1] + j;
	}
	else 
	    rank = MPI_PROC_NULL;
    }
    return rank;
}

/* It is important to find the distribution of the measured times */
int FindTimes( double intime, 
	       double *maxTime, double *minTime, double *avgTime )
{
    double tavg;
    int size;

    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Allreduce( &intime, maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &intime, minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &intime, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    *avgTime = tavg / size;
    
    return 0;
}

/* This provides a quick check that the computed communication pattern is
   consistent; that is, that for every send there is a matching receive.
   Returns 1 on failure, 0 on success.
*/
int ConfirmPattern( HaloElement *halo )
{
    int csize, crank, i, *msgvec, nMsgs;
    MPI_Comm_size( halo->comm, &csize );
    MPI_Comm_rank( halo->comm, &crank );
    msgvec = (int *)calloc( csize, sizeof(int) );
    if (!msgvec) {
	fprintf( stderr, "Unable to allocate %d ints\n", csize );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    nMsgs = 0;
    for (i=0; i<halo->nNeighbors ; i++) {
	int p = halo->partners[i];
	if (p >= 0) {
	    if (p >= csize) {
		fprintf( stderr, "Invalid partner process %d for partner %d, must be in [0,%d]\n", p, i, csize-1 );
		MPI_Abort( MPI_COMM_WORLD, 1 );
	    }
	    msgvec[p]++;
	    nMsgs++;
	}
    }
    msgvec[crank] = -nMsgs;
    
    MPI_Allreduce( MPI_IN_PLACE, msgvec, csize, MPI_INT, MPI_SUM, halo->comm );
    nMsgs = 0;
    for (i=0; i<csize; i++) 
	if (msgvec[i] != 0) nMsgs++;
    free( msgvec );
    if (nMsgs != 0) {
	fprintf( stderr, "Found %d orphaned messages\n", nMsgs );
    }
    return nMsgs != 0;
}

/* Return the min and max values of v across all processes in communicator 
   comm.  Return the value at rank == commrank unless commanke < 0, in 
   which case return the values at all processes */
void IntFindMinMax( int v, int *vmin, int *vmax, MPI_Comm comm, int commrank )
{
    int tmp[2], rank;
    tmp[0] = v;
    tmp[1] = -v;
    if (commrank >= 0) {
	MPI_Comm_rank( comm, &rank );
	MPI_Reduce( (rank == commrank) ? MPI_IN_PLACE : tmp, tmp, 2, 
		    MPI_INT, MPI_MAX, commrank, comm );
	if (rank == commrank) {
	    *vmin = - tmp[1];
	    *vmax =   tmp[0];
	}
    }
    else {
	MPI_Allreduce( MPI_IN_PLACE, tmp, 2, MPI_INT, MPI_MAX, comm );
	*vmin = -tmp[1];
	*vmax =  tmp[0];
    }
}

/* */
void printMeshNodes( FILE *fp )
{
    MPI_Comm ncomm;
    int color, rank;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    color = (rank == myCliqueRanks[0]);
    MPI_Comm_split( MPI_COMM_WORLD, color, rank, &ncomm );
    if (color != 1) MPI_Comm_free( &ncomm );
    else {
        int size, *coordbuf, coords[5], qtorus[5], ndim=5, maxndim;
        MPI_Comm_size( ncomm, &size );
	topoMeshCoords(tinfo, &ndim, coords, qtorus);
	maxndim = ndim;
	MPI_Allreduce(MPI_IN_PLACE, &maxndim, 1, MPI_INT, MPI_MAX, ncomm);
/*
        TopoMyMeshCoords( coords+0, coords+1, coords+2 );
*/
        coordbuf = (int *)malloc( size * maxndim * sizeof(int) );
        MPI_Gather(coords, maxndim, MPI_INT,
		   coordbuf, maxndim, MPI_INT, 0, ncomm );
        MPI_Comm_rank(ncomm, &rank);
        if (rank == 0) {
	    int i, j, *c=coordbuf;
             for (i=0; i<size; i++) {
		 printf("[");
		 for (j=0; j<maxndim; j++) {
		     printf("%d%s", c[j], (j<maxndim-1)?",":"]\n");
		 }
		 c += maxndim;
             }
             fflush(stdout);
        }
        free( coordbuf );
        MPI_Comm_free( &ncomm );
    }
}
