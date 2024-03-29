/* test program to find out how much buffering a system supplies */

#include "mpi.h"
#include "mpptestconf.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <stdio.h>

static int verbose = 0;

int TestBuffering( int bufsize, int myid, int other );

int main( int argc, char *argv[])
{
    int  myid, numprocs;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  bufsize, other, done, i;
    int  bufsizeHigh, bufsizeLow;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    if ((numprocs % 2) != 0) {
	fprintf( stderr, "buflimit requires an even number of processes\n" );
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Output processor names in rank order */
    MPI_Get_processor_name(processor_name,&namelen);
    if (myid > 0) 
	MPI_Recv( MPI_BOTTOM, 0, MPI_INT, myid - 1, 5, MPI_COMM_WORLD, 
		  &status );
    fprintf(stderr,"Process %d on %s\n", myid, processor_name);
    fflush( stderr );
    if (myid + 1 < numprocs)
	MPI_Send( MPI_BOTTOM, 0, MPI_INT, myid + 1, 5, MPI_COMM_WORLD );


    bufsize = 1024;
    other   = (myid + 1) % 2;
    done    = 0;

    /* 
     * This simply doubles bufsize until the communication appears delayed
     * until the send is issued.  An improvement is to ensure that the
     * bufsize at which this occurs is in fact the minimum size, and that the
     * result is stable (e.g., the implementation does not adapt to the
     * communication pattern).
     */

    /* Find a bufsize that blocks */
    while (!done && bufsize < 1024*1024*16) {
	done = TestBuffering( bufsize, myid, other );
	i = done;
	MPI_Allreduce( &i, &done, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	if (!done) bufsize *= 2;
    }
    bufsizeHigh = bufsize;
    bufsizeLow  = bufsize/2;
    bufsize     = bufsizeLow;
    while (bufsizeHigh - bufsizeLow > 4) {
	done = TestBuffering( bufsize, myid, other );
	i = done;
	MPI_Allreduce( &i, &done, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	if (done) {
	    bufsizeHigh = bufsize;
	}
	else {
	    bufsizeLow = bufsize;
	}
	bufsize = (bufsizeHigh + bufsizeLow) / 2;
    }
    if (myid == 0) 
	printf( "MPI_Send blocks with buffers of size %d\n", bufsizeHigh );

    MPI_Finalize();
    return 0;
}


int TestBuffering( int bufsize, int myid, int other )
{
    char *buf;
    double t1, t2, tbase;
    MPI_Status status;
    int done = 0;
            
    if ((buf = (char *) malloc (bufsize)) == NULL) {
	fprintf(stderr, "%d could not malloc %d bytes\n", myid, bufsize );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    /* fprintf(stderr,"%d sending %d to %d\n", myid, bufsize, other ); */
    if ((myid % 2) == 0) {
	MPI_Send( MPI_BOTTOM, 0, MPI_INT, other, 1, MPI_COMM_WORLD );
	MPI_Recv( MPI_BOTTOM, 0, MPI_INT, other, 2, MPI_COMM_WORLD, 
		  &status );
	/* Compute a time to send when the receive is waiting */
	t1 = MPI_Wtime();
	MPI_Send( buf, bufsize, MPI_CHAR, other, 100, MPI_COMM_WORLD );
	t2 = MPI_Wtime();
	tbase = t2 - t1;
	MPI_Recv( MPI_BOTTOM, 0, MPI_INT, other, 2, MPI_COMM_WORLD, 
		  &status );
	/* Compute a time when the receive is NOT waiting */
	t1 = MPI_Wtime();
	MPI_Send( buf, bufsize, MPI_CHAR, other, 100, MPI_COMM_WORLD );
	t2 = MPI_Wtime();
	if (t2 - t1 > 1.5 && t2 - t1 > 2.0 * tbase) {
	    if (verbose) printf( "MPI_Send blocks with buffers of size %d\n", 
		    bufsize );
	    done = 1;
	}
    }
    else {
	MPI_Recv( MPI_BOTTOM, 0, MPI_INT, other, 1, MPI_COMM_WORLD, 
		  &status );
	t1 = MPI_Wtime();
	MPI_Send( MPI_BOTTOM, 0, MPI_INT, other, 2, MPI_COMM_WORLD );
	MPI_Recv( buf, bufsize, MPI_CHAR, other, 100, MPI_COMM_WORLD, 
		  &status );
	MPI_Send( MPI_BOTTOM, 0, MPI_INT, other, 2, MPI_COMM_WORLD );
	while (MPI_Wtime() - t1 < 2.0) ;
	MPI_Recv( buf, bufsize, MPI_CHAR, other, 100, MPI_COMM_WORLD, 
		  &status );
    }
    if (verbose) fprintf(stderr,"%d received %d fr %d\n", myid, bufsize, other );
    free( buf );

    return done;
}
