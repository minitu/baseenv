/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include <stdio.h>

#include "mpi.h"
#include "mpptest.h"
#include "getopts.h"

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
/* Include string.h for memset */
#ifdef HAVE_STRING_H
#include <string.h>
#endif

/* This test generates some data on the overhead of operations that 
   may be involved in MPI RMA.  This was created as a result of running on 
   one system where MPI_Alloc_mem took over 1 second (!) */

#define MAX_SIZE 4096

int main( int argc, char *argv[] )
{
    int i, len, n, rank;
    double t, tall;
    MPI_Win win;
    void *buf;

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    /* First, check for catastrophically long alloc_mem time */
    len = MAX_SIZE;
    t = MPI_Wtime();
    MPI_Alloc_mem( len, MPI_INFO_NULL, &buf );
    t = MPI_Wtime() - t;
    MPI_Allreduce( &t, &tall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Free_mem( buf );

    n = 100 * MPI_Wtick() / tall;
    if (n > 1) {
	t = MPI_Wtime();
	for (i=0; i<n; i++) {
	    MPI_Alloc_mem( len, MPI_INFO_NULL, buf );
	    MPI_Free_mem( buf );
	}
	t = MPI_Wtime() - t;
	t = t / n;
	MPI_Allreduce( &t, &tall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    }

    if (rank == 0) {
	printf( "Alloc_mem\t%f\n", tall );
    }

    /* Do the same thing for Win_create */
    len = MAX_SIZE;
    MPI_Alloc_mem( len, MPI_INFO_NULL, &buf );
    MPI_Barrier( MPI_COMM_WORLD );
    t = MPI_Wtime();
    MPI_Win_create( buf, len, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win );
    MPI_Win_free( &win );
    t = MPI_Wtime() - t;
    MPI_Allreduce( &t, &tall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

    n = 100 * MPI_Wtick() / tall;
    if (n > 1) {
	t = MPI_Wtime();
	for (i=0; i<n; i++) {
	    MPI_Win_create( buf, len, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win );
	    MPI_Win_free( &win );
	}
	t = MPI_Wtime() - t;
	t = t / n;
	MPI_Allreduce( &t, &tall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    }
    MPI_Free_mem( buf );

    if (rank == 0) {
	printf( "Win_create\t%f\n", tall );
    }

    MPI_Finalize();

    return 0;
}
