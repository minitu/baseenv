#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

void myAbort(int, const char *);
void myAbort(int err, const char *msg)
{
    char errmsg[MPI_MAX_ERROR_STRING];
    int rlen;
    MPI_Error_string(err, errmsg, &rlen);
    fprintf(stderr, "%s: MPI error message = %s\n", msg, errmsg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char *argv[])
{
    int         err, outlen;
    MPI_File    fh;
    MPI_Request rq;
    MPI_Status  st;
    char        buf[1024];
    int         buflen = 1024;
    int         rank, size;
    MPI_Offset  disp;
    char        filename[] = "testfile";
    int         errs = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Must run with at least two processes */
    if (size < 2) {
	/* There is only one process */
	fprintf(stderr, "%s requires at least 2 processes\n", argv[0]);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }
    /* Create a test file of known length */
    /* First, make sure that the testfile does not already exist */
    if (rank == 0) {
	err = MPI_File_delete(filename, MPI_INFO_NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    err = MPI_File_open(MPI_COMM_WORLD, filename,
			MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
			MPI_INFO_NULL, &fh);
    if (err) myAbort(err, "Open for writing");
    if (rank == 0) {
	err = MPI_File_write(fh, buf, buflen, MPI_CHAR, &st);
	if (err) myAbort(err, "First write");
	err = MPI_File_write(fh, buf, buflen/2, MPI_CHAR, &st);
	if (err) myAbort(err, "Second write");
    }
    err = MPI_File_close(&fh);
    if (err) myAbort(err, "Close after writing");

    MPI_Barrier(MPI_COMM_WORLD);

    /* Open the file for reading */
    err = MPI_File_open(MPI_COMM_WORLD, filename,
			MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
			MPI_INFO_NULL, &fh);
    if (err) { fprintf(stderr, "Open failed with err = %d\n", err);
	MPI_Abort(MPI_COMM_WORLD,1);
    }

    disp = rank * buflen;

    /* With at least two processes, each process should read as follows:
       rank 0: buflen
       rank 1: buflen / 2
       rank > 1: 0
       Anything else is an error */
/* Note that the test fails with both nonblocking reads on Blue Waters */
#ifndef BLOCKING
    err = MPI_File_iread_at(fh, disp, buf, buflen, MPI_CHAR, &rq);
    if (err) myAbort(err, "Iread");
    MPI_Wait(&rq, &st);
#else
    err = MPI_File_read_at(fh, disp, buf, buflen, MPI_CHAR, &st);
    if (err) myAbort(err, "read");
#endif
    MPI_Get_count(&st, MPI_CHAR, &outlen);
    switch (rank) {
    case 0:
	if (outlen != buflen) {
	    errs++;
	    fprintf(stderr, "Expected %d, read %d\n", buflen, outlen);
	}
	break;
    case 1:
	if (outlen != buflen/2) {
	    errs++;
	    fprintf(stderr, "Expected %d, read %d\n", buflen/2, outlen);
	}
	break;
    default:
	if (outlen != 0) {
	    errs++;
	    fprintf(stderr, "Expected %d, read %d\n", 0, outlen);
	}
	break;
    }

    err = MPI_File_close(&fh);
    if (err) myAbort(err, "Close for read");
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
	err = MPI_File_delete(filename, MPI_INFO_NULL);
	if (err) myAbort(err, "Delete");
    }

    MPI_Allreduce(MPI_IN_PLACE, &errs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
	if (errs == 0) printf(" No Errors\n");
	else           printf(" Found %d errors\n", errs);
    }
    MPI_Finalize();
    return 0;
}
