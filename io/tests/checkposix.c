/*
 * This is a program to check one aspect of correct POSIX behavior.
 * The program uses MPI to manage multiple processes, which should be running
 * on separate nodes (not just cores) if possible.
 */

#include "baseenv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
/* Prepare to have lustre open as an option */
#ifdef HAVE_LUSTRE_LUSTREAPI_H
#include <lustre/lustreapi.h>
#endif

/*
 * Approach:
 * Each process either writes or read-modify-writes separate byte ranges
 * in the file (separated in either space or time, using MPI_Barrier to
 * separate in time).  The entire file is then checked by one process,
 * then by all of the other processes (concurrently).  Errors are summarized.
 *
 * Note that POSIX requires that this program succeed without errors, unless
 * there is no room to write the relatively small file.  A file system that
 * does not work with program cannot be called a POSIX file system, even if
 * it offers the POSIX API.  Semantics are important.
 *
 * To run this program:
 *
 * cd to a directory that writes to a parallel file system.  If the file
 * must be striped to get parallelism in the I/O operations, make sure that
 * the directory causes new files to be striped (there is no standard POSIX
 * call to open a striped file).
 *
 * mpiexec -n xx -ppn yy ./checkposix [-f filename] [-v] [-b bsize] [-g gap]
 *       [-n nb]
 *
 * Parameters for mpiexec
 *   -n xx -ppn yy - this syntax depends on the MPI implementation; this means
 *                   xx total processes and yy processes per node
 * Parameters for checkposix
 *   -f filenamem  - name of the file to use.  iotest.txt is the default.
 *                   This file should be in a parallel file system
 *   -v            - Turn on verbose messaging about activity
 *   -b bsize      - Read and write blocks of bsize integers.  Default is 128
 *   -g gap        - Most processes skip this many ints before the next block
 *                   Default is 64
 *   -n nb         - Number of major blocks that each process writes (see
 *                   below for an explanation of blocks and big blocks.
 *   Only if llapi_file_open is available:
 *   -l            - Use Lustre (and the llapi_file_open routine)
 *   -c c          - Use c stripes
 *   -s s          - Use a stripe size of s * 64K bytes (e.g., -s 16 gives
 *                   1MB stripes
 *
 * I/O pattern
 * Let there be wsize processes (size of MPI_COMM_WORLD).  This should be the
 * same as "xx" in the example above (e.g., mpiexec -n xx ...).
 *
 *   The file consists of nb consecutive runs of wsize blocks that consist
 *   of bsize + gap integers.
 *
 * Thus, the entire file is of size nb*wsize*(bsize+gap)*sizeof(int)
 *
 * Each process adds 1 to every element in the "bsize" blocks.  At the end of
 * the run, the code checks that the value for each of those elements is wsize.
 *
 * (bsize+gap)*sizeof(int) should be different than the natural blocksize of
 * the file system - the goal is to have multiple processes reading and
 * writing from the same disk block.  Thus bsize and gap should be relatively
 * small.  nb*(bsize+gap)*sizeof(int) should be larger than a single stripe.
 *
 */

void myAbort(int wrank, const char *msg);
void myAbortErrno(int wrank, const char *msg);
void clearbuf(int *buf, int len, int val);
long blocknum(off_t offset, int bsize, int gap);

static int maxErrors = 30;
static int verbose   = 0;

int main(int argc, char *argv[])
{
    const char *filename;
    int wsize, wrank, fd, rank, i, j, k, err, ww, nerr=0;
    int *buf, *gapbuf, bsize, gap, nb;
    double t, tread, twrite, tseek, t0, tbody;
    off_t  offset;
    FILE *outfile, *errfile;
#ifdef HAVE_LLAPI_FILE_OPEN
    int use_lustre = 0;
    unsigned long long stripesize=65536*16;
    int stripecount = 16;
#endif

    tread = twrite = tseek = 0.0;

    /* Set defaults */
    bsize    = 128;
    gap      = 64;
    nb       = 1024;
    filename = "iotest.txt";
    outfile  = stdout;
    errfile  = stderr;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Use to compute time-so-far in verbose output */
    t0 = MPI_Wtime();

    /* Process args */
    for (i=1; i<argc; i++) {
        if (strcmp(argv[i], "-f") == 0) {
            filename = argv[++i];
        }
        else if (strcmp(argv[i], "-b") == 0) {
            errno = 0;
            bsize = (int)strtol(argv[++i], NULL, 10);
            if (errno || (errno == 0 && bsize <= 0)) {
                myAbort(wrank, "-b argument invalid");
            }
        }
        else if (strcmp(argv[i], "-g") == 0) {
            errno = 0;
            gap   = (int)strtol(argv[++i], NULL, 10);
            if (errno || (errno == 0 && gap <= 0)) {
                myAbort(wrank, "-g argument invalid");
            }
        }
        else if (strcmp(argv[i], "-n") == 0) {
            errno = 0;
            nb = (int)strtol(argv[++i], NULL, 10);
            if (errno || (errno == 0 && nb <= 0)) {
                myAbort(wrank, "-n argument invalid");
            }
        }
        else if (strcmp(argv[i], "-v") == 0) {
            verbose = 1;
        }
#ifdef HAVE_LLAPI_FILE_OPEN
        else if (strcmp(argv[i], "-l") == 0) {
            use_lustre = 1;
        }
        else if (strcmp(argv[i], "-c") == 0) {
            errno = 0;
            stripecount = (int)strtol(argv[++i], NULL, 10);
            if (errno || (errno == 0 && nb <= 0)) {
                myAbort(wrank, "-c argument invalid");
            }
        }
        else if (strcmp(argv[i], "-s") == 0) {
            errno = 0;
            stripesize = 65536L * strtol(argv[++i], NULL, 10);
            if (errno || (errno == 0 && nb <= 0)) {
                myAbort(wrank, "-s argument invalid");
            }
        }
#endif
        else {
            myAbort(wrank, "invalid argument");
        }
    }


    buf    = (int *)malloc(bsize * sizeof(int));
    gapbuf = (int *)malloc(gap * sizeof(int));
    if (!buf || !gapbuf) myAbort(wrank, "Memory allocation");
    for (k=0; k<gap; k++) gapbuf[k] = wrank + k;

    if (verbose && wrank == 0) {
        fprintf(outfile, "filename %s, block size %d, gap %d, nblocks %d\n",
		filename, bsize, gap, nb);
#ifdef HAVE_LLAPI_FILE_OPEN
        if (use_lustre) {
            fprintf(outfile, "Lustre: stripecount %d, size %llu\n",
		    stripecount, stripesize);
        }
#endif
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
#ifdef HAVE_LLAPI_FILE_OPEN
    if (use_lustre) {
        /* Only one process can open with O_CREATE.  In fact, and
           contrary to the documentation, llapi_file_open *fails* if
           the file exists, even when O_CREAT is not specified. */
        if (wrank == 0) {
            fd = llapi_file_open(filename, O_RDWR | O_CREAT | O_TRUNC,
                                 S_IRUSR | S_IWUSR,
                                 stripesize, -1, stripecount, 0);
        if (fd < 0) myAbortErrno(wrank, "llapi_file_open for create");
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else {
            MPI_Barrier(MPI_COMM_WORLD);
            fd = open(filename, O_RDWR, S_IRUSR | S_IWUSR);
            if (fd < 0) myAbortErrno(wrank, "open in lustre branch for create");
        }
    }
    else
#endif
    /* An alternative is to have one rank create the file, barrier, and
       then have the others open the file with just O_RDWR */
    {
        fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
        if (fd < 0) myAbortErrno(wrank, "open for create");
    }

    /* Read-modify-write in the file */
    for (rank=0; rank<wsize; rank++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (verbose && wrank == 0) {
            fprintf(outfile, "RMW for outer step %d of %d starting at time %.2e\n",
		    rank, wsize, MPI_Wtime()-t0);
            fflush(outfile);
        }
	tbody = MPI_Wtime();
        for (j=0; j<nb; j++) {
            offset = (((wrank + rank) % wsize) + j * wsize) *
                     (bsize + gap)*sizeof(int);
            if (rank == 0) {
                for (k=0; k<bsize; k++)
                    buf[k] = 0;
            }
            else {
                t = MPI_Wtime();
                err = lseek(fd, offset, SEEK_SET);
                if (err == -1) myAbortErrno(wrank, "lseek for read");
                tseek += MPI_Wtime() - t;
		t = MPI_Wtime();
		err = read(fd, buf, bsize*sizeof(int));
		if (err == -1) myAbortErrno(wrank, "read");
		tread += MPI_Wtime() - t;
	    }
	    for (k=0; k<bsize; k++)
		buf[k] += 1;
	    t = MPI_Wtime();
	    err = lseek(fd, offset, SEEK_SET);
	    if (err == -1) myAbortErrno(wrank, "lseek for write");
	    tseek += MPI_Wtime() - t;
	    t = MPI_Wtime();
	    err = write(fd, buf, bsize*sizeof(int));
	    if (err == -1) myAbortErrno(wrank, "write");
	    twrite += MPI_Wtime() - t;
	    if (rank == 0) {
		t = MPI_Wtime();
		err = write(fd, gapbuf, gap*sizeof(int));
		if (err == -1) myAbortErrno(wrank, "write gap");
		twrite += MPI_Wtime() - t;
	    }
	}
	/* Provide some info about performance over this step */
	tbody = MPI_Wtime() - tbody;
	if (verbose && wrank == 0) {
	    double rate = (bsize + gap)*sizeof(int)*nb / tbody;
	    fprintf(outfile, "Rate %.2eMB/sec\n", rate * 1.0e-6);
	    fflush(outfile);
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (verbose && wrank == 0) {
	fprintf(outfile, "Finished writing; twrite = %.2e\n", twrite);
	fflush(outfile);
    }
    err = close(fd);
    if (err == -1) myAbortErrno(wrank, "close");
    if (verbose && wrank == 0) {
	fprintf(outfile, "Closed file\n");
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* Verify file */
    if (wrank == 0) {
	fd = open(filename, O_RDONLY);
	if (fd < 0) myAbortErrno(wrank, "open for check");
	for (ww=0; ww<wsize; ww++) {
            if (verbose) {
                fprintf(outfile, "Read check: ww=%d at time %.2e\n",
			ww, MPI_Wtime()-t0);
                fflush(outfile);
            }
	    for (rank=0; rank<wsize; rank++) {
		for (j=0; j<nb; j++) {
		    offset = (((ww + rank) % wsize) + j * wsize) *
			(bsize + gap)*sizeof(int);
		    err = lseek(fd, offset, SEEK_SET);
		    if (err == -1) myAbortErrno(wrank, "lseek for readcheck");
		    clearbuf(buf, bsize, -1);
		    err = read(fd, buf, bsize*sizeof(int));
		    if (err == -1) myAbortErrno(wrank, "readcheck");
		    clearbuf(gapbuf, gap, -1);
		    err = read(fd, gapbuf, gap*sizeof(int));
		    if (err == -1) myAbortErrno(wrank, "readcheck gap");
		    /* Verify data */
		    for (k=0; k<bsize; k++) {
			if (buf[k] != wsize) {
			    if (nerr < maxErrors) {
				fprintf(errfile, "%d: offset %ld (block %ld) buf[%d] = %d, expected %d\n",
					wrank, (long)offset,
					blocknum(offset, bsize, gap), k,
					buf[k], wsize);
				fflush(errfile);
			    }
			    nerr++;
			}
		    }
		    if (rank == 0) {
			for (k=0; k<gap; k++) {
			    if (gapbuf[k] != ww + k) {
				if (nerr < maxErrors) {
				    fprintf(errfile, "%d: offset %ld (block %ld) gap[%d] = %d, expected %d\n",
					    wrank, (long)offset,
					    blocknum(offset, bsize, gap), k,
					    gapbuf[k], ww+k);
				    fflush(errfile);
				}
				nerr++;
			    }
			}
		    }
		}
	    }
	}
	if (nerr != 0) {
	    fprintf(outfile, " Found %d errors.  File system is not POSIX\n", nerr);
	}
	else {
	    fprintf(outfile, " No errors.\n");
	}
	fprintf(outfile, "Time seek = %.2e, read = %.2e, write = %.2e\n",
		tseek, tread, twrite);
	fflush(outfile);
    }

    free(buf);
    free(gapbuf);

    MPI_Finalize();
    return 0;
}

void myAbort(int wrank, const char *msg)
{
    fprintf(stderr, "%d: %s\n", wrank, msg);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

/* Errno must not be explicitly declared (e.g., may use thread-local storage) */
#include <errno.h>
void myAbortErrno(int wrank, const char *msg)
{
    fprintf(stderr, "%d: %s: %s\n", wrank, msg, strerror(errno));
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

void clearbuf(int *buf, int len, int val)
{
    int i;
    for (i=0; i<len; i++) buf[i] = val;
}

long blocknum(off_t offset, int bsize, int gap)
{
    return offset / ((bsize+gap)*sizeof(int));
}
