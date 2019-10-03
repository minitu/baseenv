#include "faststream.h"

/* Test program */
int main(int argc, char *argv[])
{
    int      wsize, wrank;
    TextFile fh;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    /* Example writing to a file */
    TextFileWriteOpen(MPI_COMM_WORLD, 0, "testfile", MPI_INFO_NULL, &fh);
    /* We also test having one process perform no I/O */
    if (wrank != 4)
	TextFileWritePrintf(fh, "I am %d of %d\n", wrank, wsize);
    if (wrank == 2)
	TextFileWritePrintf(fh, "Rank 2 has an extra line of output\n");
    TextFileWriteFlush(fh);
    TextFileWriteClose(fh);

    /* Generate some performance output.  If running at large scale, just
       print summary. */
    if (wsize <= 16) {
	TextFileWritePrintStats(fh, stdout);
    }
    else {
	double     t[10], maxt[10], maxnoroot[10], z[10];
	MPI_Offset written;
	int        nt = 10;

	TextFileWriteGetStats(fh, &nt, t, &written);

	MPI_Reduce((wrank == 0) ? MPI_IN_PLACE : &written,
		   &written, 1, MPI_OFFSET, MPI_SUM,
		   0, MPI_COMM_WORLD);
	MPI_Reduce(t, maxt, nt, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (wrank == 0) {
	    int i;
	    for (i=0; i<nt; i++) z[i] = 0;
	}
	MPI_Reduce((wrank == 0) ? z : t, maxnoroot,
		   nt, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (wrank == 0) {
	    int i;
	    printf("Bytes written\t%ld\n", (long)written);
	    for (i=0; i<nt; i++) {
		int len = 40;
		char label[40];
		TextFileWriteNameStats(i, len, label);

		if (maxt[i] > 0) {
		    printf("%20s\t%.2e\n", label, maxt[i]);
		}
		if (maxnoroot[i] > 0 && maxnoroot[i] != maxt[i]) {
		    printf("root-%15s\t%.2e\n", label, maxnoroot[i]);
		}
	    }
	}
    }
    TextFileWriteFree(&fh);

    /* Example writing to stdout */
    TextFileWriteOpenFromFH(MPI_COMM_WORLD, 0, stdout, MPI_INFO_NULL, &fh);
    TextFileWritePrintf(fh, "I am %d of %d\n", wrank, wsize);
    if (wrank == 3)
	TextFileWritePrintf(fh, "Rank 3 has an extra line of output\n");
    TextFileWriteFlush(fh);
    TextFileWriteClose(fh);
    TextFileWriteFree(&fh);

    MPI_Finalize();
    return 0;
}
