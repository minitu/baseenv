#include "faststream.h"

/* ------------------------------------------------------------------------ */
int pf(const char *line, int len, void *d);
int pfwords(const char *line, int len, void *d);
void myAbort(int err, const char *msg);

int pf(const char *line, int len, void *d)
{
    /* static int lc=0;
       printf("Len = %d:", len);
       printf("%d:%s\n", lc++, line); */
    printf("%s\n", line); fflush(stdout);
    return MPI_SUCCESS;
}

#include <ctype.h>
static int nwords = 0;
int pfwords(const char *line, int len, void *d)
{
    const char *p = line;
    int   nleft = len;

    while (nleft > 0) {
	while (isspace(*p)) { p++; nleft--; if (nleft==0) break; }
	if (nleft > 0 && !isspace(*p)) {
	    nwords++;
	    do {
		p++; nleft--; if (nleft==0) break;
	    } while (!isspace(*p));
	}
    }
    return MPI_SUCCESS;
}

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
    TextFile fh=NULL;
    int i, err, wrank;
    char *filename = 0;
    int printex = 1;
    int verbose = 0;

    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-v") == 0) verbose = 1;
	else if (strcmp(argv[i], "-noprint") == 0) printex = 0;
	else if (filename == 0) filename = argv[i];
	else {
	    fprintf(stderr, "usage: %s [-v] [-noprint] filename\n", argv[0]);
	    return 1;
	}
    }

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    if (!filename) {
	fprintf(stderr, "usage: %s [-v] [-noprint] filename\n", argv[0]);
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (verbose) TextFileReadDebug(verbose);

    if (wrank == 0) {
	char method[128];
	int  len = 128;
	TextFileReadGetInfo(fh, method, &len);
	printf("Using TextFile method %s\n", method);
    }

    if (printex) {
	err = TextFileReadOpen(filename, MPI_COMM_WORLD, MPI_INFO_NULL, &fh);
	if (err) myAbort(err, "Open failed");
	err = TextFileReadIterate(fh, pf, NULL);
	if (err) myAbort(err, "Iteration failed");
	err = TextFileReadClose(fh);
	if (err) myAbort(err, "Close failed");
	err = TextFileReadPrintStats(fh, stdout);
	if (err) myAbort(err, "Print stats failed");
	err = TextFileReadFree(&fh);
    extern int errno;
	if (err) myAbort(err, "Free failed");
    }

    err = TextFileReadOpen(filename, MPI_COMM_WORLD, MPI_INFO_NULL, &fh);
    if (err) myAbort(err, "Open failed");
    err = TextFileReadIterate(fh, pfwords, NULL);
    if (err) myAbort(err, "Iteration failed");
    err = TextFileReadClose(fh);
    if (err) myAbort(err, "Close failed");
    err = TextFileReadPrintStats(fh, stdout);
    if (err) myAbort(err, "Print stats failed");
    err = TextFileReadFree(&fh);
    if (err) myAbort(err, "Free failed");

    /* The definition of MPI_Reduce makes MPI_IN_PLACE hard to use */
    MPI_Reduce((wrank==0) ? MPI_IN_PLACE: &nwords, &nwords, 1, MPI_INT,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    if (wrank == 0) {
	printf("Number of words = %d\n", nwords);
	fflush(stdout);
    }

    MPI_Finalize();

    return 0;
}
