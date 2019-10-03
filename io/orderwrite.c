/*
 * Routines to provide ordered writing to stdout or to a file
 *
 * This first version simply forwards data to a single process that
 * then outputs the data.
 *
 * For output to a file that this software opens (rather than stdout, for
 * example), it could use parallel I/O, though the file sizes should rarely
 * need that.
 *
 */

#include "faststream.h"

/* -- end of public header -- */

#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

/* Many of these are not yet in use.  However, if a parallel write
   option is provided, these will be needed */
typedef enum {
    TIME_OPEN=0, TIME_VIEW=1, TIME_FILE_IO=2, TIME_IO_WAIT=3,
    TIME_EDGE_SR=4, TIME_ITERATE=5, TIME_TOTAL=6, TIME_CLOSE=7,
    TIME_LAST=8 } TimeNames;

typedef struct TextFileData {
    MPI_Comm comm;
    int      root, commrank, commsize;
    char     *lbuf;                    /* Local buffer */
    int      lbufAllocLen;
    int      lbufLen;
    char     *wbuf;                    /* Buffer used to accumulate writes, */
    int      wbufAllocLen;             /* valid on root only */
    int      wbufLen;                  /* Provided to allow buffering of
					  wbuf data.  UNUSED */
    int      fd;                       /* Unix for for output; valid at root */
    int      fdFromFH;                 /* True if fd provided */
    MPI_Offset bytesWritten;
    double   t[TIME_LAST];             /* Record performance information */
} TextFileData;

static int wbufAllocLenDefault = 65536;
static int lbufAllocLenDefault = 4096;

/* Verbose values > 0 cause extra printing.
 */
static int verbose          = 0;
#define DEBUG_MSG(_str) do { if (verbose) { printf("DBG:%s\n",_str);\
					    fflush(stdout);}}while(0)
#define DEBUG_MSG_I1(_fmt,_i1) do { if (verbose) { printf("DBG:"_fmt,_i1);\
						   fflush(stdout);}}while(0)
#define DEBUG_MSG_I2(_fmt,_i1,_i2) do { if (verbose) { \
	    printf("DBG:"_fmt,_i1,_i2); fflush(stdout);}}while(0)
#define DEBUG_MSG_I3(_fmt,_i1,_i2,_i3) do { if (verbose) { \
	    printf("DBG:"_fmt,_i1,_i2,_i3); fflush(stdout);}}while(0)

int TextFileiInitFH(MPI_Comm comm, int root, MPI_Info info, TextFile *fh_p);

/*@
  TextFileWriteOpen - Open a file for writing text output from a parallel
  program

  Input Parameters:
+ comm - MPI Communicator of all processes that will perform text writes
. root - Rank of process (in 'comm') that will perform I/O to the file
. fname - Name of the file to open.
- info - MPI Info object that may be used for opening the file.  Use
 'MPI_INFO_NULL' for no special behavior.

  Output Parameter:
. fh_p - Pointer to a 'TextFile'

  Notes:
  The 'TextFileWriteXxx' routines provide a way to ensure that multiple
  parallel processes correctly write output in a rank-ordered way.  Using
  output from individual processes, even if they are synchronized (e.g.,
  with 'MPI_Barrier' or the 'seqBegin' and 'seqEnd' routines), does not
  guarantee ordered output because the aggregation of data sent by I/O routines
  on each process to a file server does not guarantee any ordering.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.
  @*/
int TextFileWriteOpen(MPI_Comm comm, int root, const char *fname, MPI_Info info,
		      TextFile *fh_p)
{
    TextFile     fh;
    int          err;
    double       t;

    err = TextFileiInitFH(comm, root, info, fh_p);
    if (err) return err;

    /* Specific to each open style */
    fh = *fh_p;
    fh->bytesWritten = 0;
    fh->fdFromFH = 0;
    if (fh->commrank == fh->root) {
	/* Decide how to handle mask.  fopen uses the values below (0666) */
	fh->t[TIME_TOTAL] = -MPI_Wtime();  /* We'll add Wtime at the end */
	t = MPI_Wtime();
	fh->fd = open(fname, O_WRONLY | O_CREAT,
	      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
	if (fh->fd < 0) {
	    return MPI_ERR_OTHER;
	}
	fh->t[TIME_OPEN] = MPI_Wtime() - t;
    }

    *fh_p = fh;
    return MPI_SUCCESS;
}

/*@
  TextFileWriteOpenFromFH - Associate a TextFile with an open FILE handle

Input Parameters:
+ comm - MPI Communicator of all processes that will perform text writes
. root - Rank of process (in 'comm') that will perform I/O to the file
. fp   - Already open for writing FILE.  A common use is with 'stdout'
- info - MPI Info object that may be used for opening the file.  Use
 'MPI_INFO_NULL' for no special behavior.

  Output Parameter:
. fh_p - Pointer to a 'TextFile'

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.

  @*/
int TextFileWriteOpenFromFH(MPI_Comm comm, int root, FILE *fp, MPI_Info info,
			    TextFile *fh_p)
{
    TextFile     fh;
    int          err;

    err = TextFileiInitFH(comm, root, info, fh_p);
    if (err) return err;

    /* Specific to each open style */
    fh = *fh_p;
    if (fh->commrank == root) {
	fh->fd       = fileno(fp);
	fh->fdFromFH = 1;
    }

    return MPI_SUCCESS;
}

/*@
  TextFileWriteOrdered - Write ordered output to a TextFile

Input Parameters:
+ fh  - TextFile previously opened
. buf - Data to write
- len - Number of characters to write

Notes:
This routine writes to the file (or FILE handle) previously opened.  The
output is ordered with respect to the rank in the communicator associated
with the TextFile.  This is a collective routine; internally, it calls
'TextFileWriteFlush' to output data.

Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.
  @*/
int TextFileWriteOrdered(TextFile fh, const char *buf, size_t len)
{
    if (len > fh->lbufAllocLen - fh->lbufLen) {
	int reqlen = len + fh->lbufLen;
	fh->lbuf = (char *)realloc(fh->lbuf, reqlen);
	if (!fh->lbuf) return MPI_ERR_NO_MEM;
	fh->lbufAllocLen = reqlen;
        DEBUG_MSG_I1("WriteOrdered: Realloc %d\n",reqlen);
    }
    memcpy(fh->lbuf+fh->lbufLen, buf, len);
    fh->lbufLen += len;
    TextFileWriteFlush(fh);

    return MPI_SUCCESS;
}
/*@
  TextFileWritePrintf - Write ordered output to a file, using printf-style

Input Parameters:
+ fh  - TextFile previously opened
. format - 'printf' style format
- ... - Parameters for the 'format'.

Notes:
This routine writes to the file (or FILE handle) previously opened.  The
output is ordered with respect to the rank in the communicator associated
with the TextFile.  The output is buffered in this routine; a collective
call to 'TextFileWriteFlush' is requried to ensure that the data is output.
The usage is nearly identical to 'fprintf'; the only difference is that
instead of a 'FILE' as the first parameter, a 'TextFile' must be provided.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.

  @*/
int TextFileWritePrintf(TextFile fh, const char *restrict format, ...)
{
    int     n, availLen;;
    va_list ap;

    availLen = fh->lbufAllocLen - fh->lbufLen;
    if (availLen <= 0) {
	/* Make sure at least some space is available */
        int len = fh->lbufAllocLen + 128;
	fh->lbuf = (char *)realloc(fh->lbuf, len);
	if (!fh->lbuf) return MPI_ERR_NO_MEM;
	fh->lbufAllocLen = len;
        DEBUG_MSG_I1("WritePrintf: Realloc %d\n", len);
        availLen = fh->lbufAllocLen - fh->lbufLen;
    }

    /* va_start and va_end are fragile and need to be immediately before
       and after the routine that needs them. */
    va_start(ap, format);
    n = vsnprintf(fh->lbuf+fh->lbufLen, availLen, format, ap);
    va_end(ap);
    if (n + 1 > availLen) {
	int len = n + 1 + fh->lbufLen;
	/* The return value from vsnprintf does not include the null */
	DEBUG_MSG_I1("vsnprintf returned %d\n", n);
	fh->lbuf = (char *)realloc(fh->lbuf, len);
	if (!fh->lbuf) return MPI_ERR_NO_MEM;
	fh->lbufAllocLen = len;
        DEBUG_MSG_I1("WritePrintf: Realloc %d\n", len);
	availLen = fh->lbufAllocLen - fh->lbufLen;
	va_start(ap, format);
	n = vsnprintf(fh->lbuf+fh->lbufLen, availLen, format, ap);
	va_end(ap);
	DEBUG_MSG_I1("Second vsnprintf returned %d\n", n);
	if (n + 1 > availLen) {
	    printf("PANIC: vsnprintf failed! n = %d, avail = %d\n", n, availLen);
	}
    }
    fh->lbufLen += n;

    return MPI_SUCCESS;
}

/*
  Notes on TextFileWriteFlush.
  This is a very simple implementation where all processes other than the
  root (rank 0 by default) send their buffered output to the root process;
  that process writes out the data immediately as it is received.  This
  should be sufficient for modest output, even at scale.  However, the
  following updates may improve the performance.

  1. Buffering at the root.  This works for all output files.  Rather
  than immediately write the output that is received from another
  process, the root process buffers it until the next process would
  exceed the remaining buffer availability.  When that happens, the
  current buffer is written out, and the length is reset to zero.

  2. Parallel output by more processes.  This only works for output files
  that are either opened with MPI_File_open or for which POSIX lseek works
  (e.g., this does not apply to output to stdout).  In addition, this should
  only be necessary if the total file size is very large (at least 10s of
  gigabytes) In this case, at the time TextFileWriteFlush is called,

      a. All processes first perform an MPI_Exscan with the length of their
      local buffered data.  The total length is distributed to all processes.

      b. A subset of processes is determined either based on the total data
      length or on a predetermined set.  These processes will write the data.
      In addition, those processes will write block-aligned data, if possible,
      with the same size as a disk block or stripe.  This eliminates the need
      for read-modify-write to those blocks, and the required synchronization
      between server and clients (this may let the file system perform better).
      Normally, there should be no more than one process per file stripe or
      separate file block (having multiple processes write to the same file
      server will not provide useful parallelism).

      c. The processes that will write must collect data from the other
      processes.  This is done explicitly by communicating the data.  In
      principle, this could be done by setting up the File view for each
      process, or using MPI_File_iwrite_all (or the older, two-phase
      collective I/O), but implementations are unlikely to be efficient or
      optimized for this communication pattern.

      d. Using either an MPI_File_set_view with MPI_File_iwrite or
      MPI_File_iwrite_at (depending on whether the non-blocking iwrite
      correctly works with views, which is *not* true for most versions
      of MPICH), the selected processes (in step b above) write to the
      output file.

      e. Note that if TextFileWriteFlush is called multiple times, after the
      first time, the output from the process with rank 0 will (almost always)
      not be aligned at a start of a file block or stripe.  This is also true
      if there are other ways to write to the output file (e.g., to have
      a single process write a header).
 */

/*@
   TextFileWriteFlush - Flush all pending output to a TextFile

Input Parameter:
. fh - TextFile handle

Notes:
This is a collective routine that must be called by every process in the
communicator that was used to open the 'TextFile' handle 'fh'.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.

  @*/
int TextFileWriteFlush(TextFile fh)
{
    double t;

    if (fh->commrank == fh->root) {
	int i;
	/* Eventually, can buffer output into wbuf and write only
	   when full, using write */
	for (i=0; i<fh->commsize; i++) {
	    if (i == fh->root) {
		int n;
		t = MPI_Wtime();
		do {
		    n = write(fh->fd, fh->lbuf, fh->lbufLen);
		} while (n == -1 && errno == EINTR);
		if (n != fh->lbufLen) {
		    /* A better solution is to ensure all of the other
		       processes know that the write has failed. */
		    MPI_ERR_OTHER;
		}
		fh->lbufLen = 0;
		fh->bytesWritten += n;
		fh->t[TIME_FILE_IO] += MPI_Wtime() - t;
	    }
	    else {
		MPI_Status st;
		int        len, n;
		MPI_Probe(i, 0, fh->comm, &st);
		MPI_Get_count(&st, MPI_CHAR, &len);
		if (len > fh->wbufAllocLen) {
		    fh->wbuf         = (char *)realloc(fh->wbuf, len);
		    fh->wbufAllocLen = len;
                    DEBUG_MSG_I1("WriteFlush: Realloc %d\n", len);
		}
		MPI_Recv(fh->wbuf, len, MPI_CHAR, i, 0, fh->comm,
			 MPI_STATUS_IGNORE);
		t = MPI_Wtime();
		do {
		    n = write(fh->fd, fh->wbuf, len);
		} while (n == -1 && errno == EINTR);
		if (n != len) {
		    /* A better solution is to ensure all of the other
		       processes know that the write has failed. */
		    MPI_ERR_OTHER;
		}
		fh->wbufLen = 0;
		fh->bytesWritten += n;
		fh->t[TIME_FILE_IO] += MPI_Wtime() - t;
	    }
	}
    }
    else {
	/* Use synchronous send to avoid flooding root process.  This
	   relies on the MPI implementation using rendezvous for this,
	   which is the common choice */
	MPI_Ssend(fh->lbuf, fh->lbufLen, MPI_CHAR, fh->root, 0, fh->comm);
    }

    return MPI_SUCCESS;
}

/*@
  TextFileWriteClose - Close a TextFile

Input Parameter:
. fh - TextFile handle previously opened

Notes:
  In the case where the 'TextFile' was opened from an existing file handle,
  this routine does `not` close the file associated with the handle.  However,
  this routine should still be called to ensure that any other operations
  needed to close the 'TextFile' are completed.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.
  @*/
int TextFileWriteClose(TextFile fh)
{
    if (fh->root == fh->commrank) {
	double t;
	fh->t[TIME_CLOSE] = -MPI_Wtime();
	close(fh->fd);
	t = MPI_Wtime();
	fh->t[TIME_CLOSE] += t;
	fh->t[TIME_TOTAL] += t; /* Set to -Wtime in open */
    }
    return MPI_SUCCESS;
}

/*@
  TextFileWriteFree - Free all resources associated with a TextFile

Input Parameter:
. fh_p - Pointer to a 'TextFile'

Notes:
  This routine should be called after 'TextFileWriteClose'.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.

@*/
int TextFileWriteFree(TextFile *fh_p)
{
    TextFile fh = *fh_p;

    /* TODO: Error check - make sure that the TextFile is closed first */

    free(fh->lbuf);
    if (fh->wbuf) free(fh->wbuf);
    free(fh);
    *fh_p = 0;

    return MPI_SUCCESS;
}

/*
 * Common code to setup the TextFile for Writes.  This is an internal routine
 * and should not be called by users.
 */
int TextFileiInitFH(MPI_Comm comm, int root, MPI_Info info, TextFile *fh_p)
{
    TextFile     fh = malloc(sizeof(TextFileData));
    int          i;

    if (!fh) {
	/* Panic */
	return MPI_ERR_NO_MEM;
    }

    /* Zero the time records */
    for (i=0; i<TIME_LAST; i++) fh->t[i] = 0;

    MPI_Comm_dup(comm, &fh->comm);
    MPI_Comm_rank(fh->comm, &fh->commrank);
    MPI_Comm_size(fh->comm, &fh->commsize);
    fh->root         = root;

    if (fh->commrank == fh->root) {
	fh->wbufAllocLen = wbufAllocLenDefault;
	fh->wbuf         = (char *)malloc(fh->wbufAllocLen);
	fh->wbufLen      = 0;
	if (!fh->wbuf) return MPI_ERR_NO_MEM;
        DEBUG_MSG_I1("iInitFH: malloc %d\n",fh->wbufAllocLen);
    }
    else {
	fh->wbuf = 0;
    }
    fh->lbufAllocLen = lbufAllocLenDefault;
    fh->lbuf         = (char *)malloc(fh->lbufAllocLen);
    fh->lbufLen      = 0;
    if (!fh->lbuf) return MPI_ERR_NO_MEM;

    fh->fd = -1;

    *fh_p = fh;

    return MPI_SUCCESS;
}

/*@ TextFileWritePrintStats - Print performance statistics about a TextFile

Input Parameters:
+ fh - TextFile for which statistics should be printed
- fp - File handle for output (stdout is common)

Notes:
This routine is collective (it computes some min and max values over all
processes in fh).  It makes no attempt to order the output; use
TextFileWriteGetStats to access and control the output of this performance
data.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.
@*/
int TextFileWritePrintStats(TextFile fh, FILE *fp)
{
    double rate, t;
    long long b = (long long)fh->bytesWritten;
    double tvals[2];

    /* FIXME:
       Produce aggregate times over all processes, using fh->comm
    */
    t    = fh->t[TIME_FILE_IO] + fh->t[TIME_IO_WAIT];
    tvals[0] = t;
    tvals[1] = -t;
    MPI_Allreduce(MPI_IN_PLACE, tvals, 2, MPI_DOUBLE, MPI_MAX,
		  MPI_COMM_WORLD);
    tvals[1] = -tvals[1]; /* Min value */
    t = tvals[0];         /* Use maximum time for rate computation */
    if (t > 0) {
	rate = b / t;
	fprintf(fp, "Bytes\t%lld\t%.2e\t%.2e\n", b, t, rate);
    }
    else {
	fprintf(fp, "Bytes\t%lld\t0.0\t--\n", b);
    }
    fprintf(fp, "Open       \t%.2e\n", fh->t[TIME_OPEN]);
#if 0
#ifdef FILE_SET_VIEW_WORKS
    fprintf(fp, "View       \t%.2e\n", fh->t[TIME_VIEW]);
#endif
#ifdef FILE_IWRITE_WORKS
    fprintf(fp, "Iwrite     \t%.2e\n", fh->t[TIME_FILE_IO]);
    fprintf(fp, "Wait Iwrite\t%.2e\n", fh->t[TIME_IO_WAIT]);
#else
    /* Use IO_WAIT for the blocking write */
    fprintf(fp, "Write      \t%.2e\n", fh->t[TIME_IO_WAIT]);
#endif
#endif
    fprintf(fp, "Block SR   \t%.2e\n", fh->t[TIME_EDGE_SR]);

    return MPI_SUCCESS;
}

/*@
  TextFileWriteGetStats - Return performance statistics about a TextFile

Input Parameters:
+ fh - TextFile
- tsize - Size of array 't'

Output Parameters:
+ t - Time values for different operations
- bytesWritten - Number of bytes written by this process

Notes:
 This is not a collective routine - each process can read its own data
 independently.

  Return value:
  An error code.  Currently, 'MPI_SUCCESS' is returned on success and an
  MPI error class is returned on failure.  A future implementation may
  return a user-defined MPI error code.
  @*/
int TextFileWriteGetStats(TextFile fh, int *tsize, double t[],
			  MPI_Offset *bytesWritten)
{
    int i, n;
    *bytesWritten = fh->bytesWritten;
    n = *tsize;
    if (n > TIME_LAST) n = TIME_LAST;
    for (i=0; i<n; i++) t[i] = fh->t[i];
    *tsize = n;
    return MPI_SUCCESS;
}

int TextFileWriteNameStats(int idx, int len, char label[])
{
    switch (idx) {
    case TIME_OPEN:    strncpy(label, "Open", len); break;
    case TIME_VIEW:    strncpy(label, "View", len); break;
    case TIME_FILE_IO: strncpy(label, "Iwrite", len); break;
    case TIME_IO_WAIT:
#ifdef FILE_IWRITE_WORKS
	strncpy(label, "Wait Iwrite", len);
#else
	strncpy(label, "Write", len);
#endif
    case TIME_EDGE_SR: strncpy(label, "Block SR", len); break;
    case TIME_ITERATE: label[0] = 0; break;
    case TIME_TOTAL:   strncpy(label, "Total", len); break;
    case TIME_CLOSE:   strncpy(label, "Close", len); break;
    case TIME_LAST:    strncpy(label, "Dummy entry", len); break;
    default: return MPI_ERR_OTHER;
    }
    return MPI_SUCCESS;
}

