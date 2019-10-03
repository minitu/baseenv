/*
 * Read and process text files with variable length records
 *
 * Operations:
 *    Open text file  - Set up buffering, including whether the entire file
 *                      should be preloaded
 *    Close text file - Close the file and free up an memory
 *    ProcessLines    - Allows multiple processes to work on parts of
 *                      a file
 *    GetLineByNumber - Return a particular line from the file
 *
 *  Internal or advanced routines
 *    Create index    - Create the index for the file lines, based on a
 *                      specified separator (e.g., \n or NULL)
 *    Set memory usage - Size of the internal memory to use for the file.
 *                       May replicate data to improve performance
 */

/*
 * ToDo:
 * 1. Add good test cases: large (enough) file, with distinct lines,
 * with a way to process and ensure all data read and processed.
 *
 * 3. Update stats to return data; separate routine to print.  Include
 * global data.
 *
 */

/* ------------------------------------------------------------------------- */
#include "faststream.h"

#define DEBUG 1

/* ---------------------------------------------------------------------- */
/* Change this to 4MiB for real use */
#ifdef DEBUG
static int defaultBlocksize = 1024;
#else
static int defaultBlocksize = 4194304; /* 4MiB */
#endif
/* Maximum to read is only 32GiB */
static MPI_Offset defaultMaxBytesToRead =
    (MPI_Offset)(1024*1024*1024)*(MPI_Offset)(32);
static int defaultLineLen   = 256;
static char defaultTermchar = '\n';
/* Verbose values > 0 cause extra printing.
 */
static int verbose          = 0;

typedef enum {
    TIME_OPEN=0, TIME_VIEW=1, TIME_FILE_IO=2, TIME_IO_WAIT=3,
    TIME_EDGE_SR=4, TIME_ITERATE=5, TIME_TOTAL=6, TIME_CLOSE=7,
    TIME_LAST=8 } TimeNames;

typedef enum { TEXTFILE_LINE_TAG=0, TEXTFILE_EOF_TAG=1} TagNames;

typedef struct TextFileData {
    MPI_File    mf;                  /* File handle for MPI IO */
    MPI_Comm    comm;                /* Private communicator for send/recv */
    MPI_Request frq,                 /* Used to read data in the background */
	prq;                         /* Used to update lines between blocks */
    MPI_Offset  nextOffset,          /* Location of the next read for this
					process */
	incrOffset,                  /* Increments between reads */
	bytesRead;                   /* Total bytes read by this process */
    /* Having encountered many errors in the MPI implementations, we
       include a panic value for bytes read */
    MPI_Offset  maxBytesToRead;      /* Abort if more than this much data
					is read */
    int         atEOF;               /* True if end of file seen (there may
					be data in curbuf) */
    int         sentEOF,             /* True if notified partner that EOF
					seen */
	recvEOF;                     /* True if notified that EOF seen */
    int         firstread;           /* Reading the very first block in the
					file on rank 0 */
    char        *curbuf;             /* Active block (double buffered with
					nextbuf) */
    char        *nextbuf;            /* Double buffer the data read from the
					file.  Each of size blocksize */
    int         firstterm, lastterm; /* Index of first and last terminator
					in curbuf */
    char        *firstline;          /* Special memory used for a complete
					first line of data for each block */
    char        *lastline;           /* Special memory used to hold the
					beginning of a line that will need
					to be combined with the next block */
    int         blocksize,           /* Size of block to read */
	curlen,                      /* Number of bytes in curbuf */
	firstlen,                    /* Number of bytes in firstline */
	lastlen,                     /* Number of bytes in lastline */
	firstAllocLen,               /* Size of firstline */
        lastAllocLen;                /* Size of lastline */

    int         nextproc, prevproc;  /* Ranks of next and previous processes,
				        in a ring */
    int         commrank;            /* My rank in the communicator */
    int         commsize;            /* Number of processes in comm */

    char        termchar;            /* char that terminates a line */
    /* Q: Add an escape char - so escapechar-termchar pairs do not terminate
       a line.  This permits things like <backslash><newline> to embed a
       <newline> */

    double      t[TIME_LAST];        /* Record performance information */
} TextFileData;

int TextFileiGetNextBlock(TextFile fh);
int TextFileiFindInterior(TextFile fh);
int TextFileiForwardLastLine(TextFile fh, int step);
int TextFileiProcessFirstLine(TextFile fh, int step);

#define DEBUG_MSG(_str) do { if (verbose) { printf("DBG:%s\n",_str);\
	    fflush(stdout);}}while(0)
#define DEBUG_MSG_I1(_fmt,_i1) do { if (verbose) { printf("DBG:"_fmt,_i1);\
	    fflush(stdout);}}while(0)
#define DEBUG_MSG_I2(_fmt,_i1,_i2) do { if (verbose) { \
	    printf("DBG:"_fmt,_i1,_i2); fflush(stdout);}}while(0)
#define DEBUG_MSG_I3(_fmt,_i1,_i2,_i3) do { if (verbose) { \
	    printf("DBG:"_fmt,_i1,_i2,_i3); fflush(stdout);}}while(0)

/* ---------------------------------------------------------------------- */
/* Need to make this a read-file - pass access-mode, and have future
   link for write, and disallow Read-write? */
/*@
  TextFileReadOpen - Open a file for reading and processing in parallel

  Input Parameters:
+ filename - File to open
. comm     - Communicator of processes that will open the file
- info     - MPI_Info value passed to MPI_File_open

  Output Paramter:
. fh_p     - TextFile handle to be used with other TextFile routines

  @*/
int TextFileReadOpen(const char filename[], MPI_Comm comm, MPI_Info info,
		     TextFile *fh_p)
{
    TextFile     fh = malloc(sizeof(TextFileData));
    char         cval[128];
    int          i, rank, size, err, flag, blocksize = defaultBlocksize;
    double       t;
#ifdef FILE_SET_VIEW_WORKS
    MPI_Offset   disp;
    MPI_Datatype dt, dt2;
    MPI_Aint     viewsize;
#endif

    /* Get the blocksize to use from the info */
    if (info != MPI_INFO_NULL) {
	MPI_Info_get(info, "TextFileBlocksize", sizeof(cval), cval, &flag);
	if (flag) {
	    /* Note: This should use strtol to detect errors in cval */
	    /* Note: I could check for a K, M (or Ki, Mi) at the end and
	       multiply by the appropriate value */
	    blocksize = atoi(cval);
	}
    }
    fh->blocksize = blocksize;

    fh->termchar = defaultTermchar;

    /* Zero the time records */
    for (i=0; i<TIME_LAST; i++) fh->t[i] = 0;

    /* Create a private communicator.  We might consider a 1-d Cartesian
       topology here */
    MPI_Comm_dup(comm, &fh->comm);

    /* Open the file.  MPI may return an error here */
    fh->t[TIME_TOTAL] = -MPI_Wtime();  /* We'll add Wtime at the end */
    t = MPI_Wtime();
    err = MPI_File_open(comm, filename,
			MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
			info, &(fh->mf));
    if (err) return err;
    fh->atEOF   = 0;
    fh->sentEOF = 0;
    fh->recvEOF = 0;
    fh->maxBytesToRead = defaultMaxBytesToRead;
    fh->t[TIME_OPEN] = MPI_Wtime() - t;

    /* A natural approach here is to create a view using MPI_File_set_view
       that exposes only the part of the file for this process.
       Unfortunately, as of MPICH 3.2, the transfer count set in the
       status return, which is the only way to detect and EOF, does not
       work; it always returns the requested read size rather than the
       number actually read */
    MPI_Comm_rank(fh->comm, &rank);
    MPI_Comm_size(fh->comm, &size);
    fh->commsize = size;
    fh->commrank = rank;
#ifdef FILE_SET_VIEW_WORKS
    /* Create a view for this process.
       If the viewsize is too large for an MPI_Aint, we could use
       MPI_File_read_at with an offset instead of using a view.
     */
    t = MPI_Wtime();
    disp = rank * (MPI_Offset)blocksize;
    MPI_Type_contiguous(blocksize, MPI_CHAR, &dt);
    viewsize = (MPI_Aint)size*(MPI_Aint)blocksize;
    if (viewsize < 0 || viewsize/size != blocksize) {
	/* The blocksize * number of processes must be less than the maximum
	   value of an MPI_Aint, typically 2GiB */
	return MPI_ERR_ARG;
    }
    MPI_Type_create_resized(dt, 0, viewsize, &dt2);
    MPI_Type_commit(&dt2);
    if (verbose) {
	printf("DBG: disp = %lld, viewsize = %ld\n",
	       (long long)disp, (long)viewsize);
    }
    err = MPI_File_set_view(fh->mf, disp, MPI_CHAR, dt2, "native",
			    MPI_INFO_NULL);
    if (err) return err;
    MPI_Type_free(&dt);
    MPI_Type_free(&dt2);
    fh->t[TIME_VIEW] = MPI_Wtime() - t;
#else
    fh->nextOffset = rank * (MPI_Offset)blocksize;
    fh->incrOffset = size * (MPI_Offset)blocksize;
#endif

    /* Allocate the buffers for double-buffered reads */
#ifdef HAVE_POSIX_MEMALIGN
    /* Using aligned data may improve I/O performance if MPI-IO uses RDMA
       to move data from the filesystem */
    posix_memalign(&fh->curbuf, ALIGN_SIZE, blocksize);
    posix_memalign(&fh->nextbuf, ALIGN_SIZE, blocksize);
#else
    fh->curbuf  = (char *)malloc(blocksize);
    fh->nextbuf = (char *)malloc(blocksize);
#endif
    if (!fh->curbuf || !fh->nextbuf) {
	return MPI_ERR_NO_SPACE;
    }

    /* Initiate read */
    fh->bytesRead = 0;
    t = MPI_Wtime();
#ifdef FILE_SET_VIEW_WORKS
    err = MPI_File_iread(fh->mf, fh->nextbuf, blocksize, MPI_CHAR, &fh->frq);
#else
#ifdef FILE_IREAD_WORKS
    err = MPI_File_iread_at(fh->mf, fh->nextOffset, fh->nextbuf, blocksize,
			    MPI_CHAR, &fh->frq);
    fh->nextOffset += fh->incrOffset;
#else
    fh->frq = MPI_REQUEST_NULL;  /* We'll perform the blocking read */
    err = MPI_SUCCESS;
#endif
#endif
    fh->t[TIME_FILE_IO] += MPI_Wtime() - t;
    if (err != MPI_SUCCESS) return err;
    /* We could initiate the second read, but that complicates some of this
       code */

    /* Note that no send of the last characters in the line has been
       initiated (another thing we could do here but adds complexity) */
    fh->prq = MPI_REQUEST_NULL;

    /* Preallocate first and lastline buffers */
    fh->firstAllocLen = defaultLineLen;
    fh->firstline     = (char *)malloc(fh->firstAllocLen);
    if (!fh->firstline) return MPI_ERR_NO_SPACE;
    if (rank == size-1) {
	/* Only the end of the ring needs a lastline buffer */
	fh->lastAllocLen = defaultLineLen;
	fh->lastline     = (char *)malloc(fh->lastAllocLen);
	if (!fh->lastline) return MPI_ERR_NO_SPACE;
    }
    else {
	fh->lastline     = 0;
	fh->lastAllocLen = 0;
    }

    /* Setup the ring of processes */
    fh->firstread = (rank == 0);
    fh->nextproc  = (rank + 1) % size;
    fh->prevproc  = (rank + size - 1) % size;

    *fh_p = fh;
    return err;
}

/*@
  TextFileReadClose - Close a text file

  Input Parameter:
. f - File handle from TextFileOpen

  @*/
int TextFileReadClose(TextFile fh)
{
    int      err;
    double   t;

    DEBUG_MSG("Text close");

    /* Add sanity check that no communication is pending */
    if (fh->frq != MPI_REQUEST_NULL || fh->prq != MPI_REQUEST_NULL) {
	if (verbose) {
	    fprintf(stderr, "Some communication not complete before close\n");
	}
	if (fh->frq != MPI_REQUEST_NULL)
	    MPI_Cancel(&fh->frq);
	if (fh->prq != MPI_REQUEST_NULL)
	    MPI_Cancel(&fh->prq);
    }

    fh->t[TIME_CLOSE] = -MPI_Wtime();
    err = MPI_File_close(&fh->mf);
    t = MPI_Wtime();
    fh->t[TIME_CLOSE] += t;
    fh->t[TIME_TOTAL] += t; /* Set to -Wtime in open */

    return err;
}

/*@
  TextFileReadFree - Free a text file and associated resources

  Input Parameter:
. f - File handle from TextFileOpen

  @*/
int TextFileReadFree(TextFile *fh_p)
{
    int      err=MPI_SUCCESS;
    TextFile fh = *fh_p;

    DEBUG_MSG("Text free");

    free(fh->curbuf);
    free(fh->nextbuf);
    if (fh->firstline) free(fh->firstline);
    if (fh->lastline)  free(fh->lastline);
    MPI_Comm_free(&fh->comm);
    free(fh);
    *fh_p = NULL;

    return err;
}

/*@
  TextFileReadIterate - Process all lines in a file

  Input Parameters:
+ fh - TextFile handle
. processline - Routine to process an input line
- extra_data - Pointer to extra data for processline

  @*/
int TextFileReadIterate(TextFile fh,
			int (*processline)(const char *, int len, void *),
			void *extra_data)
{
    int j, step;
    int rank;

    MPI_Comm_rank(fh->comm, &rank);
    step = 0;
    do {
	if (verbose == 2 && rank == 0) {
	    printf("."); fflush(stdout);
	}
	DEBUG_MSG("Starting processing of next block...");

	/* Get the next block */
	TextFileiGetNextBlock(fh);

	/* Find the interior */
	TextFileiFindInterior(fh);

	/* Forward the last partial line */
	TextFileiForwardLastLine(fh, step);

	/* Process the interior (if any) */
	/* Check that there IS an interior block */
	if (fh->firstterm < fh->lastterm) {
	    /* Process the interior of the block */
	    j = fh->firstterm + 1;
	    while (j < fh->lastterm) {
		int len;
		char *p = fh->curbuf + j;
		len = j;
		/* will terminate at lastterm */
		while (fh->curbuf[j] != fh->termchar) j++;
		fh->curbuf[j] = 0;    /* Null terminate the string */
		len = j - len;        /* Does not include the null */
		(*processline)(p, len, extra_data);
		j++;                  /* move to character after term char
					 location */
	    }
	}

	/* Receive and process the partial line at the beginning
	   of this block */
	TextFileiProcessFirstLine(fh, step);
	if (fh->firstlen > 0) {
	    (*processline)(fh->firstline, fh->firstlen, extra_data);
	}

    } while (!fh->sentEOF || !fh->recvEOF);

    /* Handle any pending EOF message */
    if (fh->prq != MPI_REQUEST_NULL && fh->sentEOF) {
	DEBUG_MSG("About to wait on prq");
	MPI_Wait(&fh->prq, MPI_STATUS_IGNORE);
    }

    if (verbose == 2 && rank == 0) {
	printf("\n"); fflush(stdout);
    }
    DEBUG_MSG("Finished reading file");

    return MPI_SUCCESS;
}

/*@
  TextFileReadDebug - Set debugging options for TextFile routines

Input Parameter:
. v - If nonzero, turns on verbose output to stdout

Notes:
This routine is for use by developers of TextFile and experts.  Understanding
the output requires access to the source code for TextFile.

  @*/
int TextFileReadDebug(int v)
{
    verbose = v;
    return MPI_SUCCESS;
}

int TextFileReadGetInfo(TextFile fh, char *str, int *len)
{
    int nlen;
#ifdef FILE_SET_VIEW_WORKS
    const char *method = "File set view";
#else
#ifdef FILE_IREAD_WORKS
    const char *method = "Iread_at";
#else
    const char *method = "Blocking read_at";
#endif
#endif
    nlen = strlen(method);
    if (nlen > *len) {
	strncpy(str, method, *len-1);
	str[*len] = 0;
    }
    else {
	strncpy(str, method, nlen+1);  /* Copy the null as well */
	*len = nlen;
    }

 /* Could return more information about the methods, such as allocation
    and buffer sizes */

    return MPI_SUCCESS;
}

/* ----------------------------------------------------------------------- */

int TextFileiGetNextBlock(TextFile fh)
{
    char       *swap;
    double     t;
    MPI_Status st;
    int        err = MPI_SUCCESS;

    DEBUG_MSG("GetNextBlock");
    /* Wait for the buffer from previous iread to be loaded */
    t = MPI_Wtime();
#ifdef FILE_IREAD_WORKS
    MPI_Wait(&fh->frq, &st);
#else
    err = MPI_File_read_at(fh->mf, fh->nextOffset, fh->nextbuf, fh->blocksize,
			    MPI_CHAR, &st);
    fh->nextOffset += fh->incrOffset;
#endif
    fh->t[TIME_IO_WAIT] += MPI_Wtime() - t;

    swap = fh->curbuf; fh->curbuf = fh->nextbuf; fh->nextbuf = swap;

    MPI_Get_count(&st, MPI_CHAR, &fh->curlen);
    DEBUG_MSG_I1("Read %d bytes from file\n", fh->curlen);
    fh->bytesRead += fh->curlen;
    if (fh->maxBytesToRead > 0 && fh->bytesRead > fh->maxBytesToRead) {
	if (verbose) {
	    fprintf(stderr, "Exceeded maximum number of bytes to read = %ld\n",
		    (long)fh->maxBytesToRead);
	    fflush(stderr);
	}
	return MPI_ERR_OTHER;
    }
    if (fh->curlen != fh->blocksize) fh->atEOF = 1;
    else {
	/* post the next buffer read */
	DEBUG_MSG("Post next buffer read");
	t = MPI_Wtime();
#ifdef FILE_SET_VIEW_WORKS
	err = MPI_File_iread(fh->mf, fh->nextbuf, fh->blocksize, MPI_CHAR,
			     &fh->frq);
#else
#ifdef FILE_IREAD_WORKS
	err = MPI_File_iread_at(fh->mf, fh->nextOffset, fh->nextbuf,
				fh->blocksize, MPI_CHAR, &fh->frq);
	fh->nextOffset += fh->incrOffset;
#else
	fh->frq = MPI_REQUEST_NULL;  /* We'll perform the blocking read */
#endif
#endif
	fh->t[TIME_FILE_IO] += MPI_Wtime() - t;
    }

    return err;
}

int TextFileiFindInterior(TextFile fh)
{
    char *p;
    int  i;

    /* Find the location of the first and last terminator */
    DEBUG_MSG("Find location of first and last terminator");
    fh->firstterm = -2;
    fh->lastterm  = -2;
    p             = fh->curbuf;
    /* Special case.  At the beginning of a file, there is a "virtual"
       terminator right before the block */
    if (fh->firstread) {
	fh->firstterm = -1;
    }
    else {
	for (i=0; i<fh->curlen; i++) {
	    if (p[i] == fh->termchar) {
		fh->firstterm = i;
		break;
	    }
	}
    }
    /* Special case.  At the end of a file, there is a "virtual"
       terminator right after the block. */
    if (fh->atEOF && fh->curlen > 0) {
	fh->lastterm = fh->curlen;
	if (fh->firstterm == -2) fh->firstterm = fh->lastterm;
    }
    else {
	for (i=fh->curlen-1; i>0; i--) {
	    if (p[i] == fh->termchar) {
		fh->lastterm = i;
		break;
	    }
	}
    }

    /* Sanity check: Must have found at least one terminator */
    if (fh->curlen > 0 && (fh->firstterm == -2 || fh->lastterm == -2)) {
	if (verbose) {
	    fprintf(stderr, "No terminator in block\n");
	    fflush(stderr);
	}
	return MPI_ERR_OTHER;
    }
    DEBUG_MSG_I3("first term index = %d, last = %d, len = %d\n",
	       fh->firstterm, fh->lastterm, fh->curlen);
    if (verbose) {
	char *b1 = 0;
	/* Print the leading and trailing partial lines */
	if (fh->firstterm > 0) {
	    b1 = (char *)malloc(fh->firstterm);
	    strncpy(b1, fh->curbuf, fh->firstterm-1);
	    b1[fh->firstterm-1] = 0;
	    printf("DBG: partial beginning = %s\n", b1);
	    free(b1);
	}
	if (fh->lastterm > fh->firstterm) {
	    int len = fh->curlen - fh->lastterm + 1;
	    b1 = (char *)malloc(len);
	    strncpy(b1, fh->curbuf+fh->lastterm+1, len-1);
	    b1[len-1] = 0;
	    printf("DBG: partial end = %s\n", b1);
	    free(b1);
	}
	fflush(stdout);
    }
    return MPI_SUCCESS;
}

int TextFileiForwardLastLine(TextFile fh, int step)
{
    double t;
    int    err=MPI_SUCCESS;

    /* Handle lines that cross blocks.  Different code is needed for the
       last process.  Note that if this process saw an EOF, then all of the
       processes after this one (e.g., the nextproc) will have also seen
       an EOF */
    if (!fh->sentEOF) {
	int  len, stag;
	char *linebuf;

	if (fh->curlen == 0)
	    len = 0;
	else
	    len = fh->curlen - fh->lastterm - 1;

	if (len > 0) {
	    if (fh->nextproc == 0) {
		/* This will be data for the next phase - copy it */
		if (fh->lastAllocLen < len) {
		    fh->lastAllocLen = len;
		    fh->lastline = (char *)realloc(fh->lastline, fh->lastAllocLen);
		    if (!fh->lastline) return MPI_ERR_NO_MEM;
		}
		fh->lastlen = len;
		strncpy(fh->lastline, fh->curbuf+fh->lastterm+1, len);
		linebuf = fh->lastline;
	    }
	    else {
		/* Note that we can't guarantee a null at the end */
		linebuf = fh->curbuf+fh->lastterm+1;
		DEBUG_MSG_I1("Sending %d from curbuf\n", len);
	    }
	}
	else {
	    /* It could be that there is no data to send (e.g., at EOF,
	       or a termchar at the very end of the block */
	    linebuf = NULL;
	    len     = 0;
	}
	if (fh->atEOF) {
	    stag        = TEXTFILE_EOF_TAG;
	    fh->sentEOF = 1;
	}
	else
	    stag = TEXTFILE_LINE_TAG;
	t = MPI_Wtime();
	MPI_Isend(linebuf, len, MPI_CHAR, fh->nextproc,
		  stag, fh->comm, &fh->prq);
	fh->t[TIME_EDGE_SR] += MPI_Wtime() - t;
	DEBUG_MSG_I3("Sent %d to %d with tag %d\n", len, fh->nextproc, stag);
    }
    else {
	/* SentEOF is true, so nothing more to send */
	fh->prq = MPI_REQUEST_NULL;
    }

    return err;
}

int TextFileiProcessFirstLine(TextFile fh, int step)
{
    int        len, err = MPI_SUCCESS;
    MPI_Status st;

    DEBUG_MSG("Entering ProcessFirstLine");

    /* Check to see if we're already done */
    if (fh->recvEOF) {
	fh->firstlen = 0;
	return MPI_SUCCESS;
    }
    if (fh->firstread) {
	/* The very first block read does not receive data (there is
	   no "previous" block) */
	fh->firstread = 0;
	return MPI_SUCCESS;
    }

    /* In all other cases, we will receive a block of data from the
       previous process.  In the case of rank 0, the data was sent
       from the previous step. */

    DEBUG_MSG_I1("About to probe from %d\n", fh->prevproc);
    MPI_Probe(fh->prevproc, MPI_ANY_TAG, fh->comm, &st);
    if (st.MPI_TAG == TEXTFILE_EOF_TAG) {
	fh->recvEOF = 1;
	DEBUG_MSG("Saw EOF tag");
    }
    MPI_Get_count(&st, MPI_CHAR, &len);
    DEBUG_MSG_I3("msg len = %d, total space needed = %d, available = %d\n",
		 len, len+fh->firstterm, fh->firstAllocLen);

    /* Ensure there is enough space to hold the message */
    if (len + fh->firstterm > fh->firstAllocLen) {
	DEBUG_MSG_I1("About to realloc %d\n", len+fh->firstterm);
	fh->firstline = (char *)realloc(fh->firstline, len + fh->firstterm);
	fh->firstAllocLen = len+fh->firstterm;
	if (!fh->firstline) return MPI_ERR_NO_MEM;
    }

    /* Receive the message and form the first line */
    MPI_Recv(fh->firstline, len, MPI_CHAR, fh->prevproc, st.MPI_TAG,
	     fh->comm, MPI_STATUS_IGNORE);
    if (len > 0) {
	fh->firstlen = len;
	if (fh->firstterm > 0) {
	    /* We may have no data on this process if we've seen and EOF */
	    strncpy(fh->firstline+len, fh->curbuf, fh->firstterm);
	    fh->firstlen += fh->firstterm;
	}
	fh->firstline[fh->firstlen] = 0;
	if (verbose) {
	    printf("DBG: Received %d for firstline of length %d\n", len,
		   fh->firstlen);
	    printf("DBG: firstline len = %d\n", (int)strlen(fh->firstline));
	    printf("DBG: firstline = %s\n", fh->firstline);
	    fflush(stdout);
	}
    }
    else {
	/* len == 0 and fh->firstterm > 0 can happen when a line
	   ends for the end of a block on the previous process */
	if (fh->firstterm > 0) {
	    strncpy(fh->firstline, fh->curbuf, fh->firstterm+1);
	    fh->firstlen = fh->firstterm;
	}
	else {
	    fh->firstlen     = 0;
	}
	fh->firstline[fh->firstlen] = 0;
    }

    /* Return the first line.  It will be processed in iterate */
    return err;
}

/*@ TextFileReadPrintStats - Print performance statistics about a TextFile

Input Parameters:
+ fh - TextFile for which statistics should be printed
- fp - File handle for output (stdout is common)

Notes:
This routine is collective (it computes some min and max values over all
processes in fh).  It makes no attempt to order the output; use
TextFileReadGetStats to access and control the output of this performance
data.

@*/
int TextFileReadPrintStats(TextFile fh, FILE *fp)
{
    double rate, t;
    long long b = (long long)fh->bytesRead;
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
    fprintf(fp, "Open      \t%.2e\n", fh->t[TIME_OPEN]);
#ifdef FILE_SET_VIEW_WORKS
    fprintf(fp, "View      \t%.2e\n", fh->t[TIME_VIEW]);
#endif
#ifdef FILE_IREAD_WORKS
    fprintf(fp, "Iread     \t%.2e\n", fh->t[TIME_FILE_IO]);
    fprintf(fp, "Wait Iread\t%.2e\n", fh->t[TIME_IO_WAIT]);
#else
    /* Use IO_WAIT for the blocking read */
    fprintf(fp, "Read      \t%.2e\n", fh->t[TIME_IO_WAIT]);
#endif
    fprintf(fp, "Block SR  \t%.2e\n", fh->t[TIME_EDGE_SR]);
    fprintf(fp, "Close     \t%.2e\n", fh->t[TIME_CLOSE]);

    return MPI_SUCCESS;
}

/*@
  TextFileReadGetStats - Return performance statistics about a TextFile

Input Parameters:
+ fh - TextFile
- tsize - Size of array 't'

Output Parameters:
+ t - Time values for different operations
- bytesRead - Number of bytes read by this process

Notes:
 This is not a collective routine - each process can read its own data
 independently.

  @*/
int TextFileReadGetStats(TextFile fh, int tsize, double t[],
			 MPI_Offset *bytesRead)
{
    int i, n;
    *bytesRead = fh->bytesRead;
    n = tsize;
    if (n > TIME_LAST) n = TIME_LAST;
    for (i=0; i<n; i++) t[i] = fh->t[i];
    return MPI_SUCCESS;
}

