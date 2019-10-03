#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "mpi.h"
#define MPUTIL_PARALLEL
#include "mputil.h"

static int wrank = -1, wsize = -1;
static MPI_Comm mycomm = MPI_COMM_NULL;

static void outputtable(FILE *);
static void topoftable(void);
static void tablefree(void);

static int debug = 0;

void MPUTIL_Init(int thread_level)
{
    int requested = MPI_THREAD_SINGLE, provided;
    if (thread_level == 1) requested = MPI_THREAD_FUNNELED;
    else if (thread_level == 2) requested = MPI_THREAD_MULTIPLE;
    MPI_Init_thread(0, 0, requested, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    /* Default communicator */
    mycomm = MPI_COMM_WORLD;
}

void MPUTIL_Finalize(void)
{
    outputtable(stdout);
    tablefree();
    MPI_Finalize();
}

int MPUTIL_NumProcs(void)
{
    return wsize;
}

int MPUTIL_MyRank(void)
{
    return wrank;
}

int  MPUTIL_IsMaster(void)
{
    return wrank == 0;
}

static int cursize = 1, stride = 1;

int MPUTIL_IterBegin(void)
{
    if (debug && wrank == 0) {
	printf("Begining IterBegin for cursize %d\n", cursize);
    }

    if (cursize > wsize) return 0;
    MPI_Comm_split(MPI_COMM_WORLD,
		   wrank < cursize ? 1 : MPI_UNDEFINED, wrank,
		   &mycomm);
    topoftable();
    return 1;
}
void MPUTIL_IterEnd(void)
{
    if (debug && wrank == 0) {
	printf("Begining IterEnd for cursize %d\n", cursize);
    }

    if (mycomm != MPI_COMM_NULL && mycomm != MPI_COMM_WORLD)
	MPI_Comm_free(&mycomm);
    cursize += stride;
    MPI_Barrier(MPI_COMM_WORLD);
}
int MPUTIL_InIter(void)
{
    return (mycomm != MPI_COMM_NULL);
}
void MPUTIL_Sync(void)
{
    MPI_Barrier(mycomm);
}
int MPUTIL_Cursize(void)
{
    return cursize;
}

/* Combine data across processes */
double MPUTIL_GetDmin(const double t)
{
    double ttmp = t;
    /* In case this is called outside of the Iter block */
    if (mycomm != MPI_COMM_NULL) {
	MPI_Allreduce(MPI_IN_PLACE, &ttmp, 1, MPI_DOUBLE, MPI_MIN, mycomm);
    }
    return ttmp;
}

double MPUTIL_GetDmax(const double t)
{
    double ttmp = t;
    /* In case this is called outside of the Iter block */
    if (mycomm != MPI_COMM_NULL) {
	MPI_Allreduce(MPI_IN_PLACE, &ttmp, 1, MPI_DOUBLE, MPI_MAX, mycomm);
    }
    return ttmp;
}

int MPUTIL_SetDebug(int flag)
{
    int oldval = debug;
    debug = flag;
    return oldval;
}

void MPUTIL_Abort(int rc)
{
    MPI_Abort(MPI_COMM_WORLD, rc);
}
/* ------------------------------------------------------------------------- */
/* output table handling */
#define TMPBUF_INCR 512
typedef struct {
    char *line;
    int  maxline, curpos;
} textline;
static textline *outtable=0;
static int      currow = -1, maxrows=0, incrrows=16, incrline=128, ateol=1;
/* FIXME: More robust code will allow tmpbuf to grow */
static char     *tmpbuf = 0;
static int      tmpbuflen=0;

static void growtable(int);
static void growtablerow(int, int);
static void appendtext(const char *, int);

/* We can't use a va_list to pass variable args through a routine */
static void appendtext(const char *str, int slen)
{
    int     haseol, len, i;

    if (debug) {
	printf("%d: Appending %s\n", wrank, str);
    }
    /* Determine length to add; if necessary, expand row */
    haseol = 0;
    for (len=0; len < slen && str[len]; len++) {
	if (str[len] == '\n') haseol = 1;
    }

    /* if (currow < 0) printf("PANIC\n");*/
    if (outtable[currow].curpos + len > outtable[currow].maxline) {
	growtablerow(currow, len < incrline ? incrline : len);
    }

    for (i=0; str[i] && i<slen; i++) {
	/* Append, but leave out any newline at the end
	   Query: Convert newline into another character, e.g., tab or
	   comma */
	if (str[i] == '\n' && str[i+1] == 0) break;
	outtable[currow].line[outtable[currow].curpos+i] = str[i];
    }
    outtable[currow].curpos += i;
    outtable[currow].line[outtable[currow].curpos] = 0;

    if (debug)
	printf("%d: Done Appending\n", wrank);
}

void MPUTIL_OutLabel(const char *fmt, ...)
{
    va_list ap;
    int     actlen;

    /* We always assume that the label starts a new row */
    currow ++;
    ateol = 0;

    if (currow >= maxrows) {
	growtable(maxrows + incrrows);
    }

    /* The label is only generated for the first run and by rank 0 */
    if (cursize > 1 || wrank > 0) return;

    if (debug) {
	printf("%d: outlabel currow = %d %s\n",
	       wrank, currow, ateol?"At Eol":"In line");
    }

    if (tmpbuflen == 0) {
	tmpbuflen = TMPBUF_INCR;
	tmpbuf = (char *)malloc(tmpbuflen);
	if (!tmpbuf) abort();
    }
    do {
	va_start(ap, fmt);
	actlen = vsnprintf(tmpbuf, tmpbuflen, fmt, ap);
	va_end(ap);
	if (actlen < tmpbuflen) {
	    break;
	}
	tmpbuflen += TMPBUF_INCR;
	tmpbuf = (char *)realloc(tmpbuf, tmpbuflen);
	if (!tmpbuf) abort();
    } while (1);

    appendtext(tmpbuf, tmpbuflen);

    if (debug) {
	printf("%d: Done: outlabel currow = %d\n",
	       wrank, currow);
    }

}

void MPUTIL_OutApp(const char *fmt, ...)
{
    va_list ap;
    int     actlen;

    /* The output is only performed on rank 0 */
    if (wrank > 0) return;

    if (debug) {
	printf("%d: outapp currow = %d %s\n",
	       wrank, currow, ateol?"At Eol":"In line");
    }
    /* Only increment row if there was no label */
    if (ateol) {
	currow ++;
	ateol = 0;
    }

    if (currow >= maxrows) {
	growtable(maxrows + incrrows);
    }

    if (tmpbuflen == 0) {
	tmpbuflen = TMPBUF_INCR;
	tmpbuf = (char *)malloc(tmpbuflen);
	if (!tmpbuf) abort();
    }
    do {
	va_start(ap, fmt);
	actlen = vsnprintf(tmpbuf, tmpbuflen, fmt, ap);
	va_end(ap);
	if (actlen < tmpbuflen) {
	    break;
	}
	tmpbuflen += TMPBUF_INCR;
	tmpbuf = (char *)realloc(tmpbuf, tmpbuflen);
	if (!tmpbuf) abort();
    } while (1);

    appendtext(tmpbuf, tmpbuflen);

    if (debug) {
	printf("%d: Done: outapp currow = %d\n",
	       wrank, currow);
    }
}

static void growtable(int nsize)
{
    int      i;

    outtable = (textline *)realloc(outtable, sizeof(textline)*(maxrows+nsize));
    if (!outtable) abort();
    for (i=0; i<nsize; i++)
	outtable[maxrows+i].line = 0;
    maxrows += nsize;
}

static void growtablerow(int row, int incr)
{
    if (!outtable[row].line) {
	outtable[row].line = (char *)malloc(incr * sizeof(char));
	if (!outtable[row].line) abort();
	outtable[row].line[0] = 0;
	outtable[row].maxline = incr;
	outtable[row].curpos  = 0;
    }
    else {
	outtable[row].line = (char *)realloc(outtable[row].line,
					     outtable[row].maxline+incr);
	if (!outtable[row].line) abort();
	outtable[row].line[outtable[row].maxline] = 0;
	outtable[row].maxline += incr;
	/* curpos is left unchanged */
    }
}

static void outputtable(FILE *fp)
{
    if (wrank == 0) {
	int i;
	for (i=0; i<=currow; i++) {
	    fprintf(fp, "%s\n", outtable[i].line);
	}
	fflush(fp);
    }
}
static void topoftable(void)
{
    currow = -1;
    ateol  = 1;
}
static void tablefree(void)
{
    int i;
    for (i=0; i<maxrows; i++)
	if (outtable[i].line) free(outtable[i].line);
    maxrows = 0;
    free(outtable);
    if (tmpbuf) free(tmpbuf);
}
/* ------------------------------------------------------------------------- */
