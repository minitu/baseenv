/* -*- Mode: C; c-basic-offset:4 ; -*- */
#include <stdio.h>

#include "mpi.h"
#include "mpptest.h"
#include "getopts.h"
/* For sqrt */
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
/* Include string.h for memset */
#ifdef HAVE_STRING_H
#include <string.h>
#endif

/*
 * This file provides a communication pattern similar to "halo" or ghost
 * cell exchanges.
 *
 */

/* Local structure 
   
   Special notes:
   The "readRecv" is used to force the received data into the receiving process.
   This may be important on machines with shared memory and on with a lazy
   caching strategy.  
 */
#define MAX_PARTNERS 64
typedef enum { WAITALL, WAITANY } HaloWaitKind;
typedef struct _HaloData {
    int          n_partners;
    int          partners[MAX_PARTNERS];
    HaloWaitKind kind;
    int          flags;                 /* MPI_Win_xxx flags, as appropriate */
    int          debug_flag;            /* Provide printed debug information */
    int          readRecv;              /* Touch received data */
    /* These next allow us to reuse any allocated buffer */
    int          maxLen;
    void         *sbuffer, *rbuffer;
} HaloData;

/* This macro can be used to touch (at least) one item per cache line.  In
   the absence of an easy way to determine the cache line, we use 32 bytes,
   since few systems have shorter cache lines.  The call at the end
   (which never is made because of the other test) */
/* This is externally visible so that the compiler cannot know that the
   value is 0 */
int touchParm = 0;
#define TOUCH(_ptr,_bytes) {			\
	int _i, *_p = (int *)(_ptr), _r=0;      \
	const int _n = (_bytes) / sizeof(int);  \
	const int _stride = 32 / sizeof(int) ;  \
	if (_r) {\
	    for (_i=0; _i<_n; _i+=_stride ) {	\
		_r=_r+_p[_i];			\
	    }					\
	}					\
	if (touchParm && _r==0) { MPI_Abort(MPI_COMM_WORLD,1); }\
    }

static void halo_set_buffers( int len, HaloData *ctx, 
			      char *sbuffer[], char *rbuffer[] );
static void halo_free_buffers( HaloData *ctx, 
			       char *sbuffer[], char *rbuffer[] );
double halo_nb( int reps, int len, HaloData *ctx );
double halo_persistnb( int reps, int len, HaloData *ctx );

#ifdef HAVE_MPI_PUT
void haloGetRMABuffers( int, HaloData *, char **, char ** );
void haloFreeRMABuffers( HaloData * );

double halo_put( int reps, int len, HaloData *ctx );
#ifdef HAVE_MPI_GET
double halo_get( int reps, int len, HaloData *ctx );
#endif
#define FLAG_FENCE 0x1
#define FLAG_INFO  0x2
#define FLAG_ALLOC 0x4
#define FLAG_FREERMAMEM 0x40
#endif
#ifdef HAVE_MPI_PSCW
double halo_putpscw( int reps, int len, HaloData *ctx );
#define FLAG_NOCHECK 0x10
#endif
#ifdef HAVE_MPI_PASSIVERMA
double halo_putlock( int reps, int len, HaloData *ctx );
#ifdef HAVE_MPI_GET
double halo_getlock( int reps, int len, HaloData *ctx );
#endif
#define FLAG_LOCK_SHARED 0x8
#define FLAG_LOCK_NOBARRIER 0x20
#endif
int HaloSent( int, HaloData * );

/* Set up the initial buffers */
static void halo_set_buffers( int len, HaloData *ctx, 
			      char *sbuffer[], char *rbuffer[] )
{
    int  i;

    if (ctx->debug_flag) { 
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	printf( "[%d] len = %d, npartners = %d:",
		rank, len, ctx->n_partners );
	for (i=0; i<ctx->n_partners; i++) {
	    printf( ",%d", ctx->partners[i] );
	}
	puts( "" ); fflush( stdout );
    }

    /* Allocate send and receive buffers */
    if (len == 0) len += sizeof(int); 
    for (i=0; i<ctx->n_partners; i++) {
	sbuffer[i] = (char *)malloc( len );
	rbuffer[i] = (char *)malloc( len );
	if (!sbuffer[i] || !rbuffer[i]) {
	    fprintf( stderr, "Could not allocate %d bytes\n", len );
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	memset( sbuffer[i], -1, len );
	memset( rbuffer[i], 0, len );
    }
}
static void halo_free_buffers( HaloData *ctx, 
			       char *sbuffer[], char *rbuffer[] )
{
    int i;

    for (i=0; i<ctx->n_partners; i++) {
	free(sbuffer[i]);
	free(rbuffer[i]);
    }
}

double halo_nb( int reps, int len, HaloData *ctx )
{
    double elapsed_time;
    int    i, j, n_partners, n2_partners;
    double t0, t1;
    MPI_Request req[2*MAX_PARTNERS], *rq;
    MPI_Status  status[2*MAX_PARTNERS];
    char *(sbuffer[MAX_PARTNERS]), *(rbuffer[MAX_PARTNERS]);

    halo_set_buffers( len, ctx, sbuffer, rbuffer );

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    n2_partners  = 2 * n_partners; 
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    for(i=0;i<reps;i++){
	/*printf( "rep %d\n", i ); fflush(stdout);  */
	rq = req;
	for (j=0; j<n_partners; j++) {
	    MPI_Irecv( rbuffer[j], len, MPI_BYTE, ctx->partners[j], i,
		       MPI_COMM_WORLD, rq++ );
	    MPI_Isend( sbuffer[j], len, MPI_BYTE, ctx->partners[j], i, 
		       MPI_COMM_WORLD, rq++ );
	}
	if (ctx->kind == WAITALL)  {
	    MPI_Waitall( n2_partners, req, status );
	    if (ctx->readRecv) {
		for (j=0; j<n_partners; j++) {
		    TOUCH(rbuffer[j], len);
		}
	    }
	}
	else {
	    int idx;
	    for (j=0; j<n2_partners; j++) {
		MPI_Waitany( n2_partners, req, &idx, status );
		if (ctx->readRecv) {
		    if (! (idx & 0x1)) {
			/* Receives have even index */
			idx >>= 1;
			TOUCH(rbuffer[idx],len);
		    }
		}
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );
    halo_free_buffers( ctx, sbuffer, rbuffer );

    return(elapsed_time);
}

double halo_persistnb( int reps, int len, HaloData *ctx )
{
    double elapsed_time;
    int    i, j, n_partners, n2_partners;
    double t0, t1;
    MPI_Request req[2*MAX_PARTNERS], *rq;
    MPI_Status  status[2*MAX_PARTNERS];
    char *(sbuffer[MAX_PARTNERS]), *(rbuffer[MAX_PARTNERS]);

    halo_set_buffers( len, ctx, sbuffer, rbuffer );

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    n2_partners  = 2 * n_partners; 
    rq = req;
    for (j=0; j<n_partners; j++) {
	MPI_Recv_init( rbuffer[j], len, MPI_BYTE, ctx->partners[j], 0,
			 MPI_COMM_WORLD, rq++ );
	MPI_Send_init( sbuffer[j], len, MPI_BYTE, ctx->partners[j], 0, 
			 MPI_COMM_WORLD, rq++ );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    for(i=0;i<reps;i++){
	/*printf( "rep %d\n", i ); fflush(stdout);  */
	MPI_Startall( n2_partners, req );
	if (ctx->kind == WAITALL) {
	    MPI_Waitall( n2_partners, req, status );
	    if (ctx->readRecv) {
		for (j=0; j<n_partners; j++) {
		    TOUCH(rbuffer[j], len);
		}
	    }
	}
	else {
	    int idx;
	    for (j=0; j<n2_partners; j++) {
		MPI_Waitany( n2_partners, req, &idx, status );
		if (ctx->readRecv) {
		    if (! (idx & 0x1)) {
			/* Receives have even index */
			idx >>= 1;
			TOUCH(rbuffer[idx],len);
		    }
		}
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 
		   MPI_COMM_WORLD );
    halo_free_buffers( ctx, sbuffer, rbuffer );

    return(elapsed_time);
}

/*
 * There are a number of optimizations that we want to consider:
 *
 *  Win_fence asserts
 *  Win_create info (nolocks)
 *
 * These are selected with the -putoption string, where string is one of
 *    fence - Use fence asserts
 *    info  - Use nolocks info on create
 *    alloc - Use Alloc_mem instead of malloc
 * Multiple options may be used
 */
#ifdef HAVE_MPI_PUT
double halo_put( int reps, int len, HaloData *ctx )
{
    double  elapsed_time;
    int     i, j, n_partners;
    double  t0, t1;
    char    *sbuffer, *rbuffer;
    MPI_Win win;
    MPI_Aint offset;
    int      alloc_len;
    int      beginFence = 0, endFence=0;
    MPI_Info info;

    alloc_len = len * ctx->n_partners;
    if (alloc_len == 0) alloc_len = sizeof(double);

    if (ctx->debug_flag) {
	static int notcalled=1;
	if (notcalled) {
	    fprintf( stdout, "Flags = %x\n", ctx->flags ); fflush(stdout);
	    notcalled = 0;
	}
    }

    if (ctx->flags & FLAG_FENCE) {
	beginFence = MPI_MODE_NOPRECEDE;
	endFence   = MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED;
    }

    haloGetRMABuffers( alloc_len, ctx, &sbuffer, &rbuffer );

    if (ctx->flags & FLAG_INFO) {
	MPI_Info_create( &info );
	MPI_Info_set( info, "no_locks", "true" );
    }
    else {
	info = MPI_INFO_NULL;
    }
    MPI_Win_create( rbuffer, alloc_len, 1, info, MPI_COMM_WORLD, &win );
    if (info != MPI_INFO_NULL) {
	MPI_Info_free( &info );
    }

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    /* Note that another option is to use a single fence, since 
       fences separate access/exposure epochs. However, the easiest way
       to use the MPI RMA correctly is to begin/end with a fence, using
       the asserts to optimize for the actual operations */
    for(i=0;i<reps;i++){
	MPI_Win_fence( beginFence, win );
	/*printf( "rep %d\n", i ); fflush(stdout);  */
	offset = 0;
	for (j=0; j<n_partners; j++) {
	    if (ctx->partners[j] != MPI_PROC_NULL) {
		/* Fix for broken MPI implementations */
		MPI_Put( sbuffer+offset, len, MPI_BYTE, ctx->partners[j], 
			 offset, len, MPI_BYTE, win );
	    }
	    offset += len;
	}
	MPI_Win_fence( endFence, win );
	if (ctx->readRecv) {
	    for (j=0; j<n_partners; j++) {
		TOUCH(rbuffer+j*len, len);
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );

    MPI_Win_free( &win );
    if (ctx->flags & FLAG_FREERMAMEM) {
	haloFreeRMABuffers( ctx );
    }
    return(elapsed_time);
}

#ifdef HAVE_MPI_PSCW
/* Post-Start-Complete-Wait version of the halo exchange */
double halo_putpscw( int reps, int len, HaloData *ctx )
{
    double  elapsed_time;
    int     i, j, n_partners;
    double  t0, t1;
    char    *sbuffer, *rbuffer;
    MPI_Win win;
    MPI_Aint offset;
    int      alloc_len;
    int      postAssert = 0, startAssert = 0;
    MPI_Info info;
    MPI_Group wgroup;
    static MPI_Group group = MPI_GROUP_NULL;

    alloc_len = len * ctx->n_partners;
    if (alloc_len == 0) alloc_len = sizeof(double);

    if (ctx->debug_flag) {
	static int notcalled=1;
	if (notcalled) {
	    fprintf( stdout, "Flags = %x\n", ctx->flags ); fflush(stdout);
	    notcalled = 0;
	}
    }

    /* Get the communication group.  This uses a symmetric pattern.  We
       do this only once, since in this code, the same set of partners
       is used for every test */
    if (group == MPI_GROUP_NULL) {
	int      ranks[64], ng;
	MPI_Comm_group( MPI_COMM_WORLD, &wgroup );
	/* We must avoid MPI_PROC_NULL in the list of ranks, so 
	   we can't just use ctx->partners */
	ng = 0;
	for (i=0; i<ctx->n_partners; i++) {
	    if (ctx->partners[i] == MPI_PROC_NULL) continue;
	    if (ng >= (64-1)) {
		fprintf( stderr, "Too many neighbors for the put/pscw routine\n" );
		MPI_Abort( MPI_COMM_WORLD, 1 );
	    }
	    ranks[ng++] = ctx->partners[i];
	}
	MPI_Group_incl( wgroup, ng, ranks, &group );
	MPI_Group_free( &wgroup );
    }

    /* The only really useful assert for the PSCW case is NOCHECK, and
       that's really like a "ready send".  This requires more synchronization
       than is present in this current code */
    if (ctx->flags & FLAG_NOCHECK) {
	postAssert  = MPI_MODE_NOCHECK;
	startAssert = MPI_MODE_NOCHECK;
    }

    haloGetRMABuffers( alloc_len, ctx, &sbuffer, &rbuffer );

    if (ctx->flags & FLAG_INFO) {
	MPI_Info_create( &info );
	MPI_Info_set( info, "no_locks", "true" );
    }
    else {
	info = MPI_INFO_NULL;
    }
    MPI_Win_create( rbuffer, alloc_len, 1, info, MPI_COMM_WORLD, &win );
    if (info != MPI_INFO_NULL) {
	MPI_Info_free( &info );
    }

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    for(i=0;i<reps;i++){
	/*printf( "rep %d\n", i ); fflush(stdout);  */
	MPI_Win_post( group, postAssert, win );
	MPI_Win_start( group, startAssert, win );
	offset = 0;
	for (j=0; j<n_partners; j++) {
	    if (ctx->partners[j] != MPI_PROC_NULL) {
		/* Fix for broken MPI implementations */
		MPI_Put( sbuffer+offset, len, MPI_BYTE, ctx->partners[j], 
			 offset, len, MPI_BYTE, win );
	    }
	    offset += len;
	}
	MPI_Win_complete( win );
	MPI_Win_wait( win );
	if (ctx->readRecv) {
	    for (j=0; j<n_partners; j++) {
		TOUCH(rbuffer+j*len, len);
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );

    MPI_Win_free( &win );

    if (ctx->flags & FLAG_FREERMAMEM) {
	haloFreeRMABuffers( ctx );
    }
    return(elapsed_time);
}
#endif

/* This version uses locks instead of the fence synchronization */
double halo_putlock( int reps, int len, HaloData *ctx )
{
    double  elapsed_time;
    int     i, j, n_partners;
    double  t0, t1;
    char    *sbuffer, *rbuffer;
    MPI_Win win;
    MPI_Aint offset;
    int      alloc_len;
    int      lockType = MPI_LOCK_EXCLUSIVE;

    alloc_len = len * ctx->n_partners;
    if (alloc_len == 0) alloc_len = sizeof(double);

    if (ctx->debug_flag) {
	static int notcalled=1;
	if (notcalled) {
	    fprintf( stdout, "Flags = %x\n", ctx->flags ); fflush(stdout);
	    notcalled = 0;
	}
    }
    if (ctx->flags & FLAG_LOCK_SHARED) {
	lockType = MPI_LOCK_SHARED;
    }

    /* Passive requires that memory be allocated */
    ctx->flags |= FLAG_ALLOC;
    haloGetRMABuffers( alloc_len, ctx, &sbuffer, &rbuffer );

    MPI_Win_create( rbuffer, alloc_len, 1, MPI_INFO_NULL, 
		    MPI_COMM_WORLD, &win );

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    for(i=0;i<reps;i++){
        /*printf( "rep %d\n", i ); fflush(stdout);  */
	offset = 0;
	for (j=0; j<n_partners; j++) {
	    if (ctx->partners[j] != MPI_PROC_NULL) {
		/* Fix for broken MPI implementations */
		MPI_Win_lock( lockType, ctx->partners[j], 0, win );
		MPI_Put( sbuffer+offset, len, MPI_BYTE, ctx->partners[j], 
			 offset, len, MPI_BYTE, win );
		MPI_Win_unlock( ctx->partners[j], win );
	    }
	    offset += len;
	}
	/* To ensure that the data is available, we may need a barrier.
	   However, if other algorithmic reasons allow us to ensure that
	   the data we use is available, then this barrier is not
	   necessary.  We allow either choice in our testing. */
	if (! (ctx->flags & FLAG_LOCK_NOBARRIER)) {
	    MPI_Barrier( MPI_COMM_WORLD );
	    if (ctx->readRecv) {
		for (j=0; j<n_partners; j++) {
		    TOUCH(rbuffer+j*len, len);
		}
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );

    MPI_Win_free( &win );
    if (ctx->flags & FLAG_FREERMAMEM) {
	haloFreeRMABuffers( ctx );
    }
    return(elapsed_time);
}

/* Only bother with the MPI_Get tests if we also have put */
#ifdef HAVE_MPI_GET
double halo_get( int reps, int len, HaloData *ctx )
{
    double  elapsed_time;
    int     i, j, n_partners;
    double  t0, t1;
    char    *sbuffer, *rbuffer;
    MPI_Win win;
    MPI_Aint offset;
    int      alloc_len;
    int      beginFence = 0, endFence=0;
    MPI_Info info;

    alloc_len = len * ctx->n_partners;
    if (alloc_len == 0) alloc_len = sizeof(double);

    if (ctx->debug_flag) {
	static int notcalled=1;
	if (notcalled) {
	    fprintf( stdout, "Flags = %x\n", ctx->flags ); fflush(stdout);
	    notcalled = 0;
	}
    }

    if (ctx->flags & FLAG_FENCE) {
	beginFence = MPI_MODE_NOPRECEDE;
	endFence   = MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED;
    }

    haloGetRMABuffers( alloc_len, ctx, &sbuffer, &rbuffer );

    if (ctx->flags & FLAG_INFO) {
	MPI_Info_create( &info );
	MPI_Info_set( info, "no_locks", "true" );
    }
    else {
	info = MPI_INFO_NULL;
    }
    MPI_Win_create( rbuffer, alloc_len, 1, info, MPI_COMM_WORLD, &win );
    if (info != MPI_INFO_NULL) {
	MPI_Info_free( &info );
    }

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    /* Note that another option is to use a single fence, since 
       fences separate access/exposure epochs. However, the easiest way
       to use the MPI RMA correctly is to begin/end with a fence, using
       the asserts to optimize for the actual operations */
    for(i=0;i<reps;i++){
	MPI_Win_fence( beginFence, win );
	/*printf( "rep %d\n", i ); fflush(stdout);  */
	offset = 0;
	for (j=0; j<n_partners; j++) {
	    if (ctx->partners[j] != MPI_PROC_NULL) {
		/* Fix for broken MPI implementations */
		MPI_Get( sbuffer+offset, len, MPI_BYTE, ctx->partners[j], 
			 offset, len, MPI_BYTE, win );
	    }
	    offset += len;
	}
	MPI_Win_fence( endFence, win );
	if (ctx->readRecv) {
	    for (j=0; j<n_partners; j++) {
		TOUCH(rbuffer+j*len, len);
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );

    MPI_Win_free( &win );
    if (ctx->flags & FLAG_FREERMAMEM) {
	haloFreeRMABuffers( ctx );
    }
    return(elapsed_time);
}
/* This version uses locks instead of the fence synchronization */
double halo_getlock( int reps, int len, HaloData *ctx )
{
    double  elapsed_time;
    int     i, j, n_partners;
    double  t0, t1;
    char    *sbuffer, *rbuffer;
    MPI_Win win;
    MPI_Aint offset;
    int      alloc_len;
    int      lockType = MPI_LOCK_EXCLUSIVE;

    alloc_len = len * ctx->n_partners;
    if (alloc_len == 0) alloc_len = sizeof(double);

    if (ctx->debug_flag) {
	static int notcalled=1;
	if (notcalled) {
	    fprintf( stdout, "Flags = %x\n", ctx->flags ); fflush(stdout);
	    notcalled = 0;
	}
    }
    if (ctx->flags & FLAG_LOCK_SHARED) {
	lockType = MPI_LOCK_SHARED;
    }

    /* Passive requires that memory be allocated */
    ctx->flags |= FLAG_ALLOC;
    haloGetRMABuffers( alloc_len, ctx, &sbuffer, &rbuffer );

    MPI_Win_create( rbuffer, alloc_len, 1, MPI_INFO_NULL, 
		    MPI_COMM_WORLD, &win );

    elapsed_time = 0;
    n_partners   = ctx->n_partners;
    MPI_Barrier( MPI_COMM_WORLD );
    t0 = MPI_Wtime();
    for(i=0;i<reps;i++){
        /*printf( "rep %d\n", i ); fflush(stdout);  */
	offset = 0;
	for (j=0; j<n_partners; j++) {
	    if (ctx->partners[j] != MPI_PROC_NULL) {
		/* Fix for broken MPI implementations */
		MPI_Win_lock( lockType, ctx->partners[j], 0, win );
		MPI_Get( sbuffer+offset, len, MPI_BYTE, ctx->partners[j], 
			 offset, len, MPI_BYTE, win );
		MPI_Win_unlock( ctx->partners[j], win );
	    }
	    offset += len;
	}
	/* To ensure that the data is available, we may need a barrier.
	   However, if other algorithmic reasons allow us to ensure that
	   the data we use is available, then this barrier is not
	   necessary.  We allow either choice in our testing. */
	if (! (ctx->flags & FLAG_LOCK_NOBARRIER)) {
	    MPI_Barrier( MPI_COMM_WORLD );
	    if (ctx->readRecv) {
		for (j=0; j<n_partners; j++) {
		    TOUCH(rbuffer+j*len, len);
		}
	    }
	}
    }
    t1 = MPI_Wtime();
    elapsed_time = t1 - t0;
    /* Use the max since in the non-periodic case, not all processes have the
       same number of partners (and rate is scaled by max # or partners) */
    MPI_Allreduce( MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD );

    MPI_Win_free( &win );
    if (ctx->flags & FLAG_FREERMAMEM) {
	haloFreeRMABuffers( ctx );
    }
    return(elapsed_time);
}

#endif /* HAVE_MPI_GET */

#endif /* HAVE_MPI_PUT */

TimeFunction GetHaloFunction( int *argc_p, char *argv[], void *MsgCtx, 
			      char *title )
{
  HaloData *new;
  int      rank, size;
  int      i;
  int      is_periodic;

  new = (HaloData *)malloc( sizeof(HaloData) );
  if (!new) return 0;

  new->sbuffer    = 0;
  new->rbuffer    = 0;
  new->maxLen     = 0;
  new->flags      = 0;
  new->n_partners = 2;
  SYArgGetInt( argc_p, argv, 1, "-npartner", &new->n_partners );

  is_periodic = SYArgHasName( argc_p, argv, 1, "-periodic" );
  if (new->n_partners > MAX_PARTNERS) {
    fprintf( stderr, "Too many halo partners specified (%d); max is %d\n",
	     new->n_partners, MAX_PARTNERS );
  }
  /* Define the basic title (with number of partners) */
  sprintf( title, "Halo exchange (%d)", new->n_partners );
  if (is_periodic) {
      sprintf( title, "Halo exchange (%d-periodic mesh)", new->n_partners );
  }
  else {
      sprintf( title, "Halo exchange (%d)", new->n_partners );
  }


  new->debug_flag = SYArgHasName( argc_p, argv, 0, "-debug" );

  new->readRecv = SYArgHasName( argc_p, argv, 0, "-touch" );

  /* Compute partners.  We assume only exchanges here.  We use a simple rule
     to compute partners: we use
     rank (+-) 1, (+-)sqrt(size), (+-) sqrt(size) (+-)1, */

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  /* This will work as long as # partners is a multiple of 2 */
  if (new->n_partners > 1 && (new->n_partners & 0x1)) {
    fprintf( stderr, "Number of partners must be even\n" );
    return 0;
  }
  if (new->n_partners == 1) {
      if (rank & 0x1) new->partners[0] = rank - 1;
      else            new->partners[0] = rank + 1;
      if (is_periodic) {
	  if (new->partners[0] >= size) 
	      new->partners[0] -= size;
	  else if (new->partners[0] < 0) 
	      new->partners[0] += size;
      }
      else {
	  if (new->partners[0] >= size) 
	      new->partners[0] = MPI_PROC_NULL;
	  else if (new->partners[0] < 0)
	      new->partners[0] = MPI_PROC_NULL;
      }
  }
  else {
      int s1;
      /* First, load up partners with the distance to the partner */
      new->partners[0] = 1;
      new->partners[1] = -1;
      s1 = sqrt((double)size);
      new->partners[2] = s1;
      new->partners[3] = -s1;
      /* We must be careful to avoid multiple identical partners */
      if (s1 > 2) {
	  new->partners[4] = s1 + 1;
	  new->partners[5] = s1 - 1;
	  new->partners[6] = -s1 + 1;
	  new->partners[7] = -s1 - 1;
      }
      else if (s1 == 2) {
	  new->partners[4] = 3;
	  new->partners[5] = 4;
	  new->partners[6] = -4;
	  new->partners[7] = -3;;
      }

      for (i=0; i<new->n_partners; i++) {
	  if (is_periodic) 
	      new->partners[i] = (rank + new->partners[i] + size) % size;
	  else {
	      int partner;
	      partner = rank + new->partners[i];
	      if (partner >= size || partner < 0) partner = MPI_PROC_NULL;
	      new->partners[i] = partner;
	  }
      }
  }
  *(void **)MsgCtx = (void *)new;

#ifdef HAVE_MPI_PUT
  if (SYArgHasName( argc_p, argv, 1, "-put"  )) {
      char value[128];
      /* Get the options and set the halo flags */
      while (SYArgGetString( argc_p, argv, 1, "-putoption", 
			     value, sizeof(value))) {
	  if (strcmp( value, "fence" ) == 0) {
	      new->flags |= FLAG_FENCE;
	  }
	  else if (strcmp( value, "info" ) == 0) {
	      new->flags |= FLAG_INFO;
	  }
	  else if (strcmp( value, "alloc" ) == 0) {
	      new->flags |= FLAG_ALLOC;
	  }
	  else {
	      fprintf(stderr, "Unrecognized -putoption %s\n", value );
	  }
      }
      strcat( title, " -put/fence" );
      return (TimeFunction) halo_put;
  }
#else
  if (SYArgHasName( argc_p, argv, 1, "-put" )) {
      fprintf( stderr, "This MPI does not support MPI_Put\n" );
      MPI_Abort( MPI_COMM_WORLD, 1 );
      return 0;
  }
#endif
#ifdef HAVE_MPI_GET
  if (SYArgHasName( argc_p, argv, 1, "-get"  )) {
      char value[128];
      /* Get the options and set the halo flags */
      while (SYArgGetString( argc_p, argv, 1, "-putoption", 
			     value, sizeof(value))) {
	  if (strcmp( value, "fence" ) == 0) {
	      new->flags |= FLAG_FENCE;
	  }
	  else if (strcmp( value, "info" ) == 0) {
	      new->flags |= FLAG_INFO;
	  }
	  else if (strcmp( value, "alloc" ) == 0) {
	      new->flags |= FLAG_ALLOC;
	  }
	  else {
	      fprintf(stderr, "Unrecognized -putoption %s\n", value );
	  }
      }
      strcat( title, " -get/fence" );
      return (TimeFunction) halo_get;
  }
#else
  if (SYArgHasName( argc_p, argv, 1, "-get" )) {
      fprintf( stderr, "This MPI does not support MPI_Get\n" );
      MPI_Abort( MPI_COMM_WORLD, 1 );
      return 0;
  }
#endif
#ifdef HAVE_MPI_PSCW
  if (SYArgHasName( argc_p, argv, 1, "-putpscw"  )) {
      char value[128];
      /* Get the options and set the halo flags */
      while (SYArgGetString( argc_p, argv, 1, "-putoption", 
			     value, sizeof(value))) {
	  if (strcmp( value, "nocheck" ) == 0) {
	      new->flags |= FLAG_NOCHECK;
	  }
	  else if (strcmp( value, "info" ) == 0) {
	      new->flags |= FLAG_INFO;
	  }
	  else if (strcmp( value, "alloc" ) == 0) {
	      new->flags |= FLAG_ALLOC;
	  }
	  else {
	      fprintf(stderr, "Unrecognized -putoption %s\n", value );
	  }
      }
      strcat( title, " -putpscw" );
      return (TimeFunction) halo_putpscw;
  }
#else
  if (SYArgHasName( argc_p, argv, 1, "-putpscw" )) {
      fprintf( stderr, "This MPI does not support MPI_Win_post/start/complete/wait\n" );
      MPI_Abort( MPI_COMM_WORLD, 1 );
      return 0;
  }
#endif
#ifdef HAVE_MPI_PASSIVERMA
  else if (SYArgHasName( argc_p, argv, 1, "-putlock"  )) {
      char value[128];
      /* Get the options and set the halo flags */
      while (SYArgGetString( argc_p, argv, 1, "-putoption", 
			     value, sizeof(value))) {
	  if (strcmp( value, "shared" ) == 0) {
	      new->flags |= FLAG_LOCK_SHARED;
	  }
	  else if (strcmp( value, "nobarrier" ) == 0) {
	      new->flags |= FLAG_LOCK_NOBARRIER;
	  }
	  else {
	      fprintf(stderr, "Unrecognized -putoption %s\n", value );
	  }
      }
      strcat( title, " -put/lock" );
      return (TimeFunction) halo_putlock;
  }
  else if (SYArgHasName( argc_p, argv, 1, "-getlock"  )) {
      char value[128];
      /* Get the options and set the halo flags */
      while (SYArgGetString( argc_p, argv, 1, "-putoption", 
			     value, sizeof(value))) {
	  if (strcmp( value, "shared" ) == 0) {
	      new->flags |= FLAG_LOCK_SHARED;
	  }
	  else if (strcmp( value, "nobarrier" ) == 0) {
	      new->flags |= FLAG_LOCK_NOBARRIER;
	  }
	  else {
	      fprintf(stderr, "Unrecognized -putoption %s\n", value );
	  }
      }
      strcat( title, " -get/lock" );
      return (TimeFunction) halo_getlock;
  }
#else
  if (SYArgHasName( argc_p, argv, 1, "-putlock" )) {
      fprintf( stderr, "This MPI does not support MPI_Win_lock/unlock\n" );
      MPI_Abort( MPI_COMM_WORLD, 1 );
      return 0;
  }
#endif
  new->kind = WAITALL;
  if (SYArgHasName( argc_p, argv, 1, "-waitany" )) {
      new->kind = WAITANY;
      strcat( title, " - waitany" );
  }
  else {
      /* This is the default */
      strcat( title, " - waitall" );
  }
  if (SYArgHasName( argc_p, argv, 1, "-sendinit" )) {
      strcat( title, " (sendinit)" );
      return (TimeFunction) halo_persistnb;
  }
  return (TimeFunction) halo_nb;
}

/* This also needs a function to compute the amount of communication */
int HaloSent( int len, HaloData *ctx )
{
    int i, totlen = 0;

    /* This needs to be adjusted for partners that are PROC_NULL */
    for (i=0; i<ctx->n_partners; i++) {
	if (ctx->partners[i] != MPI_PROC_NULL) totlen += len;
    }
    return totlen;
}

/* This also needs a function to compute the amount of communication */
int GetHaloPartners( void *ctx )
{
    int i, totlen = 0;
    HaloData *hctx = (HaloData *)ctx;

    /* This needs to be adjusted for partners that are PROC_NULL */
    for (i=0; i<hctx->n_partners; i++) {
	if (hctx->partners[i] != MPI_PROC_NULL) totlen ++;
    }
    return totlen;
}

void PrintHaloHelp( void )
{
    fprintf( stderr, "\
   Special options for -halo:\n\
   -npartner n  - Specify the number of partners\n\
   -waitany     - Use a loop over waitany instead of a single waitall\n\
   -periodic    - Use periodic mesh partners\n\
   -debug       - Provide some debugging information\n\
   -touch       - Touch the data on the target end to check for cache effects\n\
                  on shared-memory systems\n\
" );
}

/* ---------------------------------------------------------------------- */
/* Provide common code for allocating memory used in Win_create.  To
 * work around extremely slow implementations of MPI_Alloc_mem, memory is
 * reused if possible (thus reducing the number of memory allocations)  */
/* ---------------------------------------------------------------------- */
#ifdef HAVE_MPI_PUT
void haloGetRMABuffers( int alloc_len, HaloData *ctx, 
			char **sbuffer_ptr, char **rbuffer_ptr )
{
    char *sbuffer, *rbuffer;

    /* Make alloc_len a little larger if it is too small.  In particular,
       avoid allocating a 0-byte buffer, as that sometime causes problems */
    if (alloc_len == 0) alloc_len = 64;
    /* If the buffer isn't large enough, reallocate it */
    if (ctx->maxLen < alloc_len) {
	haloFreeRMABuffers( ctx );
	
	/* Because some systems have costly alloc mem routines, allocate
	   at least twice the previous size */
	if (alloc_len < 2 *ctx->maxLen) 
	    alloc_len = 2 * ctx->maxLen;
#if defined(HAVE_SHMALLOC) && !defined(HAVE_MPI_ALLOC_MEM)
	sbuffer = (char *)shmalloc((unsigned)(alloc_len));
	rbuffer = (char *)shmalloc((unsigned)(alloc_len));
#else
	if (ctx->flags & FLAG_ALLOC) {
	    MPI_Alloc_mem(alloc_len,MPI_INFO_NULL,&sbuffer);
	    MPI_Alloc_mem(alloc_len,MPI_INFO_NULL,&rbuffer);
	}
	else {
	    sbuffer = (char *)malloc((unsigned)(alloc_len));
	    rbuffer = (char *)malloc((unsigned)(alloc_len));
	}
#endif
	if (!sbuffer || !rbuffer) {
	    fprintf( stderr, "Could not allocate %d bytes\n", alloc_len );
	    exit(1 );
	}
	ctx->sbuffer = sbuffer;
	ctx->rbuffer = rbuffer;
	memset( sbuffer, 0, alloc_len );
	memset( rbuffer, 0, alloc_len );
	ctx->maxLen  = alloc_len;
    }
    
    *sbuffer_ptr = ctx->sbuffer;
    *rbuffer_ptr = ctx->rbuffer;
}
void haloFreeRMABuffers( HaloData *ctx )
{
    if (ctx->maxLen > 0) {
	if (ctx->sbuffer) {
#if defined(HAVE_SHMALLOC) && !defined(HAVE_MPI_ALLOC_MEM)
	    shfree( ctx->sbuffer );
	    shfree( ctx->rbuffer );
#else
	    if (ctx->flags & FLAG_ALLOC) {
		MPI_Free_mem( ctx->sbuffer );
		MPI_Free_mem( ctx->rbuffer );
	    }
	    else {
		free( ctx->sbuffer );
		free( ctx->rbuffer );
	    }
#endif
	}
	ctx->sbuffer = 0;
	ctx->rbuffer = 0;
	ctx->maxLen = 0;
    }
}
#endif
