/* Sequential execution code */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"

static int verbose = 0;

typedef struct {
  MPI_Comm lcomm;
  int      prevRank, nextRank;
} seqInfo;

static int seqKeyval = MPI_KEYVAL_INVALID;

static int seqDelFn( MPI_Comm comm, int keyval, void *attr, void *estate )
{
  seqInfo *myinfo = (seqInfo *)attr;

  /* To output something about the communicator, we need a common
     representation for a communicator. */
#if 0
  if (verbose) printf( "About to free communicator %d\n", myinfo->lcomm );
#endif

  MPI_Comm_free( &myinfo->lcomm );
  free( myinfo );
  return 0;
}

#if 0
/* We simply specify the NULL copy function when the keyval is created */
static int seqCopyFn( MPI_Comm comm, int keyval, void *attr, void *estate )
{
  seqInfo *myinfo, *oldinfo = (seqInfo *)attr;
  myinfo = (seqInfo *)malloc(sizeof(seqInfo));
  MPI_Comm_dup(oldinfo->lcomm, &myinfo->lcomm);
  myinfo->prevRank = oldinfo->prevRank;
  myinfo->nextRank = oldinfo->nextRank;
  /* Not done!! */
  return 0;
}
#endif
/*@
  seqBegin - Begin a sequential section

  Input Parameter:
. comm - Communicator of all processes involved in sequential section

  Notes:
  This routine is collective over the 'comm'; only one process at a time
  executes the code between the call to this routine and to 'seqEnd';
  the others wait.  The processes execute in rank order (in 'comm')

  @*/
void seqBegin( MPI_Comm comm )
{
  int      flag, mysize, myrank;
  seqInfo  *info;

  if (seqKeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval( MPI_COMM_NULL_COPY_FN, seqDelFn, &seqKeyval, NULL );
  }

  MPI_Comm_get_attr( comm, seqKeyval, &info, &flag );
  if (!flag) {
    info = (seqInfo *)malloc( sizeof(seqInfo) );
    MPI_Comm_dup( comm, &info->lcomm );
    MPI_Comm_rank( info->lcomm, &myrank );
    MPI_Comm_size( info->lcomm, &mysize );
    info->prevRank = myrank - 1;
    if (info->prevRank < 0)   info->prevRank = MPI_PROC_NULL;
    info->nextRank = myrank + 1;
    if (info->nextRank >= mysize) info->nextRank = MPI_PROC_NULL;
    if (verbose) {
      printf( "seqbegin: prev = %d, next = %d\n",
	      info->prevRank, info->nextRank );
    }
    MPI_Comm_set_attr( comm, seqKeyval, info );
  }
  MPI_Recv( NULL, 0, MPI_INT, info->prevRank, 0, info->lcomm,
	    MPI_STATUS_IGNORE );
}

/*@
  seqEnd - End a sequential section

  Input Parameter:
. comm - Communicator of all processes involved in sequential section

  @*/
void seqEnd( MPI_Comm comm )
{
  seqInfo *info;
  int     flag;

  /* Sanity check */
  if (seqKeyval == MPI_KEYVAL_INVALID)
    MPI_Abort( MPI_COMM_WORLD, 1 );
  MPI_Comm_get_attr( comm, seqKeyval, &info, &flag );
  if (!info || !flag)
    MPI_Abort( MPI_COMM_WORLD, 1 );
  if (verbose) {
    printf( "seqend: prev = %d, next = %d\n",
	    info->prevRank, info->nextRank );
  }
  MPI_Send( NULL, 0, MPI_INT, info->nextRank, 0, info->lcomm );

  /* Make everyone wait until all have completed their send */
  MPI_Barrier( info->lcomm );
}

/*@
  seqChangeOrder - Change the order in which processes proceed through a
  sequential section

  Input Parameters:
+ comm - The communicator containing all processes for the sequential section
. prev - The rank of the previous process
- next - The rank of the next process

  Notes:
  This routine changes the order in which the processes in 'comm' proceed in
  a sequential section begun with 'seqBegin'.  While not a collective routine,
  it is the users'' responsibility to ensure that the specification of 'prev'
  and 'next' contains on cycles and all processes are reachable.

  The process for which 'prev' is 'MPI_PROC_NULL' goes first; the last process
  must have 'MPI_PROC_NULL' for the 'next' value.
  @*/
void seqChangeOrder( MPI_Comm comm, int prev, int next )
{
  seqInfo *info;
  int flag;
  /* Sanity check */
  if (seqKeyval == MPI_KEYVAL_INVALID)
    MPI_Abort( MPI_COMM_WORLD, 1 );
  MPI_Comm_get_attr( comm, seqKeyval, &info, &flag );
  if (!info || !flag)
    MPI_Abort( MPI_COMM_WORLD, 1 );
  /* Update the order */
  info->prevRank = prev;
  info->nextRank = next;
}
