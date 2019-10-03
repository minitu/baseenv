#ifndef HAVE_MPUTIL_H
#define HAVE_MPUTIL_H 1

#ifdef MPUTIL_PARALLEL

void MPUTIL_Init(int);
void MPUTIL_Finalize(void);
void MPUTIL_Abort(int);
int  MPUTIL_IsMaster(void);
double MPUTIL_GetDmin(const double t);
double MPUTIL_GetDmax(const double t);
int MPUTIL_IterBegin(void);
void MPUTIL_IterEnd(void);
int MPUTIL_InIter(void);
int MPUTIL_Cursize(void);
void MPUTIL_OutLabel(const char *fmt, ...);
void MPUTIL_OutApp(const char *fmt, ...);
void MPUTIL_Sync(void);
int MPUTIL_NumProcs(void);
int MPUTIL_MyRank(void);
int MPUTIL_SetDebug(int);

/*M
  MPUTIL_INIT - Initialize the MPUTIL timing environment

Input Parameter:
. thread_level - Indicate the desired level of threadedness.  See below.

Notes:
 When benchmarking codes on SMP nodes, it is common to run the benchmark
 on a single core. However, in many applications, this is not the configuration
 that is important for the application.  It is common for the application to
 be running the same (or similar) code on all cores at the same time.
 The 'MPUTIL_XXX' macros provide a relatively easy way to run the same,
 single-core benchmark on 1, 2, 3, ..., k cores for a node with k total cores.
 Additional macros makes it easy to generate simple tables where there is
 a separate column for each number of active cores for the benchmark.
 See 'MPUTIL_LABEL' and 'MPUTIL_OUTAPP' for output.

 A typical program that uses these macros will have the following structure\:
.vb
   ...
   MPUTIL_INIT(0);
   ...any code that is executed only once, independent of the number of cores
   MPUTIL_BEGIN;
   ...any code that each active core must execute, such as initialization
   MPUTIL_LABEL("label text for row, in printf format");
   MPUTIL_SYNC;
   ... benchmark code, including timer calls.  Total time in tval
   MPUTIL_OUTAPP("\t%.2e\n",tval);
   MPUTIL_END;
   MPUTIL_FINALIZE;
.ve

 The 'thread_level' is used only when creating the parallel version, which
 uses MPI processes to execute the benchmark in parallel.  To simplify the
 interface, this does not use the MPI-defined levels (thus, the benchmark
 never need inlude mpi.h).  The valid values are '0' for no threads used
 ('MPI_THREAD_SINGLE'), '1' for threads used in loops ('MPI_THREAD_FUNNELED'),
 and '2' for benchmarks that need 'MPI_THREAD_MULTIPLE'.  Most regular
 benchmarks will use '0', but OpenMP benchmarks will need to set
 'thread_level' to '1'.
  M*/
#define MPUTIL_INIT(_tl) MPUTIL_Init(_tl)
/*M
  MPUTIL_BEGIN - Begin a block of code to be timed

Notes:
 In the parallel case, the block is executed for a number of processes
 ranging from 1 to the size of 'MPI_COMM_WORLD'.  The block is terminated
 by 'MPUTIL_END'
  M*/
#define MPUTIL_BEGIN while (MPUTIL_IterBegin()) { if (MPUTIL_InIter()){
/*M
  MPUTIL_END - End a block of code to be timed

 Notes:
 This macro ends the block of code to be timed that was started with
 'MPUTIL_BEGIN'.
  M*/
#define MPUTIL_END }MPUTIL_IterEnd();}
/*M
  MPUTIL_SYNC - Perform a barrier for all active cores

 Notes:
 This provides a barrier for the active cores.  'MPUTIL_SYNC' is only valid
 within a block that begins with 'MPUTIL_BEGIN' and ends with 'MPUTIL_END'.
  M*/
#define MPUTIL_SYNC MPUTIL_Sync()
/*M
 MPUTIL_FINALIZE - Finalizes the MPUTIL multicore benchmarking tools

 Notes:
 Matches 'MPUTIL_INIT'.
 M*/
#define MPUTIL_FINALIZE MPUTIL_Finalize()
/*M
  MPUTIL_ABORT - Abort a program

  Notes:
  This uses 'exit' in the non-parallel case and 'MPI_Abort' in the parallel
  case.
  M*/
#define MPUTIL_ABORT(_r) MPUTIL_Abort(_r)
#define MPUTIL_SEQBEGIN
#define MPUTIL_SEQEND
#define MPUTIL_GETDMIN(_t,_tmin) (_tmin) = MPUTIL_GetDmin(_t)
#define MPUTIL_GETDMAX(_t,_tmax) (_tmax) = MPUTIL_GetDmax(_t)
#define MPUTIL_GETDSTAT(_t,_tv)
/*M MPUTIL_LABEL - Label a row of output

Input Parmeters:
. args - printf style arguments (varargs).

Notes:
  M*/
#define MPUTIL_LABEL(...) MPUTIL_OutLabel(__VA_ARGS__)
/*M MPUTIL_OUTAPP - Append output to a row of output

Input Parmeters:
. args - printf style arguments (varargs).

Notes:
 This macro should be called after 'MPUTIL_LABEL' to add data to the row
 of output data. A newline (''\n'') should be used for the last value
 intended for the row.

 When used to benchmark on muiltiple cores, the macros 'MPUTIL_LABEL'
 and 'MPUTIL_OUTALL' together should be used to output timing data.
 When used on a single core (the non-parallel option), each acts like
 'printf'.  In the parallel case, the output has one column containing the
 data specified by 'MPUTIL_OUTAPP' for each number of active cores; the
 leftmost column is the text provided by 'MPUTIL_LABEL'.

 The output from both 'MPUTIL_LABEL' and 'MPUTIL_OUTAPP' is generated
 when the timing block is completed with 'MPUTIL_END'.
  M*/
#define MPUTIL_OUTAPP(...) MPUTIL_OutApp(__VA_ARGS__)
/*M MPUTIL_MASTER_BEGIN - Begin a block of code to be executed by only the
  master process

 Notes:
 This macro defines a block of code to be executed by only a single process,
 even in the parallel case.  For example, it can be used to output a
 progress indicator, as in\:
.vb
 MPUTIL_BEGIN
 ...
 MPUTIL_MASTER_BEGIN
 printf("Running test with %d processes\n", MPUTIL_CURSIZE)
 MPUTIL_MASTER_END
 ...
 MPUTIL_END
.ve
  M*/
#define MPUTIL_MASTER_BEGIN do { if (MPUTIL_IsMaster()) {
/*M MPUTIL_MASTER_END - End a block of code to be executed by one process

 Notes:
 Matches 'MPUTIL_MASTER_BEGIN'.
  M*/
#define MPUTIL_MASTER_END   }} while (0)
/*M MPUTIL_NUMPROCESSES - Provide the number of available processes

Notes:
 This returns the size of 'MPI_COMM_WORLD' in the parallel case. In the
 non-parallel case, the result is always '1'.
  M*/
#define MPUTIL_NUMPROCESSES MPUTIL_NumProcs()
/*M MPUTIL_MYRANK - Provide the rank of the process

Notes:
 This returns the rank in 'MPI_COMM_WORLD' of the calling process in the
 parallel case.  In the non-parallel case, the result is always '0'.
  M*/
#define MPUTIL_MYRANK       MPUTIL_MyRank()
/*M MPUTIL_CURSIZE - Number of active processes

Notes:
 This returns the number of processes that will execute the block of
 code between 'MPUTIL_BEGIN' and 'MPUTIL_END'.
  M*/
#define MPUTIL_CURSIZE      MPUTIL_Cursize()

/*M MPUTIL_SETDEBUG - Set the MPUTIL internal debugging flag

Input Parameter:
. flag = value for debug flag.  Zero turns off debugging; non-zero enables
debugging

Return value:
 The previous value of the debug flag is returned.

Notes:
 In the non-parallel case, this routine has no effect and always returns
 zero,
  M*/
#define MPUTIL_SETDEBUG(flag) MPUTIL_SetDebug(flag)
#else

#define MPUTIL_INIT(_tl)
#define MPUTIL_BEGIN
#define MPUTIL_END
#define MPUTIL_SYNC
#define MPUTIL_FINALIZE
#define MPUTIL_ABORT(_r) exit(_r)
#define MPUTIL_SEQBEGIN do {
#define MPUTIL_SEQEND } while (0)
#define MPUTIL_GETDMIN(_t,_tmin) _tmin = (_t)
#define MPUTIL_GETDMAX(_t,_tmax) _tmax = (_t)
#define MPUTIL_GETDSTAT(_t,_tv)  _tv[0] = _tv[1] = _tv[2] = _tv[3] = (_t)
#define MPUTIL_LABEL(...) printf(__VA_ARGS__)
#define MPUTIL_OUTAPP(...) printf(__VA_ARGS__)
#define MPUTIL_MASTER_BEGIN do {
#define MPUTIL_MASTER_END } while(0)
#define MPUTIL_NUMPROCESSES 1
#define MPUTIL_MYRANK       0
#define MPUTIL_CURSIZE      1
#define MPUTIL_SETDEBUG(flag) 0

#endif /* MPUTIL_PARALLEL */
#endif /* HAVE_MPUTIL_H */
