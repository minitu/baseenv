
#include "baseenv.h"
#include "../nodecart/nodecart.h"

/* use debugFirst, debugLast to control when output is generated.  A
   typical test is
   if (dbgcnt >= debugFirst && dbgcnt <=debugLast) {
       dbgcnt++; printSoln(...);
   }
*/
/* #define DEBUG_HALO 1*/
#ifdef DEBUG_HALO
extern int debugFirst, debugLast;
#endif

#ifdef USE_ROW_MAJOR
#define ind(i,j) (j)*(bx+2)+(i)
#define NS_CONTIG 1
#elif defined(USE_COLUMN_MAJOR)
#define ind(i,j) (i)*(by+2)+(j)
#define EW_CONTIG 1
#error Not implemented
#else
#error Must pick row or column major
#endif

/* Define this to compute heat each iteration rather than once at the
   end of the iteration.  On some systems/compilers, computing heat
   along with the stencil relaxation inhibits vectorization (for example,
   the Cray compiler on Blue Waters (Cray XE6) has this behavior */
#ifndef COMPUTE_HEAT_EACH_ITERATION
#define COMPUTE_HEAT_EACH_ITERATION 0
#endif

/* Define if a type created with MPI_Type_contiguous should be used in the
   derived datatype examples even though the data is contiguous.  Only
   partially implemented */
#define USE_CONTIG_TYPE 1

typedef struct {
    int      px, py;        /* Processor in x and y dimensions */
    int      north, south, east, west; /* rank of Neighbors in comm */
    MPI_Comm comm;          /* Communicator for communication */
    int      useRMA;        /* Whether RMA should be used */
    int      useRMAdtype;   /* Whether RMA should be used with derived
			       datatypes */
    int      skipCongruent; /* Set to true to skip comm congruent to WORLD */
    int      withCart;      /* If true, do not skip Cart_create comm */
} procInfo;

#define NSOURCES 3
typedef struct {
    int isForcingTerm;           /* Set to 1 if sources applies for t>0 */
    int locnsources;             /* number of sources in my area */
    int locsources[NSOURCES][2]; /* sources local to my rank */
    double energy;               /* Energy to add at each source */
} probDesc;

/* N_TIMEVALS is the number of time values in stencilTime */
#define N_TIMEVALS 9
#define TOTAL_TIME 8
typedef struct {
    double commInit,   /* Time to initialize communication structures */
	commStart,     /* Time to initiate communication */
	commComplete,  /* Time to complete communication */
 	commPack,      /* Time to pack data */
        commUnpack,    /* Time to unpack data */
	compInterior,  /* Time to compute interior elements */
	compBndy,      /* Time to compute boundary elements */
	compUpdate,    /* Time to compute all elements (disjoint from
			  Bndy/Interior) */
	total;         /* Total time in stencil iteration */
} stencilTime;

/* Use this struct to record a collection of trials of the same experiment
   The label uses a condensed form, to simplify the tabular output:
   [B/N][P/C/R/S][-/O][U/D][vv]

   B/N - Blocking or nonblocking
   P/C/R/S - Point-2-point, collective, RMA, or shared memory
   -/O - O if structured for communication overlap
   U/D - D if datatypes used instead of user pack, U for user pack/unpack
   v   - Indicates a variety.  For example,
        NPODR - Nonblocking, point-2-point, overlapping, with datatypes, irsend
	NRODF - Nonblocking, RMA, overlapping, with datatypes, fence
	NRODFA - Like NRODF, but with fence asserts
 */
typedef struct {
    int navail,        /* Number of stencilTime available (size of array) */
        ncur;          /* Index of last in use (-1 for none) */
    stencilTime *st;
    const char *label; /* Description of the tested approach */
} stencilTrialTime;

typedef struct {
    double *aold, *anew;  /* Local mesh, may include 1-wide halo zones */
    int n,             /* Mesh is a square n x n (global size) */
	bx, by,        /* Block sizes (local mesh size) */
	offx, offy,    /* (x,y) index for upper left corner of mesh */
	rx, ry;        /* (x,y) index for process in process mesh */
    /* Some problem information */
    int niters;        /* Number of iterations */
    /* Stencil sizes to run */
    int nSizes, *nSize;
} stencilInfo;

/* The application of each approach for applying a stencil operation for the
 heat equation has two parts:
    stencil_xxx_init - Initialize the memory for the solution (aold,anew)
    stencil_xxx      - Apply the stencil for niters iterations
    stencil_xxx_free - Release any resources allocated in the init step
 Many approaches don't require any init or free steps, and provide empty
 routines.
*/


double stencil_bnb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_nb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_pnb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_nb_ddt(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_ddt_ov(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_ddt_rma(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_ddt_rma_init(procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_ddt_rma_free(procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_shmem_nb(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_rma(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_rma_init(procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_rma_free(procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_shmem_nb_init(procInfo *pi, stencilInfo *si, stencilTime *st);
void stencil_shmem_nb_free(procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_neighcolls(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);
double stencil_neighcolls_ov(probDesc *pd, procInfo *pi, stencilInfo *si, stencilTime *st);

void getArgs(int argc, char *argv[],
	     probDesc *pd, procInfo *pi, stencilInfo *si, FILE **fp, int *);
void myAbort(MPI_Comm, int, const char *);
void getDecomp(MPI_Comm comm, procInfo *pi, stencilInfo *si);
void getDecompCart(MPI_Comm comm, procInfo *pi, stencilInfo *si);
void getDecompNodecart(MPI_Comm comm, procInfo *pi, stencilInfo *si, int wsock);
void initST(stencilTrialTime *, int);
void initProbDesc(probDesc *pd, stencilInfo *si);
void initMesh(probDesc *pd, stencilInfo *si);
void initMeshBase(probDesc *pd, int bx, int by, double *restrict aold,
		  double *restrict anew);
void printarr(FILE *fp, double *a, int n, int bx);
void printSoln(stencilInfo *si, double *v, const char *name);
void allocMesh(stencilInfo *si);
void freeMesh(stencilInfo *si);

/* Routines to take advantage of the topoinfo routines */
void initTopology(int v);
void getNodeCommunicationStats(procInfo *pi, FILE *fp);
void endTopology(void);
