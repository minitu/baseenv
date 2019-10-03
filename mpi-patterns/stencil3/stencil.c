/*
 * Test different options for stencil (halo) communication
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "stencil.h"

#ifdef HAVE_TOPOINIT
#include "topoinfo.h"
#endif

static int outputSoln = 0;
#ifdef DEBUG_HALO
int debugFirst=1, debugLast=3;
#endif

#define MAX_TEST_TYPES 20
#ifndef MAX_COMM_TYPES
#define MAX_COMM_TYPES 3
#endif

/* Set a default for the number of processes per node */
/* FIXME: This needs to be determined by the topology packages and/or
   provided as an option by the user */
static int ppn = 16;

static FILE *fp;

/*
 * All tests use a similar form for the computation and for the
 * communication steps.
 */
/*D
  stencil - Measure the performance of different halo-exchange approaches in MPI for a regular stencil code

  Input Parameters:
+ nrange - Range of problem sizes.  See below
. energy - Amount of energy to input into the simulation
. niters - Number of iterations to perform
. px -     Number of processes in the x-direction
. py -     Number of processes in the y-direction
- filename - Name of output file (stdout if not provided)

  Output Results:
 The output is printed as a tab-separated table, indicating the
 time for the different operations and halo-exchange choices, as well as an
 achieved computation rate.  The default output is to stdout, but can be
 changed with a commandline argument.

Notes:
  The 'nrange' parameter provides the size of a size of the mesh, and may
  be either a single value, a comma-separated list, or a range in the form
  'start:end:increment'.

  The program should be compiled with optimization and with vectorization
  enabled if the computation rates are to be used.  Note that the computation
  code does not compute the heat per iteration unless the code is compiled
  with '-DCOMPUTE_HEAT_EACH_ITERATION=1'.  Including this code causes the
  Cray compiler (at least) to fail to vectorize the application of the stencil
  in the same loop, causing a significant loss in performance.
  D*/
int main(int argc, char *argv[])
{
    stencilInfo si;
    stencilTrialTime st[MAX_TEST_TYPES];
    procInfo    pi;
    probDesc    pd;
    double finalHeat[MAX_TEST_TYPES];
    int provided, i, j, wsize, wrank;
    int commidx, sizeidx, maxcommidx=MAX_COMM_TYPES;
    int nt, ntimes = 10;
    int topotype;
    long nwork; /* Number of floating point operations in each (complete)
		   stencil operation */

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

#if 0
    MPIX_Nodecart_cvar_set("ppn", 4);
#endif

    fp = stdout;  /* A default for the output */
    /* Get command-line arguments about problem size etc. */
    getArgs(argc, argv, &pd, &pi, &si, &fp, &maxcommidx);
    /* All tests must be run at least 2 times */
    if (ntimes < 2) ntimes = 2;
#ifdef HAVE_TOPOINIT
    initTopology(0);
#endif

    /* Header to identify the run */
    if (wrank == 0) {
	int  rlen;
	char name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(name, &rlen);
	fprintf(fp, "stencil (3D) for baseenv version %s (%s)\n",
		PACKAGE_VERSION, BASEENV_GITHASH);
	fprintf(fp, "Processor name %s\n", name);
	fprintf(fp, "WORLD size %d; processor mesh is %d x %d x %d\n",
		wsize, pi.px, pi.py, pi.pz);
	fflush(fp);
    }

    for (commidx=0; commidx<maxcommidx; commidx++) {
	for (sizeidx=0; sizeidx < si.nSizes; sizeidx++) {
	    si.n = si.nSize[sizeidx];

	    /* 6 adds, 2 multiplies, 1 add for total heat */
	    nwork = (long)si.n * (long)si.n * (long)si.n * 8;
#if COMPUTE_HEAT_EACH_ITERATION
	    /* Add the cost of summing the heat within the stencil loop */
	    nwork += (long)si.n * (long)si.n + (long)si.n;
#endif
	    nwork *= si.niters;

	    /* Compute the decomposition of the domain */
	    switch (commidx) {
	    case 0:
		getDecomp(MPI_COMM_WORLD, &pi, &si);
		break;
	    case 1:
		getDecompCart(MPI_COMM_WORLD, &pi, &si);
		break;
	    case 2:
		getDecompNodecart(MPI_COMM_WORLD, &pi, &si);
		if (wrank == 0 && sizeidx == 0) {
		    fprintf(fp, "Nodecart process decomp:");
		    MPIX_PrintNodeDecomp(fp, pi.comm);
		}
		break;

		    /* Add other decompositions here:
		       conseq on SMP
		       round robin across nodes
		       round robin across chips (if more than one chip per node)
		    */
#ifdef HAVE_TOPOINIT
		    /* If the topo routines available, use them */
#endif
		default:
		    myAbort(MPI_COMM_WORLD, 1, "Internal error in comm choice");
	    }
	    /* Check that the communicator is different that COMM_WORLD */
	    if (commidx != 0 && pi.skipCongruent) {
		int flag;
		MPI_Comm_compare(MPI_COMM_WORLD, pi.comm, &flag);
		if (flag == MPI_CONGRUENT) {
		    MPI_Topo_test(pi.comm, &topotype);
		    if (!pi.withCart ||
			(pi.withCart &&
		       (topotype != MPI_CART && topotype != MPI_DIST_GRAPH))) {
			/* Free the created communicator */
			MPI_Comm_free(&pi.comm);
			continue;
		    }
		}
	    }

#if 0
	    /* DEBUGING */
	    if (sizeidx == 0) {
		if (wrank == 0) {
		    fprintf(stdout, "---- comm idx %d\n", commidx);
		}
		printDecomp(stdout, pi.comm, &pi, &si);
	    }
#endif

	    /* Given the stencil, initialize the problem description */
	    initProbDesc(&pd, &si);

	    /* Initialize the timing structures */
	    for (i=0; i<MAX_TEST_TYPES; i++) initST(&st[i], ntimes);

	    allocMesh(&si);
	    /* Run the tests multiple times (ntimes) */
	    for (nt=0; nt<ntimes; nt++) {
		/* Run the various stencil tests */
		int testidx = 0;
		initMesh(&pd, &si);
		st[testidx].label = "BP-U";
		finalHeat[testidx] = stencil_bnb(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);

		testidx++;
		initMesh(&pd, &si);
		st[testidx].label = "NP-U";
		finalHeat[testidx] = stencil_nb(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);
#if 0
		testidx++;
		initMesh(&pd, &si);
		st[testidx].label = "NP-UP";
		finalHeat[testidx] = stencil_pnb(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);

		testidx++;
		initMesh(&pd, &si);
		st[testidx].label = "NP-D";
		finalHeat[testidx] = stencil_nb_ddt(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);

		testidx++;
		initMesh(&pd, &si);
		st[testidx].label = "NPOD";
		finalHeat[testidx] = stencil_ddt_ov(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);

#ifndef NO_RMA
		if (pi.useRMAdtype) {
		    testidx++;
		    stencil_ddt_rma_init(&pi, &si, &st[testidx].st[nt]);
		    initMesh(&pd, &si);
		    st[testidx].label = "NR-D";
		    finalHeat[testidx] = stencil_ddt_rma(&pd, &pi, &si, &st[testidx].st[nt]);
		    if (outputSoln && nt == 0)
			printSoln(&si, si.aold, st[testidx].label);
		    stencil_ddt_rma_free(&pi, &si, &st[testidx].st[nt]);
		}
		if (pi.useRMA) {
		    testidx++;
		    stencil_rma_init(&pi, &si, &st[testidx].st[nt]);
		    initMesh(&pd, &si);
		    st[testidx].label = "NR-U";
		    finalHeat[testidx] = stencil_rma(&pd, &pi, &si, &st[testidx].st[nt]);
		    if (outputSoln && nt == 0)
			printSoln(&si, si.aold, st[testidx].label);
		    stencil_ddt_rma_free(&pi, &si, &st[testidx].st[nt]);
		}
#endif
		MPI_Topo_test(pi.comm, &topotype);
		if (topotype == MPI_CART || topotype == MPI_DIST_GRAPH) {
		    testidx++;
		    initMesh(&pd, &si);
		    st[testidx].label = "NC-U";
		    finalHeat[testidx] = stencil_neighcolls(&pd, &pi, &si, &st[testidx].st[nt]);
		    if (outputSoln && nt == 0)
			printSoln(&si, si.aold, st[testidx].label);

		    /* Add other neighbor collective use here */
		    testidx++;
		    initMesh(&pd, &si);
		    st[testidx].label = "NCOU";
		    finalHeat[testidx] = stencil_neighcolls_ov(&pd, &pi, &si, &st[testidx].st[nt]);
		    if (outputSoln && nt == 0)
			printSoln(&si, si.aold, st[testidx].label);

		}

/* IBM BlueGene/Q systems may only have MPI-2 */
#if MPI_VERSION >= 3
		testidx++;
		stencil_shmem_nb_init(&pi, &si, &st[testidx].st[nt]);
		st[testidx].label = "NS-U";
		initMesh(&pd, &si);
		finalHeat[testidx] = stencil_shmem_nb(&pd, &pi, &si, &st[testidx].st[nt]);
		if (outputSoln && nt == 0)
		    printSoln(&si, si.aold, st[testidx].label);
		stencil_shmem_nb_free(&pi, &si, &st[testidx].st[nt]);
#endif
#endif

		/* Mark the first case for which ran no tests */
		st[testidx+1].st[0].commInit = -1.0;
	    }
	    freeMesh(&si);

	    /* Report the results */
	    if (sizeidx == 0) {
		if (wrank == 0) {
		    char commname[MPI_MAX_OBJECT_NAME];
		    int  nlen;
		    commname[0] = 0;
		    MPI_Comm_get_name(pi.comm, commname, &nlen);
		    if (nlen > 0 && commname[0])
			fprintf(fp, "Communicator\t%s\n", commname);
		}
#ifdef HAVE_TOPOINIT
		getNodeCommunicationStats(&pi, fp);
#endif
		if (wrank == 0) {
		    fprintf(fp, "n\tmethod\theat    \tcommInit\tcommStart\tcommCmpl\tcommPack\tcommUnpk\tInterior\tBndy    \tUpdate  \tTotal T \tTotal T(avg)\tGFLOP/s  \tGFLOP/s (avg)\n");
		}
	    }
/* This block from stencil.h indicates the components, in order, of the
   time values in the arrays in the timing struture */
#if 0
	    double commInit,   /* Time to initialize communication structures */
		commStart,     /* Time to initiate communication */
		commComplete,  /* Time to complete communication */
		commPack,      /* Time to pack data */
		commUnpack,    /* Time to unpack data */
		compInterior,  /* Time to compute interior elements */
		compBndy,      /* Time to compute boundary elements */
		compUpdate;    /* Time to compute all elements (disjoint from
				  Bndy/Interior) */
#endif
	    for (i=0; i<MAX_TEST_TYPES; i++) {
		double tmax[N_TIMEVALS], tmin[N_TIMEVALS], tavg[N_TIMEVALS];

		if (st[i].st[0].commInit < 0) break;
		/* First: for each of the trials (nt), compute the min, max, and
		   mean time on each process */
		for (j=0; j<N_TIMEVALS; j++) {
		    tmax[j] = 0.0; tavg[j] = 0.0; tmin[j] = 1.0e12;
		}
		/* For each test, get the min, max, and average of each
		   measured time.  We *always* skip the first test */
		for (nt=1; nt<ntimes; nt++) {
		    double *t = (double *)&st[i].st[nt];
		    for (j=0; j<N_TIMEVALS; j++) {
			if (t[j] > tmax[j]) tmax[j] = t[j];
			if (t[j] < tmin[j]) tmin[j] = t[j];
			tavg[j] += t[j];
		    }
		}
		for (j=0; j<N_TIMEVALS; j++) {
		    tavg[j] /= (ntimes-1);
		}

		/* Get the min,max,average times across all processes */
		MPI_Allreduce(MPI_IN_PLACE, tmax, N_TIMEVALS, MPI_DOUBLE,
			      MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, tmin, N_TIMEVALS, MPI_DOUBLE,
			      MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, tavg, N_TIMEVALS, MPI_DOUBLE,
			      MPI_SUM, MPI_COMM_WORLD);
		for (j=0; j<N_TIMEVALS; j++) tavg[j] /= wsize;

		if (wrank == 0) {
		    fprintf(fp, "%d\t%.8s\t%.3e\t", si.n, st[i].label,
			    finalHeat[i]);
		    /* By default, output the min times (these are the
		       most reproducible) */
		    /* FIXME: option to add average */
		    for (j=0; j<N_TIMEVALS; j++) {
			if (tmax[j] > 0) {
			    fprintf(fp, "%.2e\t", tmin[j]);
			}
			else fprintf(fp, "        \t");
		    }
		    fprintf(fp, "%.2e\t", tavg[TOTAL_TIME]);
		    /* Compute a rate based on the best time measured
		       and the total work */
		    fprintf(fp, "%.2e\t%.2e\n",
			    1.0e-9 * (double)nwork / tmin[TOTAL_TIME],
			    1.0e-9 * (double)nwork / tavg[TOTAL_TIME]);
		    fflush(fp);
		}
	    }
	}
    }

#ifdef HAVE_TOPOINIT
    endTopology();
#endif

    MPI_Finalize();

    return 0;
}
