/* -*- Mode: C; c-basic-offset:4 ; -*- */

/* This requires that the include path include the baseenv root, including
 the nodecart headers */
#include "baseenv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "nodecart.h"
#include "util/benvutil.h"

/* This code tests three different communicators for a (1,2,3)-D halo
 exchange: MPI_COMM_WORLD, the communicator from MPI_Cart_create,
 and the communicator from MPIX_Nodecart_create.
*/

#define MAX_REQUESTS 2*6
#define MAX_NTEST 20
double timeHalo(MPI_Comm comm, int nreps,
		int nr, const int ranks[], int msgsize);

static int verbose = 0;
static int outputProcessmapping = 0;
static FILE *fp = 0;

int totalpartners(MPI_Comm comm, int nr, const int ranks[], int *tot);
void getMinMaxAvg(int len, const double v[],
		  double *vmin, double *vmax, double *vavg);

#define MAX_DIMS 3
int main(int argc, char *argv[])
{
    int      wrank, wsize;
    MPI_Comm cartcomm, ncartcomm, ncartcommh;
    int      dims[MAX_DIMS], periods[MAX_DIMS], debug=0, i, k,
	ntest=MAX_NTEST, sz[2];
    int      *meshdims, ndims, ndim;   /* For the mesh dimensions (e.g., =2) */
    int      cartranks[2*MAX_DIMS], ncartranks[2*MAX_DIMS],
	ncartranksh[2*MAX_DIMS];
    MPI_Comm nodecomm, leadercomm, socketcomm;
    int      nsize, nnodes, noderank, n, nsizes;
    int      cartdim, *msgsizes = 0, nreps=32;
    int      nc, *commnodes, ncommnodes, ssize;
    MPI_Comm subsetcomm;
    int      totcr, totncr;
    int      worldonly = 0;      /* Make it possible to run only COMM_WORLD */
    double   tcart[MAX_NTEST], tncart[MAX_NTEST], tncarth[MAX_NTEST];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Process command line arguments */
    msgsizes  = BENV_GetSizes(argc, argv, "sizes", &nsizes);
    commnodes = BENV_GetSizes(argc, argv, "nodes", &ncommnodes);
    meshdims  = BENV_GetSizes(argc, argv, "dims", &ndims);
    fp        = stdout;
    for (i=1; i<argc; i++) {
	char *v = argv[i];
	/* Skip arguments already handled */
	if (v == 0) continue;
	if (*v == '-') v++;
	if (*v == '-') v++;
	if (strcmp(v, "debug") == 0) debug = 1;
	else if (strcmp(v, "worldonly") == 0) worldonly = 1;
	else if (strcmp(v, "v") == 0) verbose = 1;
	else if (strcmp(v, "nodemapping") == 0) outputProcessmapping = 1;
	else if (strcmp(v, "o") == 0) {
	    /* Only process 0 opens the file */
	    if (wrank == 0) {
		fp = fopen(argv[++i], "w");
		if (!fp) {
		    fprintf(stderr, "Unable to open file %s\n", argv[i]);
		    MPI_Abort(MPI_COMM_WORLD, 1);
		}
	    }
	    else i++;  /* Skip filename argument */
	}
	else if (strcmp(v, "ntest") == 0) {
	    ntest = atoi(argv[++i]);
	    if (ntest > MAX_NTEST && wrank == 0) {
		fprintf(stderr,
			"Maximum value of ntest is %d, given value was %d\n",
			MAX_NTEST, ntest);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
	else {
	    if (wrank == 0) {
		fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
		fprintf(stderr,
			"Usage: ncartperf --debug --v     Debugging options\n\t--nodemapping     Output the node mapping\n\t--worldonly      Use with rank mapping tools\n\t-o filename --ntest n --sizes a:b:/*c --nodes a:b:/*c --dims a:b\n");
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD,1);
	    }
	}
    }

    /* Handle defaults */
    /* Run all tests at least twice */
    if (ntest < 2) ntest = 2;
    /* Get a default range of sizes */
    if (nsizes <= 0 && !msgsizes) {
	msgsizes = BENV_GetSizesMult(128, 66000, 2.0, &nsizes);
    }
    if (ndims <= 0) {
	meshdims = BENV_GetSizesArith(2, 3, 1, &ndims);
    }
    else {
	/* Validate dimensions */
	for (i=0; i<ndims; i++) {
	    if (meshdims[i] < 0 || meshdims[i] > MAX_DIMS) {
		if (wrank == 0) {
		    fprintf(stderr,
			    "Invalid values for --dims: %d not in [0,%d]\n",
			    meshdims[i], MAX_DIMS);
		    fflush(stderr);
		    MPI_Abort(MPI_COMM_WORLD, 1);
		}
	    }
	}
    }

    /* Set the default for commnodes after we know the number of nodes */

    /* Temp: For testing on laptop, use debug option to set the node size */
    if (debug) {
	if (wsize > 32 && (wsize % 16) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 16);
	}
	else if (wsize >= 12 && (wsize % 6) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 6);
	    MPIX_Nodecart_cvar_set("debug", 2);
	}
	else if (wsize >=8 && (wsize % 4) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 4);
	}
	else if ((wsize & 0x1) == 0) {
	    MPIX_Nodecart_cvar_set("ppn", 2);
	}
    }

    /* Get the node size */
    MPIX_Nodecomm_create(MPI_COMM_WORLD, &nodecomm, &leadercomm,
			 &noderank, &nnodes);
    /* Get the socket info (if any).  This is used to get more information
       about the mapping to neigbors */
    MPIX_Socketcomm_create(nodecomm, &socketcomm);

    if (ncommnodes <= 0 || !commnodes) {
	commnodes = BENV_GetSizesMult(nnodes, nnodes, 2.0, &ncommnodes);
    }

    /* Sanity check that all nodes are the same size */
    MPI_Comm_size(nodecomm, &nsize);
    sz[0] = nsize;
    sz[1] = -nsize;
    MPI_Allreduce(MPI_IN_PLACE, sz, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (sz[0] != -sz[1]) {
	if (wrank == 0) {
	    fprintf(stderr, "Nodes must have the same number of processes!\n");
	    fprintf(stderr, "Found range between %d-%d\n", -sz[1], sz[0]);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
    }

    /* Output header with version and run info */
    if (wrank == 0) {
 	int  rlen;
	char name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(name, &rlen);
	fprintf(fp, "ncartperf for baseenv version %s (%s)\n",
		PACKAGE_VERSION, BASEENV_GITHASH);
	fprintf(fp, "SMP: nodes = %d, nodesize = %d\n", nnodes, nsize);
	fprintf(fp, "Processor name %s\n", name);
	fflush(fp);
    }

    /* Print the processes in the same node for each leader */
    if (leadercomm != MPI_COMM_NULL && outputProcessmapping) {
	int j, nodesize, *nranks, *nwranks, nodenum;
	MPI_Group gworld, ngroup;

	MPI_Comm_size(nodecomm, &nodesize);
	MPI_Comm_group(nodecomm, &ngroup);
	nranks = (int *)malloc(2*nsize*sizeof(int));
	nwranks = nranks + nodesize;
	for (j=0; j<nodesize; j++) nranks[j] = j;

	MPI_Comm_group(MPI_COMM_WORLD, &gworld);

	MPI_Group_translate_ranks(ngroup, nodesize, nranks, gworld, nwranks);

	MPI_Comm_rank(leadercomm, &nodenum);

	/* Simple sequentialization nodes */
	/* Write to stdout since all processes are performing the output */
	for (i=0; i<nnodes; i++) {
	    if (i == nodenum) {
		fprintf(stdout, "%d: process on node %d:", wrank, nodenum);
		for (j=0; j<nodesize; j++) fprintf(stdout,"%d ", nwranks[j]);
		fprintf(stdout, "\n");fflush(stdout);
	    }
	}
	free(nranks);
	MPI_Group_free(&ngroup);
    }

    for (nc=0; nc<ncommnodes; nc++) {
	if (commnodes[nc] > nnodes) {
	    if (wrank == 0) {
		fprintf(stderr, "Number of nodes given (%d) exceeds number available (%d)\n",
			commnodes[nc], nnodes);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	    break;
	}
	MPIX_NodeSubset(MPI_COMM_WORLD, leadercomm, nodecomm,
			commnodes[nc], &subsetcomm);

	if (subsetcomm == MPI_COMM_NULL) continue;
	for (ndim = 0; ndim<ndims; ndim++) {
	    cartdim = meshdims[ndim];

	    /* MPI Cartesian routine */
	    for (i=0; i<cartdim; i++) {
		dims[i]    = 0;
		periods[i] = 0;
	    }
	    MPI_Comm_size(subsetcomm, &ssize);
	    if (wrank == 0) {
		fprintf(fp,"Subset nodes,size\t%d\t%d\n", commnodes[nc], ssize);
		fflush(fp);
	    }

	    MPI_Dims_create(ssize, cartdim, dims);
	    if (wrank == 0) {
		fprintf(fp, "Mesh is ");
		for (i=0; i<cartdim; i++)
		    fprintf(fp, "%d%c", dims[i], (i!=cartdim-1)?'x':'\n');
		fflush(fp);
	    }

	    MPI_Cart_create(subsetcomm, cartdim, dims, periods, 1, &cartcomm);
	    MPI_Comm_set_name(cartcomm, "cartcomm");
	    for (i=0; i<cartdim; i++) {
		MPI_Cart_shift(cartcomm, i, 1,
			       &cartranks[2*i], &cartranks[2*i+1]);
	    }
	    totalpartners(cartcomm, 2*cartdim, cartranks, &totcr);

	    /* Check whether anything has changed in the rank order */
	    if (wrank == 0) {
		int cresult;
		MPI_Comm_compare(subsetcomm, cartcomm, &cresult);
		if (cresult == MPI_CONGRUENT) {
		    fprintf(fp, "subsetcomm and cartcomm have congruent groups\n");
		}
		else {
		    fprintf(fp, "subsetcomm and cartcom have different groups\n");
		}
		fflush(fp);
	    }
	    MPI_Barrier(subsetcomm);
	    if (socketcomm != MPI_COMM_NULL) {
		MPIX_PrintNodeCommCounts_X(fp, cartcomm, 2*cartdim, cartranks,
					   nodecomm, socketcomm);
	    }
	    else {
		MPIX_PrintNodeCommCounts(fp, cartcomm, 2*cartdim, cartranks,
					 nodecomm);
	    }


	    /* MPIX Nodecart routine */
	    MPIX_Nodecart_create(subsetcomm, cartdim, dims, periods, 1,
				 &ncartcomm);
	    MPI_Comm_set_name(ncartcomm, "ncartcomm");
	    for (i=0; i<cartdim; i++) {
		MPIX_Nodecart_shift(ncartcomm, i, 1,
				    &ncartranks[2*i], &ncartranks[2*i+1]);
	    }
	    totalpartners(ncartcomm, 2*cartdim, ncartranks, &totncr);

	    MPI_Barrier(subsetcomm);
	    if (wrank == 0) {
		fprintf(fp, "nodecart process decomp:");
		MPIX_PrintNodeDecomp(fp, ncartcomm);
	    }
	    if (socketcomm != MPI_COMM_NULL) {
		MPIX_PrintNodeCommCounts_X(fp, ncartcomm, 2*cartdim, ncartranks,
					   nodecomm, socketcomm);
	    }
	    else {
		MPIX_PrintNodeCommCounts(fp, ncartcomm, 2*cartdim, ncartranks,
					 nodecomm);
	    }

	    /* Sanity check */
	    if (wrank == 0) {
		if (totcr != totncr) {
		    fprintf(fp, "Unexpected total partners: %d vs %d\n",
			    totcr, totncr);
		    fflush(fp);
		}
	    }

	    /* MPIX Nodecart routine, for multisocket nodes */
	    MPIX_SetUseSocket(1);
	    MPIX_Nodecart_create(subsetcomm, cartdim, dims, periods, 1,
				 &ncartcommh);
	    MPIX_SetUseSocket(0);
	    MPI_Comm_set_name(ncartcommh, "ncartcomm-h");
	    for (i=0; i<cartdim; i++) {
		MPIX_Nodecart_shift(ncartcommh, i, 1,
				    &ncartranksh[2*i], &ncartranksh[2*i+1]);
	    }
	    totalpartners(ncartcommh, 2*cartdim, ncartranksh, &totncr);

	    MPI_Barrier(subsetcomm);
	    if (wrank == 0) {
		fprintf(fp, "nodecart-h process decomp:");
		MPIX_PrintNodeDecomp(fp, ncartcommh);
	    }
	    if (socketcomm != MPI_COMM_NULL) {
		MPIX_PrintNodeCommCounts_X(fp, ncartcommh,
					   2*cartdim, ncartranksh,
					   nodecomm, socketcomm);
	    }
	    else {
		MPIX_PrintNodeCommCounts(fp, ncartcommh, 2*cartdim, ncartranksh,
					 nodecomm);
	    }

	    /* Sanity check */
	    if (wrank == 0) {
		if (totcr != totncr) {
		    fprintf(fp, "Unexpected total partners: %d vs %d\n",
			    totcr, totncr);
		    fflush(fp);
		}
	    }
	    if (wrank == 0) {
		int cresult;
		MPI_Comm_compare(ncartcomm, ncartcommh, &cresult);
		if (cresult == MPI_CONGRUENT) {
		    fprintf(fp, "ncartcomm and ncartcomm-h have congruent groups\n");
		}
		else {
		    fprintf(fp, "ncartcomm and ncartcomm-h have different groups\n");
		}
		fflush(fp);
	    }


	    /* End of communicator constructions */
	    if (wrank == 0) {
		fprintf(fp, "n\tcart min\tcart avg\tcart rate\tnodec min\tnodec avg\tnodec rate\tnodech min\tnodech avg\tnodech rate\n");
		fflush(fp);
	    }
	    for (n=0; n<nsizes; n++) {
		for (k=0; k<ntest; k++) {
		    tcart[k]  = timeHalo(cartcomm,  nreps, 2*cartdim,
					 cartranks,  msgsizes[n]);
		    if (!worldonly) {
			tncart[k] = timeHalo(ncartcomm, nreps, 2*cartdim,
					     ncartranks, msgsizes[n]);
			tncarth[k] = timeHalo(ncartcommh, nreps, 2*cartdim,
					      ncartranksh, msgsizes[n]);
		    }
		}
		/* Compute the min, max, and average */
		if (wrank == 0) {
		    double rcart=0, rncart=0, rnhcart=0;
		    double tcmin, tcmax, tcavg, tncmin, tncmax, tncavg,
			tnchmin, tnchmax, tnchavg;

		    /* We ignore the first run in both cases */
		    getMinMaxAvg(ntest-1, tcart+1, &tcmin, &tcmax, &tcavg);
		    getMinMaxAvg(ntest-1, tncart+1, &tncmin, &tncmax, &tncavg);
		    getMinMaxAvg(ntest-1, tncarth+1,
				 &tnchmin, &tnchmax, &tnchavg);

		    if (tcmin > 0) {
			rcart = 2*cartdim*msgsizes[n]*sizeof(double)/tcmin;
		    }
		    if (tncmin > 0) {
			rncart = 2*cartdim*msgsizes[n]*sizeof(double)/tncmin;
		    }
		    if (tnchmin > 0) {
			rnhcart = 2*cartdim*msgsizes[n]*sizeof(double)/tnchmin;
		    }
		    fprintf(fp,
		  "%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			    msgsizes[n],
			    tcmin, tcavg, rcart, tncmin, tncavg, rncart,
			    tnchmin, tnchavg, rnhcart);
		    fflush(fp);
		}
	    }
	}
	/* Free communicators */
	MPI_Comm_free(&subsetcomm);
	MPI_Comm_free(&cartcomm);
	MPI_Comm_free(&ncartcomm);
	MPI_Comm_free(&ncartcommh);
    }

    /* Npte that the delete function COMM_WORLD will deleted the nodecomm and
       leadercomm */

    free(msgsizes);
    free(commnodes);

    MPI_Finalize();
    return 0;
}

double timeHalo(MPI_Comm comm, int nreps,
		int nr, const int ranks[], int msgsize)
{
    double      t0, tf;
    int         i, reps, wrank;
    MPI_Request rq[MAX_REQUESTS];
    double      *bufptrs[MAX_REQUESTS];

    if (verbose) {
	printf("starting time halo size %d and %d reps...\n",
	       msgsize, nreps);
	fflush(stdout);
    }
    for (i=0; i<2*nr; i++) {
	int    j;
	double *p;

	bufptrs[i] = (double *)malloc(msgsize * sizeof(double));
	if (!bufptrs[i]) {
	    fprintf(stderr, "Null pointer for bufptrs[%d]\n", i);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	/* Initialize buffer */
	p = bufptrs[i];
	for (j=0; j<msgsize; j++) p[j] = 0.0;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (reps=0; reps<nreps; reps++) {
	for (i=0; i<nr; i++) {
	    MPI_Irecv(bufptrs[i], msgsize, MPI_DOUBLE,
		      ranks[i], 0, comm, &rq[i]);
	}
	for (i=0; i<nr; i++) {
	    MPI_Isend(bufptrs[nr+i], msgsize, MPI_DOUBLE,
		      ranks[i], 0, comm, &rq[nr+i]);
	}
	MPI_Waitall(2*nr, rq, MPI_STATUSES_IGNORE);
    }
    MPI_Barrier(comm);
    tf = (MPI_Wtime() - t0) / nreps;
    if (verbose > 0 && wrank == 0) {
	printf("Irecv/Send took %.2e sec\n", tf*nreps);
	fflush(stdout);
    }
    for (i=0; i<2*nr; i++) {
	free(bufptrs[i]);
    }
    return tf;
}

/* Sanity check: count the number of valid partner ranks */
int totalpartners(MPI_Comm comm, int nr, const int ranks[], int *tot)
{
    int i, ltot, err=0, csize;

    MPI_Comm_size(comm, &csize);
    ltot = 0;
    for (i=0; i<nr; i++) {
	if (ranks[i] == MPI_PROC_NULL) continue;
	if (ranks[i] >= 0 && ranks[i] < csize) {
	    ltot++;
	}
	else {
	    /* Invalid rank value */
	    err++;
	}
    }
    MPI_Allreduce(&ltot, tot, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if (err > 0) {
	int crank;
	MPI_Comm_rank(comm, &crank);
	if (crank == 0) {
	    fprintf(stderr, "Invalid ranks found!\n");
	    fflush(stderr);
	}
	return MPI_ERR_OTHER;
    }
    else
	return MPI_SUCCESS;
}

void getMinMaxAvg(int len, const double v[],
		  double *vmin, double *vmax, double *vavg)
{
    double minval, maxval, avgval;
    int    i;

    if (len <= 0) {
	*vmin = 0;
	*vmax = 0;
	*vavg = 0;
	return;
    }
    minval = v[0];
    maxval = v[0];
    avgval = v[0];
    for (i=1; i<len; i++) {
	if (v[i] < minval) minval = v[i];
	if (v[i] > maxval) maxval = v[i];
	avgval += v[i];
    }
    *vavg = avgval / len;
    *vmin = minval;
    *vmax = maxval;
}

