/* -*- Mode: C; c-basic-offset:4 ; -*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/* Needed for restrict definition if using an old C compiler */
#include "baseenv.h"
#include "../util/benvutil.h"
#include "../nodecart/nodecart.h"

/*
 * This program compares the performance of different approaches to
 * moving strided data in MPI.  These include:
 *
 * User pack and unpack, communicating contiguous data
 * Type_vector
 * Type_create_resized
 *
 * The types of strided data are:
 *    (line) Edge of 2-d array
 *    (plane) Edge of 3-d array
 *
 * Communication is by
 *    Send/Recv
 *    Put; lock/unlock
 *    Put; fence
 *
 * Message size given by defining array size
 *
 * Initially, only send/recv used for communication, because the
 * expectation is that RMA datatype performance will be poor.
 */

double timeUser(int, int, const double *, double *, double *, double *,
		    MPI_Comm);
double timeVec(MPI_Datatype, const double *, double *, MPI_Comm);
double timeResize(MPI_Datatype, int, const double *, double *, MPI_Comm);
double timeUser3(int, int, const double *, double *, double *, double *,
		    MPI_Comm);
double timeUser3D(int, int, int, const double *, double *,
		  double *, double *, MPI_Comm);
void initbufs(int, int, double *, double *, double *, double *);
void getMinAvg(int, const double [], double *, double *);
int setNreps(int totlen, int minreps);

#define MAX_NTEST 20
static int nreps = 100;

static FILE *fp;

int main(int argc, char *argv[])
{
    int          *arrsizes, nsizes, i, dim, dim3;
    double       *sbuf, *rbuf, *scbuf, *rcbuf;
    double       tupack[MAX_NTEST], tvec[MAX_NTEST], tresize[MAX_NTEST],
	thvec[MAX_NTEST];
    MPI_Datatype vtype, rtype, htype, c3type;
    MPI_Comm     nodecomm, leadercomm, comm;
    int          nnodes, noderank, csize, crank, wrank, nt, ntests=10, max3d;

    MPI_Init(0,0);

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Get arguments */
    /* vecperf --size range [ --max3d dim ] */
    fp       = stdout; /* Default is output to stdout */
    arrsizes = BENV_GetSizes(argc, argv, "size", &nsizes);
    max3d = 1024;
    for (i=1; i<argc; i++) {
	char *v = argv[i];
	if (!v) continue;
	if (*v == '-') v++;
	if (*v == '-') v++;
	if (strcmp(v, "max3d") == 0) {
	    max3d = atoi(argv[++i]);
	}
	else if (strcmp(v, "o") == 0) {
	    fp = fopen(argv[++i], "w");
	    if (!fp) {
		fprintf(stderr, "Unable to open output file %s\n", argv[i]);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
	else {
	    if (wrank == 0) {
		fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
		fflush(stderr);
	    }
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
    }

    /* Find processes to use (internode preferred) */
    MPIX_Nodecomm_create(MPI_COMM_WORLD, &nodecomm, &leadercomm,
			 &noderank, &nnodes);
    if (nnodes > 1) comm = leadercomm;
    else            comm = nodecomm;

    MPI_Comm_size(comm, &csize);
    MPI_Comm_rank(comm, &crank);
    if (csize < 2) {
	if (wrank == 0) {
	    fprintf(stderr, "This test requires at least 2 processes\n");
	}
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Basic 2D vector */
    /* Run tests: loop over sizes, then datatype, then communication methods */
    if (crank == 0) {
	fprintf(fp, "2D Vector\n");
	fprintf(fp, "n\tuser\tuser(avg)\tvec\tvec(avg)\tresize\tresize(avg)\n");
    }
    for (i=0; i<nsizes; i++) {
	dim = arrsizes[i];
	/* Array is dim x dim */

	sbuf  = (double *)malloc(dim * dim * sizeof(double));
	rbuf  = (double *)malloc(dim * dim * sizeof(double));
	scbuf = (double *)malloc(dim * sizeof(double));
	rcbuf = (double *)malloc(dim * sizeof(double));
	if (!sbuf || !rbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d x %d\n", dim, dim);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (!scbuf || !rcbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d\n", dim);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	initbufs(dim*dim, dim, sbuf, rbuf, scbuf, rcbuf);
	nreps = setNreps(dim, 10);

	/* Create the datatypes */
	MPI_Type_vector(dim, 1, dim, MPI_DOUBLE, &vtype);
	MPI_Type_create_resized(MPI_DOUBLE, 0, dim*sizeof(double), &rtype);
	MPI_Type_commit(&vtype);
	MPI_Type_commit(&rtype);

	for (nt=0; nt<ntests; nt++) {
	    /* Communicate with user pack/unpack */
	    tupack[nt]  = timeUser(dim, dim, (const double *)sbuf, rbuf, scbuf, rcbuf, comm);

	    /* Communicate with datatypes */
	    tvec[nt]    = timeVec(vtype, (const double *)sbuf, rbuf, comm);

	    tresize[nt] = timeResize(rtype, dim, (const double *)sbuf, rbuf, comm);
	}

	if (crank == 0) {
	    double tumin, tuavg, tvmin, tvavg, trmin, travg;
	    getMinAvg(ntests, tupack, &tumin, &tuavg);
	    getMinAvg(ntests, tvec, &tvmin, &tvavg);
	    getMinAvg(ntests, tresize, &trmin, &travg);
	    fprintf(fp, "%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
		    dim, tumin, tuavg, tvmin, tvavg, trmin, travg);
	    fflush(fp);
	}

	free(sbuf);
	free(rbuf);
	free(scbuf);
	free(rcbuf);
	MPI_Type_free(&vtype);
	MPI_Type_free(&rtype);
    }

    /* 2D Vector of 3 elements */
    /* Run tests: loop over sizes, then datatype, then communication methods */
    if (crank == 0) {
	fprintf(fp, "2D Vector blocksize=3\n");
	fprintf(fp, "n\tuser\tuser(avg)\tvec\tvec(avg)\tresize\tresize(avg)\thvec\thvec(avg)\n");
    }
    for (i=0; i<nsizes; i++) {
	dim = arrsizes[i];
	/* Array is dim x dim */

	sbuf  = (double *)malloc(dim * dim * sizeof(double));
	rbuf  = (double *)malloc(dim * dim * sizeof(double));
	scbuf = (double *)malloc(3 * dim * sizeof(double));
	rcbuf = (double *)malloc(3 * dim * sizeof(double));
	if (!sbuf || !rbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d x %d\n", dim, dim);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (!scbuf || !rcbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d\n", 3*dim);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	initbufs(dim*dim, 3*dim, sbuf, rbuf, scbuf, rcbuf);
	nreps = setNreps(3*dim, 10);

	/* Create the datatypes */
	MPI_Type_vector(dim, 3, dim, MPI_DOUBLE, &vtype);
	MPI_Type_contiguous(3, MPI_DOUBLE, &c3type);
	MPI_Type_create_hvector(dim, 1, dim * sizeof(double), c3type, &htype);
	MPI_Type_create_resized(c3type, 0, dim*sizeof(double), &rtype);
	MPI_Type_commit(&vtype);
	MPI_Type_commit(&rtype);
	MPI_Type_commit(&htype);

	for (nt=0; nt<ntests; nt++) {
	    /* Communicate with user pack/unpack */
	    tupack[nt]  = timeUser3(dim, dim, (const double *)sbuf, rbuf, scbuf, rcbuf, comm);

	    /* Communicate with datatypes.  We can use the same code as
	       for the 1-element vector tests */
	    tvec[nt]    = timeVec(vtype, (const double *)sbuf, rbuf, comm);

	    thvec[nt]   = timeVec(htype, (const double *)sbuf, rbuf, comm);

	    tresize[nt] = timeResize(rtype, dim, (const double *)sbuf, rbuf, comm);
	}

	if (crank == 0) {
	    double tumin, tuavg, tvmin, tvavg, trmin, travg, thvmin, thvavg;
	    getMinAvg(ntests, tupack, &tumin, &tuavg);
	    getMinAvg(ntests, tvec, &tvmin, &tvavg);
	    getMinAvg(ntests, tresize, &trmin, &travg);
	    getMinAvg(ntests, thvec, &thvmin, &thvavg);
	    fprintf(fp, "%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
		    dim, tumin, tuavg, tvmin, tvavg, trmin, travg, thvmin,
		    thvavg);
	    fflush(fp);
	}

	free(sbuf);
	free(rbuf);
	free(scbuf);
	free(rcbuf);
	MPI_Type_free(&c3type);
	MPI_Type_free(&vtype);
	MPI_Type_free(&rtype);
	MPI_Type_free(&htype);
    }

    /* 2D plane of 3D domain */
    /* Run tests: loop over sizes, then datatype, then communication methods */
    if (crank == 0) {
	fprintf(fp, "2D y-z plane in 3D size n-2 x n-2\n");
	fprintf(fp, "n\tuser\tuser(avg)\tvec\tvec(avg)\n");
    }
    for (i=0; i<nsizes; i++) {
	dim = arrsizes[i];
	dim3 = dim;
	if (dim3 > 32) dim3 = 32;
	if (dim > max3d) continue;
	/* Array is dim x dim x dim3 */

	sbuf  = (double *)malloc(dim * dim * dim3 * sizeof(double));
	rbuf  = (double *)malloc(dim * dim * dim3 * sizeof(double));
	scbuf = (double *)malloc(dim * dim * sizeof(double));
	rcbuf = (double *)malloc(dim * dim * sizeof(double));
	if (!sbuf || !rbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d x %d x %d\n",
		    dim, dim, dim3);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (!scbuf || !rcbuf) {
	    fprintf(stderr, "Unable to allocate arrays of %d\n", dim*dim);
	    fflush(stderr);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	initbufs(dim*dim*dim3, dim*dim, sbuf, rbuf, scbuf, rcbuf);
	nreps = setNreps(dim*dim, 10);

	/* Create the datatypes */
	/* Communication the dim-2 x dim-2 piece of the face */
	MPI_Type_vector(dim-2, 1, dim3, MPI_DOUBLE, &vtype);
	MPI_Type_create_hvector(dim-2, 1, dim3 * dim * sizeof(double), vtype,
				&htype);
	MPI_Type_commit(&htype);

	for (nt=0; nt<ntests; nt++) {
	    /* Communicate with user pack/unpack */
	    tupack[nt]  = timeUser3D(dim, dim, dim3,
			       (const double *)sbuf, rbuf, scbuf, rcbuf, comm);

	    /* Communicate with datatypes.  We can use the same code as
	       for the 2D vector tests */
	    tvec[nt]    = timeVec(htype, (const double *)sbuf, rbuf, comm);
	}

	if (crank == 0) {
	    double tumin, tuavg, tvmin, tvavg;
	    getMinAvg(ntests, tupack, &tumin, &tuavg);
	    getMinAvg(ntests, tvec, &tvmin, &tvavg);
	    fprintf(fp, "%d\t%.2e\t%.2e\t%.2e\t%.2e\n",
		    dim, tumin, tuavg, tvmin, tvavg);
	    fflush(fp);
	}

	free(sbuf);
	free(rbuf);
	free(scbuf);
	free(rcbuf);
	MPI_Type_free(&vtype);
	MPI_Type_free(&htype);
    }

    MPI_Finalize();
    return 0;
}

double timeUser(int dim1, int dim2, const double *restrict sbuf,
		double *restrict rbuf, double *restrict scbuf,
		double *restrict rcbuf, MPI_Comm comm)
{
    const double * restrict ps;
    double * restrict pr;
    double t0;
    int crank, i, r;

    MPI_Comm_rank(comm, &crank);

    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (r=0; r<nreps; r++) {
	if (crank == 0) {
	    ps = sbuf;
	    for (i=0; i<dim1; i++) {
		scbuf[i] = *ps; ps += dim2;
	    }
	    MPI_Send(scbuf, dim1, MPI_DOUBLE, 1, 0, comm);
	    MPI_Recv(rcbuf, dim1, MPI_DOUBLE, 1, 0, comm, MPI_STATUS_IGNORE);
	    pr = rbuf;
	    for (i=0; i<dim1; i++) {
		*pr = rcbuf[i]; pr += dim2;
	    }
	}
	else if (crank == 1) {
	    MPI_Recv(rcbuf, dim1, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
	    pr = rbuf;
	    for (i=0; i<dim1; i++) {
		*pr = rcbuf[i]; pr += dim2;
	    }
	    ps = sbuf;
	    for (i=0; i<dim1; i++) {
		scbuf[i] = *ps; ps += dim2;
	    }
	    MPI_Send(scbuf, dim1, MPI_DOUBLE, 0, 0, comm);
	}
    }
    t0 = MPI_Wtime() - t0;

    return t0/nreps;
}

double timeUser3(int dim1, int dim2, const double *restrict sbuf,
		 double *restrict rbuf, double *restrict scbuf,
		 double *restrict rcbuf, MPI_Comm comm)
{
    const double * restrict ps;
    double * restrict pr;
    double t0;
    int crank, i, j, r;

    MPI_Comm_rank(comm, &crank);

    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (r=0; r<nreps; r++) {
	if (crank == 0) {
	    ps = sbuf;
	    j  = 0;
	    for (i=0; i<dim1; i++) {
		scbuf[j] = ps[0]; scbuf[j+1] = ps[1]; scbuf[j+2] = ps[2];
		ps += dim2;
		j += 3;
	    }
	    MPI_Send(scbuf, 3*dim1, MPI_DOUBLE, 1, 0, comm);
	    MPI_Recv(rcbuf, 3*dim1, MPI_DOUBLE, 1, 0, comm, MPI_STATUS_IGNORE);
	    pr = rbuf;
	    j  = 0;
	    for (i=0; i<dim1; i++) {
		pr[0] = rcbuf[j];
		pr[1] = rcbuf[j+1];
		pr[2] = rcbuf[j+2];
		pr += dim2;
		j  += 3;
	    }
	}
	else if (crank == 1) {
	    MPI_Recv(rcbuf, 3*dim1, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
	    pr = rbuf;
	    j  = 0;
	    for (i=0; i<dim1; i++) {
		pr[0] = rcbuf[j];
		pr[1] = rcbuf[j+1];
		pr[2] = rcbuf[j+2];
		pr += dim2;
		j  += 3;
	    }
	    ps = sbuf;
	    j  = 0;
	    for (i=0; i<dim1; i++) {
		scbuf[j] = ps[0]; scbuf[j+1] = ps[1]; scbuf[j+2] = ps[2];
		ps += dim2;
		j += 3;
	    }
	    MPI_Send(scbuf, 3*dim1, MPI_DOUBLE, 0, 0, comm);
	}
    }
    t0 = MPI_Wtime() - t0;

    return t0/nreps;
}

double timeUser3D(int dim1, int dim2, int dim3,
		  const double *sbuf, double *rbuf,
		  double *scbuf, double *rcbuf, MPI_Comm comm)
{
    const double * restrict ps;
    double * restrict pr;
    double t0;
    int crank, i, j, k, r;

    MPI_Comm_rank(comm, &crank);

    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (r=0; r<nreps; r++) {
	if (crank == 0) {
	    k = 0;
	    for (i=1; i<dim2-1; i++) {
		ps = sbuf + dim1 * (1 + i * dim3);
		for (j=1; j<dim1-1; j++) {
		    scbuf[k++] = *ps; ps += dim3;
		}
	    }
	    MPI_Send(scbuf, k, MPI_DOUBLE, 1, 0, comm);
	    MPI_Recv(rcbuf, k, MPI_DOUBLE, 1, 0, comm, MPI_STATUS_IGNORE);
	    pr = rbuf;
	    for (i=0; i<dim1; i++) {
		*pr = rcbuf[i]; pr += dim2;
	    }
	    k = 0;
	    for (i=1; i<dim2-1; i++) {
		pr = rbuf + dim1 * (1 + i * dim3);
		for (j=1; j<dim1-1; j++) {
		    *pr = rcbuf[k++]; pr += dim3;
		}
	    }
	}
	else if (crank == 1) {
	    MPI_Recv(rcbuf, (dim1-1)*(dim2-1), MPI_DOUBLE, 0, 0, comm,
		     MPI_STATUS_IGNORE);
	    k = 0;
	    for (i=1; i<dim2-1; i++) {
		pr = rbuf + dim1 * (1 + i * dim3);
		for (j=1; j<dim1-1; j++) {
		    *pr = rcbuf[k++]; pr += dim3;
		}
	    }
	    k = 0;
	    for (i=1; i<dim2-1; i++) {
		ps = sbuf + dim1 * (1 + i * dim3);
		for (j=1; j<dim1-1; j++) {
		    scbuf[k++] = *ps; ps += dim3;
		}
	    }
	    MPI_Send(scbuf, k, MPI_DOUBLE, 0, 0, comm);
	}
    }
    t0 = MPI_Wtime() - t0;

    return t0/nreps;
}


/* Communicate with datatypes */
double timeVec(MPI_Datatype vtype, const double *sbuf, double *rbuf,
		   MPI_Comm comm)
{
    int crank, r;
    double t0;
    MPI_Comm_rank(comm, &crank);

    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (r=0; r<nreps; r++) {
	if (crank == 0) {
	    MPI_Send(sbuf, 1, vtype, 1, 0, comm);
	    MPI_Recv(rbuf, 1, vtype, 1, 0, comm, MPI_STATUS_IGNORE);
	}
	else if (crank == 1) {
	    MPI_Recv(rbuf, 1, vtype, 0, 0, comm, MPI_STATUS_IGNORE);
	    MPI_Send(sbuf, 1, vtype, 0, 0, comm);
	}
    }
    t0 = MPI_Wtime() - t0;
    return t0/nreps;
}

double timeResize(MPI_Datatype rtype, int dim, const double *sbuf, double *rbuf,
		      MPI_Comm comm)
{
    int crank, r;
    double t0;
    MPI_Comm_rank(comm, &crank);

    MPI_Barrier(comm);
    t0 = MPI_Wtime();
    for (r=0; r<nreps; r++) {
	if (crank == 0) {
	    MPI_Send(sbuf, dim, rtype, 1, 0, comm);
	    MPI_Recv(rbuf, dim, rtype, 1, 0, comm, MPI_STATUS_IGNORE);
	}
	else if (crank == 1) {
	    MPI_Recv(rbuf, dim, rtype, 0, 0, comm, MPI_STATUS_IGNORE);
	    MPI_Send(sbuf, dim, rtype, 0, 0, comm);
	}
    }
    t0 = MPI_Wtime() - t0;
    return t0/nreps;
}

void initbufs(int dim1, int dim2, double *sbuf, double *rbuf,
	      double *scbuf, double *rcbuf)
{
    int i;

    for (i=0; i<dim1; i++) {
	sbuf[i] = 0;
	rbuf[i] = 0;
    }
    for (i=0; i<dim2; i++) {
	scbuf[i] = 0;
	rcbuf[i] = 0;
    }
}

void getMinAvg(int len, const double v[], double *vmin, double *vavg)
{
    double minval, avgval;
    int    i;

    if (len <= 0) {
	*vmin = 0;
	*vavg = 0;
	return;
    }
    minval = v[0];
    avgval = v[0];
    for (i=1; i<len; i++) {
	if (v[i] < minval) minval = v[i];
	avgval += v[i];
    }
    *vavg = avgval / len;
    *vmin = minval;
}

int setNreps(int totlen, int minreps)
{
    int n = minreps;
    if (totlen < 1000) n = 100;
    if (totlen < 1000) n = 20;
    return n;
}
