#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/* This is the read version of ioda */

#ifndef MYSCRATCHDIR
#ifdef _CRAYC
#define MYSCRATCHDIR "/scratch/training/trnxxx/"
#error edit the string for scratch dir to use your user name and delete this line
#else
#define MYSCRATCHDIR "./"
#error fix the scratchdir and delete this line
#endif
#endif

int verifyRead = 1;

void myAbort(int err, const char *msg);

int main(int argc, char *argv[])
{
    int iarrayOfSizes[2], iarrayOfSubsizes[2], iarrayOfStarts[2], ilocal_size;
    int nproc[2], periods[2], icoord[2];
    int m, n, i, j, wsize, wrank, crank, ndims, lrows, lcols, grow, gcol, err;
    MPI_Datatype filetype;
    MPI_File     fh;
    MPI_Comm     cartcomm;
    MPI_Info     info0, info3;
    double       t, topen, tread, tclose, wrate;
    double       *local_array;
    double       *ref_array;
    char         nstripesStr[12], stripeUnitStr[12];
    int          nstripes = -1;
    int          stripeUnit = -1;

    MPI_Init(0,0);

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    /* Get global array size */
    m = n = 128;      /* Set default size */

    /* ioda [ n ] [ m ] [ nstripes ] [ stripeunit ]*/
    if (argc > 0) {
	if (argc > 1) m = atoi(argv[1]);
	if (argc > 2) n = atoi(argv[2]);
	if (argc > 3) nstripes = atoi(argv[3]);
	if (argc > 4) stripeUnit = atoi(argv[4]);
	if (argc > 5) {
	    if (wrank == 0)
		fprintf(stderr,"Unrecognized argument %s\n", argv[5]);
	    MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
    if (wrank == 0) printf("Matrix is [%d,%d]; file dir = %s\n", m, n, MYSCRATCHDIR );

    /* The default number of stripes = totalsize/1M */
    if (nstripes < 0) {
	nstripes = n * m * sizeof(double) / (1024*1024);
	if (nstripes < 1) nstripes = 1;
    }
    if (wrank == 0) printf("nstripes = %d, stripeUnit = %d\n", nstripes, stripeUnit);

    /* Use topology routines to get decomposition and coordinates */
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    nproc[0] = 0; nproc[1] = 0;
    ndims = 2;
    MPI_Dims_create(wsize, ndims, nproc);
    periods[0] = 0; periods[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, nproc, periods, 1, &cartcomm);
    MPI_Comm_rank(cartcomm, &crank);
    MPI_Cart_coords(cartcomm, crank, ndims, icoord);

    iarrayOfSizes[0]    = m;
    iarrayOfSizes[1]    = n;
    iarrayOfSubsizes[0] = m/nproc[0];
    iarrayOfSubsizes[1] = n/nproc[1];
    iarrayOfStarts[0]   = icoord[0] * iarrayOfSubsizes[0];
    iarrayOfStarts[1]   = icoord[1] * iarrayOfSubsizes[1];

    /* Initialize my block of the data */
    ilocal_size = iarrayOfSubsizes[0] * iarrayOfSubsizes[1];
    lrows = iarrayOfSubsizes[0];
    lcols = iarrayOfSubsizes[1];
    local_array = (double *)malloc(lrows*lcols*sizeof(double));
    ref_array   = (double *)malloc(lrows*lcols*sizeof(double));
    gcol  = iarrayOfStarts[1];
    grow = iarrayOfStarts[0];
    for (i=0; i<lrows; i++) {
	for (j=0; j<lcols; j++) {
	    ref_array[j*lrows+i]   = (grow+i) + (gcol+j)*m;
	    local_array[j*lrows+i] = -1;
	}
    }

    /* Fortran order simply means the data is stored by columns */
    MPI_Type_create_subarray(ndims, iarrayOfSizes, iarrayOfSubsizes,
			     iarrayOfStarts, MPI_ORDER_FORTRAN, MPI_DOUBLE,
			     &filetype);
    MPI_Type_commit(&filetype);

    info0 = MPI_INFO_NULL;
    info3 = MPI_INFO_NULL;
    if (nstripes > 0 || stripeUnit > 0) {
	MPI_Info_create(&info0);
	if (nstripes > 0) {
	    snprintf(nstripesStr, sizeof(nstripesStr), "%d", nstripes);
	    MPI_Info_set(info0, "striping_factor", nstripesStr);
	    MPI_Info_set(info0, "cb_nodes", nstripesStr);
	}
	if (stripeUnit > 0) {
	    snprintf(stripeUnitStr, sizeof(stripeUnitStr), "%d", stripeUnit);
	    MPI_Info_set(info0, "striping_unit", stripeUnitStr);
	}
	MPI_Info_dup(info0, &info3);
	MPI_Info_set(info3, "romio_no_indep_rw", "true");

	/* Other hints to consider:
	   direct_io=true

	   The default cb_buffer_size is 16777216 , but is overridden by the
	   striping unit, which is smaller by default.
	*/
    }

    /* level - 3 */
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    err = MPI_File_open(cartcomm, MYSCRATCHDIR "testfile-3.out",
			MPI_MODE_RDONLY, info3, &fh);
    topen = MPI_Wtime() - t;
    if (err != MPI_SUCCESS) myAbort(err, "open testfile-3.out");

    MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    err = MPI_File_read_all(fh, local_array, ilocal_size, MPI_DOUBLE,
			     MPI_STATUS_IGNORE);
    tread = MPI_Wtime() - t;
    if (err != MPI_SUCCESS) myAbort(err, "collective read");

    err = MPI_File_close(&fh);
    tclose = MPI_Wtime() - t;
    if (err != MPI_SUCCESS) myAbort(err, "close testfile-3.out");

    if (verifyRead) {
	int errs = 0;
	for (i=0; i<lrows; i++) {
	    for (j=0; j<lcols; j++) {
		if (ref_array[j*lrows+i] != local_array[j*lrows+i]) {
		    errs++;
		    if (errs < 10) {
			int idx = j*lrows+i;
			printf("array[%d] = %e, expected %e\n", idx,
			       local_array[idx], ref_array[idx]);
		    }
		}
	    }
	}
    }

    MPI_Allreduce(MPI_IN_PLACE, &topen, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tread, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tclose, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (tread > 0)
	wrate = (double)m * (double)n * sizeof(double)/tread;
    if (wrank == 0)
	printf("%d\t[%d,%d]\t%d\t%.2e\t%.2e\t%.2e\t%.2e\n", wsize, m, n, nstripes, topen,
	       tread, tclose, wrate);

    /* level - 0 */
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    err = MPI_File_open(cartcomm, MYSCRATCHDIR "testfile-0.out",
			MPI_MODE_RDONLY, info0, &fh);
    topen = MPI_Wtime() - t;
    printf("err on open = %d\n", err);
    if (err != MPI_SUCCESS) myAbort(err, "open testfile-0.out");

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    gcol = iarrayOfStarts[1];
    grow = iarrayOfStarts[0];
    for (j=0; j<lcols; j++) {
	MPI_Offset offset =
	    ((MPI_Offset)(grow) + (MPI_Offset)(gcol+j)*m) * sizeof(double);
	err = MPI_File_read_at(fh, offset, local_array+j*lrows, lrows, MPI_DOUBLE,
				MPI_STATUS_IGNORE);
	if (err != MPI_SUCCESS) myAbort(err, "read at");
    }
    tread = MPI_Wtime() - t;

    err = MPI_File_close(&fh);
    tclose = MPI_Wtime() - t;
    if (err != MPI_SUCCESS) myAbort(err, "close testfile-0");

    MPI_Allreduce(MPI_IN_PLACE, &topen, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tread, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tclose, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (tread > 0)
	wrate = (double)m * (double)n * sizeof(double)/tread;
    if (wrank == 0)
	printf("%d\t[%d,%d]\t%d\t%.2e\t%.2e\t%.2e\t%.2e\n", wsize, m, n, nstripes, topen,
	       tread, tclose, wrate);

    if (info0 != MPI_INFO_NULL) {
	MPI_Info_free(&info0);
	MPI_Info_free(&info3);
    }
    free(local_array);
    MPI_Finalize();
    return 0;
}

void myAbort(int err, const char *msg)
{
    char str[MPI_MAX_ERROR_STRING];
    int  rlen;

    str[0] = 0;
    if (err != MPI_SUCCESS)
	MPI_Error_string(err, str, &rlen);
    fprintf(stderr,"%s:%s\n", msg, str);
    MPI_Abort(MPI_COMM_WORLD, 1);
}
