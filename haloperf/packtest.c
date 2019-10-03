/* */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "haloperf.h"
#include "packfunc.h"

/* 
 * Perform some simple tests about the performance of contiguous and
 * non-contiguous copies.
 */
int main( int argc, char *argv[] )
{
    int stride, nwords, bufsize, i;
    double *inbuf, *outbuf, ttest, rate;

    MPI_Init( &argc, &argv );
    
    /* Set a typical stride and size for a 1000x1000 mesh */
    stride = 3000;
    nwords = 3000;
    if (argc > 1) stride = atoi(argv[1]);
    if (argc > 2) nwords = atoi(argv[2]);

    bufsize = nwords * stride;
    inbuf  = (double *)malloc( bufsize * sizeof(double) );
    outbuf = (double *)malloc( bufsize * sizeof(double) );
    if (!inbuf || !outbuf) {
	fprintf( stderr, "Unable to allocate memory!\n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    printf("Nwords = %d, stride = %d\n", nwords, stride);
    printf("%20s\ttime\trate (MB/s)\n", "Pack approach");
    
    /* Touch the buffers */
    for (i=0; i<bufsize; i++) {
	inbuf[i] = i;
	outbuf[i] = 0;
    }
    
    ttest = pksimple( inbuf, outbuf, nwords, stride );
    rate = nwords * 8 / ttest / 1e6;
    printf( "%20s\t%.2e\t%.2e\n", "Simple", ttest, rate );

    ttest = pkconst( (const double *)inbuf, outbuf, nwords, stride );
    rate = nwords * 8 / ttest / 1e6;
    printf( "%20s\t%.2e\t%.2e\n", "Const", ttest, rate );

    ttest = pkrestrict( (const double *restrict)inbuf, 
			(double *restrict) outbuf, nwords, stride );
    rate = nwords * 8 / ttest / 1e6;
    printf( "%20s\t%.2e\t%.2e\n", "Restrict", ttest, rate );

    ttest = pksplit( (const double *restrict)inbuf, 
			(double *restrict) outbuf, nwords, stride );
    rate = nwords * 8 / ttest / 1e6;
    printf( "%20s\t%.2e\t%.2e\n", "Split", ttest, rate );

    MPI_Finalize();
    return 0;
}
