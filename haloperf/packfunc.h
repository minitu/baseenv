double pksimple( double *inbuf, double *outbuf, int count, int stride );
double pkconst( const double *inbuf, double *outbuf, int count, int stride );
double pkrestrict( const double *restrict inbuf, double *restrict outbuf, 
		   int count, int stride );
double pksplit( const double *restrict inbuf, double *restrict outbuf, 
		int count, int stride );



