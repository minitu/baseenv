#ifndef NCARTSHIM_H_INCLUDED
#define NCARTSHIM_H_INCLUDED

/* This routine allows programs that use the "shims" for MPI_Cart_create etc.
   to confirm that in fact the shim routines were called. */
int MPIX_CartShim_info(int *createcalls, int *subcalls);
#endif
