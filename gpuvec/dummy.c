#include "gpuvec.h"

#pragma acc routine (dummy) seq
#pragma acc routine (dummyab) seq
#pragma acc routine (dummyabc) seq

int dummy(BASETYPE *restrict a, BASETYPE *restrict b, BASETYPE *restrict c,
	  BASETYPE *restrict d, BASETYPE *restrict e,
	  BASETYPE aa[LEN2][LEN2], BASETYPE bb[LEN2][LEN2],
	  BASETYPE cc[LEN2][LEN2], BASETYPE s)
{
    // --  called in each loop to make all computations appear required
    return 0;
}
int dummyab(BASETYPE *restrict a, BASETYPE *restrict b)
{
    // --  called in each loop to make all computations appear required
    return 0;
}
int dummyabc(BASETYPE *restrict a, BASETYPE *restrict b, BASETYPE *restrict c)
{
    // --  called in each loop to make all computations appear required
    return 0;
}
