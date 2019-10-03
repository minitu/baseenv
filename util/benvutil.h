#ifndef BENVUTIL_H_INCLUDED
#define BENVUTIL_H_INCLUDED

int *BENV_GetSizes(int argc, char *argv[], const char *argname, int *nsizes);
int *BENV_GetSizesArith(int start, int end, int stride, int *nsizes);
int *BENV_GetSizesMult(int start, int end, double factor, int *nsizes);
int *BENV_GetSizesMultDelta(int start, int end, double factor, int delta,
			    int ndelta, int *nsizes);

int BENV_GetNrepsEst(double deltat, double latency, double bandwidth,
		     int nmsgs, int msglen, int minreps);

#endif /* BENVUTIL_H */

