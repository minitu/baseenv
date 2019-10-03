#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "benvutil.h"

typedef enum { UNKNOWN, LIST, RANGE, RANGE_ADDITIVE, RANGE_MULT,
	       RANGE_MULT_DELTA } sizeval_t;

#define MAXLISTVALS 256
static int scanInt(const char **str);
static void errSizeArg(const char *, char *, char *);
/*
 * Look through the arg list for -argname value, where value can be:
 *   a,b,...,c (comma separated list)
 *   a:b or a:b:c (range with additive stride)
 *   a:b*c (range with multiplicative stride)
 *   a:b*c+d:n (range with multiplicative stride, with arithmetic
 *              distribution around points: that is, take the values from
 *              a:b*c and add (-nd, -(n-1)d,... -d, d, ..., nd) points,
 *              staying within the interval [a,b]
 *
 * a and b make have postfix k or K for X 1024, m or M for 1024*1024.
 * On error, nsizes == -1 and return is NULL.
 * User should free(ptr) the value returned by this routine.
 *
 * Enhancement: make it possible to include either range type within
 * a list, e.g., 2,4,7:23,32.
 *
 */

#define KIBI 1024
#define MEBI (KIBI*KIBI)

/* QUERY: extract just the parse arg, return the values of a,b,c,d,n
   and the type; if list, return the list val */

/*@ BENV_GetSizes - Return an array of integers from the argument list

 Input Parameters:
+ argc - argument count
- argname - argument name

Input/Output Parameter:
. argv - argument vector on input.  If 'argname' is found, the values
 of argv for that argument and the following argument are set to null.

 Output Parameter:
. nsizes - number of values returned

 Return value:
 Pointer to an array, allocated with 'malloc', containing the integers
 On error, 'nsizes == -1' and the return value is 'NULL'.
 The user should 'free(ptr)' where 'ptr' is the value returned by this routine.

 Notes:
This routine looks for a pair of arguments of the form
.vb
      ... --argname list
.ve
where 'list' may be either a comma-separated list of values or a range with
an optional stride.  Specifically,
.vb
   a,b,...,c (comma separated list)
   a:b or a:b:c (range with additive stride)
   a:b*c (range with multiplicative stride)
   a:b*c+d:n (range with multiplicative stride, with arithmetic
              distribution around points: that is, take the values from
              a:b*c and add (-nd, -(n-1)d,... -d, d, ..., nd) points,
              staying within the interval [a,b]
.ve
Here, 'a' and 'b' may have postfix 'k' or 'K' for times 1024 and 'm' or 'M'
for time 1024*1024.

This routine sets found arguments to null to permit other routines to process
other arguments; this requires all routines to permit a 'NULL' argument for
a value in 'argv'.
  @*/
int *BENV_GetSizes(int argc, char *argv[], const char *argname, int *nsizes)
{
    sizeval_t sizetype;
    int       i, j;
    char      *p;
    int       start, end, stride=1;
    int       ival;
    int       delta=0, ndelta=0; /* for the arithmetic around multiplicative
				    strides */
    double    fval, fstride;
    int       nvals, listvals[MAXLISTVALS], *sptr;

    *nsizes = -1;

    for (i=1; i<argc; i++) {
	p = argv[i];
	if (!p || !*p) continue;
	if (*p++ != '-') continue;
	if (*p == '-') p++;
	if (strcmp(p, argname) != 0) continue;
	/* Found the argument.  Process it and exit */
	p = argv[i+1];
	sizetype = UNKNOWN;
	nvals    = 0;
	while (*p) {
	    ival = scanInt((const char **)&p);
	    fval = ival;
	    if (sizetype == RANGE_MULT && *p == '.') {
		/* Handle the case of a decimal stride for the multiplicative
		   fraction. For best precision, use strtod, but this
		   is good enough for the need here. */
		double frac = 1e-1;
		p++;
		while (*p && isdigit(*p)) {
		    fval += (*p - '0') * frac;
		    frac *= 1e-1;
		    p++;
		}
		/* No postfix (kKmM) for stride */
	    }

	    switch (sizetype) {
	    case UNKNOWN:
		if (!*p || *p == ',') {
		    if (nvals >= MAXLISTVALS) {
			errSizeArg("list", argv[i], argv[i+1]);
			return 0;
		    }
		    sizetype = LIST; listvals[nvals++] = ival;
		}
		else if (*p == ':') {
		    sizetype = RANGE;
		    start    = ival;
		}
		else {
		    errSizeArg("size", argv[i], argv[i+1]);
		    return 0;
		}
		if (*p) p++;
		break;
	    case LIST:
		if (!*p || *p == ',') {
		    if (nvals >= MAXLISTVALS) {
			errSizeArg("list", argv[i], argv[i+1]);
			return 0;
		    }
		    listvals[nvals++] = ival;
		}
		else {
		    errSizeArg("list", argv[i], argv[i+1]);
		    return 0;
		}
		if (*p) p++;
		break;
	    case RANGE:
		if (!*p || *p == ':') {
		    sizetype = RANGE_ADDITIVE;
		    end = ival;
		}
		else if (*p == '*') {
		    sizetype = RANGE_MULT;
		    end = ival;
		}
		else {
		    errSizeArg("range", argv[i], argv[i+1]);
		    return 0;
		}
		if (*p) p++;
		break;
	    case RANGE_ADDITIVE:
		if (!*p) {
		    stride = ival;
		}
		else {
		    errSizeArg("range", argv[i], argv[i+1]);
		    return 0;
		}
		break;
	    case RANGE_MULT:
		if (*p == '+') {
		    fstride  = fval;
		    p++;
		    sizetype = RANGE_MULT_DELTA;
		}
		else if (!*p) {
		    fstride = fval;
		}
		else {
		    errSizeArg("range", argv[i], argv[i+1]);
		    return 0;
		}
		break;
	    case RANGE_MULT_DELTA:
		delta   = ival;
		if (*p != ':') {
		    errSizeArg("range", argv[i], argv[i+1]);
		    return 0;
		}
		else {
		    p++;
		    ndelta = scanInt((const char **)&p);
		}
	    } /* switch */
	} /* while *p */

	switch (sizetype) {
	case LIST:
	    *nsizes = nvals;
	    sptr = (int *)malloc(nvals * sizeof(int));
	    if (!sptr) {
		fprintf(stderr, "Could not allocate %d words for size %s\n",
			nvals, argname);
		return 0;
	    }
	    for (j=0; j<nvals; j++) sptr[j] = listvals[j];
	    break;
	case RANGE:
	case RANGE_ADDITIVE:
	    sptr = BENV_GetSizesArith(start, end, stride, nsizes);
	    if (!sptr) {
		fprintf(stderr, "Could not allocate %d words for size %s\n",
			nvals, argname);
		return 0;
	    }
	    break;
	case RANGE_MULT:
	    sptr = BENV_GetSizesMult(start, end, fstride, nsizes);
	    if (!sptr) {
		fprintf(stderr, "Could not allocate %d words for size %s\n",
			nvals, argname);
		return 0;
	    }
	    break;
	case RANGE_MULT_DELTA:
	    sptr = BENV_GetSizesMultDelta(start, end, fstride, delta, ndelta,
					  nsizes);
	    if (!sptr) {
		fprintf(stderr, "Could not allocate %d words for size %s\n",
			nvals, argname);
		return 0;
	    }
	    break;
	default:
	    fprintf(stderr, "Malformed size argument %s %s\n",
			    argv[i], argv[i+1]);
	    return 0;
	    break;
	}
	argv[i]   = 0;
	argv[i+1] = 0;
	return sptr;
    }

    return 0;
}

/*
 * Routines to generate size arrays, given start, end, and increment
 * information.  Useful for providing a default is no size option provided.
 */
int *BENV_GetSizesArith(int start, int end, int stride, int *nsizes)
{
    int nvals, *sptr, j, ival;
    nvals = 1 + (end - start) / stride;
    sptr = (int *)malloc(nvals * sizeof(int));
    if (!sptr) return 0;

    ival = start;
    for (j=0; j<nvals; j++) {
	sptr[j] = ival;
	ival += stride;
    }
    *nsizes = nvals;
    return sptr;
}

int *BENV_GetSizesMult(int start, int end, double factor, int *nsizes)
{
    int nvals, *sptr, ival;

    nvals = 1;
    ival = start;
    if (factor <= 1) {
	fprintf(stderr, "Multiplicative range requires factor > 1\n");
	return 0;
    }
    if ((int)(ival * factor) == ival) {
	fprintf(stderr,
		"Multiplicative factor must satisfy factor*start>=start+1\n");
	return 0;
    }

    /* Determine number of elements */
    while (ival <= end) { ival *= factor; nvals++; }

    sptr = (int *)malloc(nvals * sizeof(int));
    if (!sptr) return 0;

    ival  = start;
    nvals = 0;
    while (ival <= end) {
	sptr[nvals++] = ival;
	ival *= factor;
    }
    *nsizes = nvals;
    return sptr;
}

int *BENV_GetSizesMultDelta(int start, int end, double factor, int delta,
			    int ndelta, int *nsizes)
{
    int nvals, *sptr, ival;

    nvals = 1;
    ival = start;
    if (factor <= 1) {
	fprintf(stderr, "Multiplicative range requires factor > 1\n");
	return 0;
    }
    if ((int)(ival * factor) == ival) {
	fprintf(stderr,
		"Multiplicative factor must satisfy factor*start>=start+1\n");
	return 0;
    }
    if (ndelta < 0 && delta < 0) {
	fprintf(stderr, "Arithmetic delta and count must be positive\n");
	return 0;
    }

    /* Determine number of elements.  This is a slight overestimate,
       because at a and b, there are only ndelta instead of 2*ndelta;
       values that are less than the current top value are also ignored. */
    while (ival <= end) { ival *= factor; nvals += 2*ndelta; }

    sptr = (int *)malloc(nvals * sizeof(int));
    if (!sptr) return 0;

    /* Add the initial value to simplify the tests */
    nvals         = 0;
    sptr[nvals++] = start;
    ival          = start;
    while (ival <= end) {
	int nd, iival;
	iival = ival-ndelta*delta;
	for (nd=-ndelta; nd<=ndelta; nd++) {
	    if (iival >= start && iival <= end) {
		/* some values might for an insert not at the end */
		if (iival > sptr[nvals-1])
		    sptr[nvals++] = iival;
		else {
		    int ll = nvals-1;
		    while (ll > 0 && sptr[ll] > iival) ll--;
		    /* Check for a duplicate */
		    if (sptr[ll] != iival) {
			/* Make room in the list */
			for (int idx=nvals; idx>ll; idx--)
			    sptr[idx+1] = sptr[idx];
			nvals++;
			sptr[ll+1] = iival;
		    }
		}
	    }
	    iival += delta;
	}
	ival *= factor;
    }
    *nsizes = nvals;
    return sptr;
}

/*
 * Scan an int, which may be followed by k,K,m, or M, which are 1024 or 1024^2
 * Update the passed argument to point to the next character.
 */
static int scanInt(const char **str)
{
    int ival = 0;
    const char *p = *str;

    while (*p && isdigit(*p)) {
	ival = (*p - '0') + 10*ival;   /* Assumes ASCII */
	p++;
    }
    if (*p) {
	if (*p == 'k' || *p == 'K') {
	    ival *= KIBI;
	    p++;
	}
	if (*p == 'm' || *p == 'M') {
	    ival *= MEBI;
	    p++;
	}
	/* Allow "i" or "I" to follow */
	if (*p == 'i' || *p == 'I') p++;
    }
    *str = p;
    return ival;
}

static void errSizeArg(const char *nm, char *a1, char *a2)
{
    fprintf(stderr, "Malformed %s argument %s %s\n", nm, a1, a2);
}


#ifdef BUILD_TEST
int main(int argc, char *argv[])
{
    int nsizes, *sptr, i;
    do {
	sptr = BENV_GetSizes(argc, argv, "size", &nsizes);
	printf("%d values:", nsizes); fflush(stdout);
	if (!sptr) break;
	for (i=0; i<nsizes; i++)
	    printf("%d%s", sptr[i], (i!=nsizes-1) ? "," : "\n");
	fflush(stdout);
	free(sptr);
    } while (1);

    return 0;
}
#endif
