/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.

   This version defines two entry points -- with
   and without appended underscores, so it *should*
   automagically link with FORTRAN */

#include <sys/time.h>

/* Since these are intended for Fortran, we define the prototypes here */
double mysecond(void);
double mysecond_(void);

double mysecond(void)
{
/* struct timeval { long        tv_sec;
            long        tv_usec;        };

struct timezone { int   tz_minuteswest;
             int        tz_dsttime;      };     */

        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

double mysecond_(void) {return mysecond();}

