AC_PREREQ(2.69)
AC_INIT([perftest],
        2.0,
        [wgropp@illinois.edu],
        [perftest],
        [http://www.cs.illinois.edu/~wgropp])

AC_CONFIG_HEADER([mpptestconf.h])

AC_ARG_ENABLE(echo,
[--enable-echo    - Turn on echoing from within the configure script],
set -x)

# Set the C and Fortran FLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for this package
: ${CFLAGS=""}
: ${FCFLAGS=""}

# Optimization should be in COPT and FCOPT
AC_SUBST([COPT])
AC_SUBST([FCOPT])

PAC_PROG_MPICC
if test -n "$MPICC" ; then CC=$MPICC ; fi
AC_PROG_CC
PAC_MPI_CCWORKS
# Use COPT for optimization-only flags for C

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/time.h stdlib.h string.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_TIME
AC_STRUCT_TM
AC_C_CONST
AC_C_RESTRICT

AC_SEARCH_LIBS([sqrt],[m])

MPICC=$CC

PAC_PROG_MPIFORT
if test -n "$MPIFORT" ; then FC=$MPIFORT ; fi
AC_PROG_FC
MPIFORT=$FC

# If we have perl, we can build the tool to extract data from the gnuplot
# output
AC_PATH_PROG(PERL,perl)

# Look for graphics programs
AC_CHECK_PROGS(GRAPHICS_PGM,cit gnuplot,"none")

#
# Select the location and name for the base environment library
have_seqbegin=0
if test -n "$BASEENV_LIB" ; then
   AC_SUBST([BASEENV_LIB])
   if test "$BUILDING_SEQ" = "yes" ; then
       if test -s ../seq/seq.h ; then
           have_seqbegin=1
	   CPPFLAGS="$CPPFLAGS -I../seq"
       fi
   else
       # Only do this if configuring separately, or if seq is already built
       if test -n "$BASEENV_LIBNAME" ; then
           AC_SEARCH_LIBS([seqBegin],[$BASEENV_LIBNAME],have_seqbegin=1)
	   if test "$have_seqbegin" = 1 ; then
	      AC_CHECK_HEADER([seq.h],[],[have_seqbegin=0])
	      # If the header is not found, consider adding the BASEENV
	      # include path to CPPFLAGS and running this again.
	   fi
       fi
   fi
fi
if test "$have_seqbegin" = 1 ; then
   AC_DEFINE([HAVE_SEQBEGIN],[1],[Define if seqBegin and seqEnd are available])
fi

# Normally, this configure should be invoked by the top-level baseenv
# builddir.  If not, this code properly defines the master_top_builddir
# and master_top_srcdir to build this as a standalone package
if test "X$master_top_builddir" = "X" ; then
   master_top_builddir=$ac_top_builddir
fi
if test "X$master_top_srcdir" = "X" ; then
   master_top_srcdir=$ac_top_srcdir
fi

AC_SUBST([master_top_builddir])
AC_SUBST([master_top_srcdir])

# Find a working mkdir
AC_PROG_MKDIR_P
# Find a working install
AC_PROG_INSTALL
# In case the documentation targets are used, find doctext and attempt to
# find the source for the doctext LaTeX style files.  Use "false" if
# doctext is not found
AC_PATH_PROG(DOCTEXT,[doctext],[false])
export DOCTEXT
AC_SUBST(DOCTEXT)


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
