AC_PREREQ(2.69)
# Additional args to init: mail address, name, website
AC_INIT([tsc],1.1)
AC_CONFIG_HEADER([tscdef.h])
AC_ARG_ENABLE([basetype],
   AC_HELP_STRING([--enable-basetype=type],
   [Set the data type for the vector computions.  This is float by default, but can be set to double]),,enable_basetype=float)
AC_ARG_ENABLE([ntimes],
   AC_HELP_STRING([--enable-ntimes],
   [Set the number of times tests are run for timing at runtime rather than fixed at compile time]),,enable_ntimes=no)

# Set the CFLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for this package
: ${CFLAGS=""}
AC_PROG_CC
AC_PROG_CC_C99
AC_C_RESTRICT
AC_C_INLINE
# Attempt to find the compiler vendor.  Currently just for informational
# purposes
AX_COMPILER_VENDOR
AC_CHECK_HEADERS([malloc.h xmmintrin.h builtins.h])
AC_CHECK_FUNCS([memalign posix_memalign])
if test "$enable_basetype" = "float" ; then
    AC_SEARCH_LIBS([sinf],m)
elif test "$enable_basetype" = "double" ; then
    AC_SEARCH_LIBS([sin],m)
fi
# posix_memalign may need this:
# _POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600

# Find a working mkdir
AC_PROG_MKDIR_P

# In case the documentation targets are used, find doctext and attempt to
# find the source for the doctext LaTeX style files.  Use "false" if
# doctext is not found
AC_PATH_PROG(DOCTEXT,[doctext],[false])
export DOCTEXT
AC_SUBST(DOCTEXT)

AC_DEFINE_UNQUOTED([BASETYPE],[$enable_basetype],[Base type for vector operations])
BASETYPE_VAL=0
# The values here must match BASETYPE_FLOAT and BASETYPE_DOUBLE in tsc.c
if test "$enable_basetype" = "float" ; then
   BASETYPE_VAL=1
elif test "$enable_basetype" = "double" ; then
   BASETYPE_VAL=2
else
    AC_MSG_ERROR([basetype must be either float or double])
fi
AC_DEFINE_UNQUOTED([BASETYPE_VAL],[$BASETYPE_VAL],[Value we can use to test in C preprocessor])
if test "$enable_ntimes" = "yes" ; then
   AC_DEFINE([RUNTIME_NTIMES],1,[Define if ntimes is set a runtime rather than compile time])
fi

# Normally, this configure should be invoked by the top-level baseenv
# builddir.  If not, this code properly defines the master_top_builddir
# and master_top_srcdir to build this as a standalone package
if test "X$master_top_builddir" = "X" ; then
   if test "X$ac_top_builddir" = "X" ; then
       master_top_builddir=$ac_pwd
   else
       master_top_builddir=$ac_top_builddir
   fi
fi
if test "X$master_top_srcdir" = "X" ; then
   if test "X$ac_top_srcdir" = "X" ; then
       master_top_srcdir=$srcdir
   else
       master_top_srcdir=$ac_top_srcdir
   fi
fi

AC_SUBST([master_top_builddir])
AC_SUBST([master_top_srcdir])

# Check for the location of Makefile.doc and Makefile.vec
if test -s $srcdir/Makefile.vec ; then
    BASEENVMAIN=$srcdir
elif test -s $master_top_srcdir/Makefile.vec ; then
    BASEENVMAIN=$master_top_srcdir
elif test -s $srcdir/../Makefile.vec ; then
    BASEENVMAIN=$srcdir/..
else
    AC_MSG_ERROR([Unable to find Makefile.vec])
fi
AC_SUBST([BASEENVMAIN])

AC_OUTPUT([Makefile])
