#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([libtopodisc], [1.0], [wgropp@illinois.edu])

# Set the C and Fortran FLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for this package
: ${CFLAGS=""}
: ${FCFLAGS=""}

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Find the compilers, if they aren't set with the environment variables
if test -z "$MPICC" ; then
    AC_MSG_CHECKING([for how to compile MPI program in C])
    if test -x $MPI_DIR/mpicc ; then
        MPICC=$MPI_DIR/mpicc
    elif test -x $MPI_DIR/bin/mpicc ; then
        MPICC=$MPI_DIR/bin/mpicc
    elif test -x $MPI_DIR/mpcc ; then
        MPICC=$MPIDIR/mpcc
    elif test -x $MPI_DIR/bin/mpcc ; then
        MPICC=$MPIDIR/bin/mpcc
    else
        saveIFS=$IFS
        IFS=":"
	for dir in $PATH ; do
	    if test -x $dir/mpicc ; then
                MPICC=$dir/mpicc
		break
            elif test -x $dir/mpcc ; then
                MPICC=$dir/mpcc
		break
            fi
	done
	IFS=$saveIFS
    fi
    if test -n "$MPICC" ; then
        AC_MSG_RESULT($MPICC)
	CC=$MPICC
    else
        AC_MSG_RESULT([cannot determine!])
    fi
fi

if test -z "$MPICXX" ; then
    AC_MSG_CHECKING([for how to compile MPI program in C++])
    if test -x $MPI_DIR/mpicxx ; then
        MPICXX=$MPI_DIR/mpicxx
    elif test -x $MPI_DIR/bin/mpicxx ; then
        MPICXX=$MPI_DIR/bin/mpicxx
    elif test -x $MPI_DIR/mpCC ; then
        MPICXX=$MPIDIR/mpCC
    elif test -x $MPI_DIR/bin/mpCC ; then
        MPICXX=$MPIDIR/bin/mpCC
    else
        saveIFS=$IFS
        IFS=":"
	for dir in $PATH ; do
	    if test -x $dir/mpicxx ; then
                MPICXX=$dir/mpicxx
		break
            elif test -x $dir/mpCC ; then
                MPICXX=$dir/mpCC
		break
            fi
	done
	IFS=$saveIFS
    fi
    if test -n "$MPICXX" ; then
        AC_MSG_RESULT($MPICXX)
	CXX=$MPICXX
    else
        AC_MSG_RESULT([cannot determine!])
    fi
fi

AM_PROG_AR
LT_INIT
AC_CONFIG_SRCDIR([findcliques.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PATH_PROG([DOCTEXT],[doctext],[false])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gethostname])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
