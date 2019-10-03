AC_PREREQ(2.69)
AC_INIT([topoinfo],
        1.1,
        [wgropp@illinois.edu],
        [topoinfo],
        [http://www.cs.illinois.edu/~wgropp/])

AC_CONFIG_HEADER([topoconf.h])

AC_ARG_WITH([hwloc],AC_HELP_STRING([--with-hwloc=dir],[location of hwloc]))

# By default, try to find the best option.  Use -disable to skip one of
# these tests.
AC_ARG_ENABLE([mpinodename],
	AC_HELP_STRING([--enable-mpinodename],[Use MPI and nodename to find processes sharing the same node]),,enable_mpinodename=default)
AC_ARG_ENABLE([rca],AC_HELP_STRING([--enable-rca],[Use Cray RCA if available]),,enable_rca=yes)
AC_ARG_ENABLE([splittype],AC_HELP_STRING([--enable-splittype],[Use MPI_Comm_split_type to find nodes]),,enable_splittype=yes)
#
# Select the location and name for the base environment library
if test -n "$BASEENV_LIB" ; then
   TOPO_LIB=$BASEENV_LIB
else
   if test -z "$TOPO_LIBNAME" ; then
      TOPO_LIBNAME="topo"
   fi
   TOPO_LIB="lib$TOPO_LIBNAME.a"
fi
export TOPO_LIB
AC_SUBST(TOPO_LIB)

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

AC_PROG_RANLIB
# Find a working mkdir
AC_PROG_MKDIR_P

# Look for the MPI 3 function to find SMP nodes
if test "$enable_splittype" = "yes" ; then
    AC_CHECK_FUNCS([MPI_Comm_split_type])
fi
# Check whether to use MPI with gethostname to find processes sharing
# a hostname (and hence likely on the same node)
if test "$enable_mpinodename" = "default" ; then
   have_func_MPI_Init=no
   AC_CHECK_FUNC([MPI_Init],[have_func_MPI_Init=yes])
   if test "$have_func_MPI_Init" = "yes" ; then
       enable_mpinodename=yes
   fi
fi
if test "$enable_mpinodename" = "yes" ; then
   # Only enable this if either request or if the default makes this available
   # Note this is only used if no better approach is available
   AC_DEFINE([HAVE_MPI_WITH_NODENAME],1,[Set if MPI and nodename used to discover nodes])
fi

# Check for hwloc
if test -d "$with_hwloc/include" ; then
   CPPFLAGSsave=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS -I$with_hwloc/include"
   foundhwloc=no
   AC_CHECK_HEADERS([hwloc.h],[foundhwloc=yes],[CPPFLAGS=$CPPFLAGSsave])
   if test "$foundhwloc" = "yes" ; then
       LDFLAGSsave=$LDFLAGS
       LDFLAGS="$LDFLAGS -L$with_hwloc/lib"
       foundhwloc=no
       AC_CHECK_LIB([hwloc],[hwloc_topology_init],[foundhwloc=yes],[LDFLAGS=$LDFLAGSsave])
       # Ick.  May need to pass -dynamic to the linker; needed on Cray,
       # as the static libraries cause their own problems.

       # Note that some installations of hwloc apparently don't properly
       # set the rpath, so running (but not linking) can fail unless
       # LD_LIBRARY_PATH is set to include $with_hwloc/lib .
       # TODO: Add a check for that.  Action - reject hwloc in that case;
       # require the user to set LD_LIBRARY_PATH before invoking this configure.
   fi
fi

# nodetopo info
have_cray_rca=no
if test "$enable_rca" == yes ; then
    AC_CHECK_LIB([rca],[rca_get_meshcoord],[have_cray_rca=yes])
fi
if test "$have_cray_rca" = yes ; then
    AC_DEFINE([HAVE_CRAY_RCA],1,[Define if the Cray RCA library is available])
fi
have_bgq_mpix=no
have_mpix_h=no
have_new_mpix=no
AC_CHECK_HEADER([mpix.h],[have_mpix_h=yes],,[#include "mpi.h"])
if test $have_mpix_h = yes ; then
    AC_CHECK_FUNC([MPIX_Init_hw],[have_bgq_mpix=yes])
    # IBM apparently made an undocumented, gratuitous change in the
    # name of this routine (!!)
    AC_CHECK_FUNC([MPIX_Hardware],[have_bgq_mpix=yes;have_new_mpix=yes])
fi
if test $have_bgq_mpix = yes ; then
   AC_DEFINE([HAVE_BGQ_MPIX],1,[Define if BG/Q mpix routines are available])
   if test $have_new_mpix = yes ; then
      AC_DEFINE([HAVE_MPIX_HARDWARE],1,[Define if MPIX_Init_hw is now named MPIX_Hardware])
   fi
fi
# End of nodetopo info tests

# Check for node/cpu information
# Linux getcpu and the more limited sched_getcpu
AC_CHECK_HEADERS([utmpx.h])
AC_CHECK_FUNCS([getcpu sched_getcpu sched_getaffinity])

# Linux numa library
# Documented on our Cray, but not present.  Left here in hopes of finding
# a system that supports these routines
havenuma=no
AC_SEARCH_LIBS([numa_available],[numa],[havenuma=yes],)
if test "$havenuma" = "yes" ; then
    AC_CHECK_FUNCS([numa_node_to_cpus numa_node_of_cpu numa_get_run_node_mask],,[havenuma=no])
fi

MPICC=$CC
AC_SUBST(MPICC)

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

# In case the documentation targets are used, find doctext and attempt to
# find the source for the doctext LaTeX style files.  Use "false" if
# doctext is not found
AC_PATH_PROG(DOCTEXT,[doctext],[false])
export DOCTEXT
AC_SUBST(DOCTEXT)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
