AC_PREREQ(2.69)
AC_INIT([haloperf],
	1.0,
	[wgropp@illinois.edu],
	[haloperf],
	[http://wgropp.cs.illinois.edu])

AC_CONFIG_HEADERS(haloperf.h)

AC_ARG_ENABLE(echo,
[--enable-echo    - Turn on echoing from within the configure script],
set -x)

# Set the C and Fortran FLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for this package
: ${CFLAGS=""}
: ${FCFLAGS=""}

PAC_PROG_MPICC
if test -n "$MPICC" ; then CC=$MPICC ; fi

AC_ARG_WITH(topodesc,
    AC_HELP_STRING([--with-topodesc=dir],[Use the topodesc library to determine process topology.  Files in dir/include and dir/lib]),,)
AC_ARG_WITH(topodesclib,
    AC_HELP_STRING([--with-topodesclib=dir],[Give the directory for libtopodesc]),,)
AC_ARG_WITH(topodescinc,
    AC_HELP_STRING([--with-topodescinc=dir],[Give the directory for libtopodesc include file]),,)

# At this point, setup the language system.  From autoconf 2.57, this requires
# calling the PROG_CC etc
AC_PROG_CC
PAC_MPI_CCWORKS
# Use COPT for optimization-only flags for C
AC_SUBST([COPT])

AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT

AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h unistd.h])


AC_CHECK_FUNCS(MPI_Init_thread)

if test -n "$BASEENV_LIB" ; then
   AC_SUBST([BASEENV_LIB])
   # Location of topo code if not yet installed
   :
fi
if test -s ../topo/topoinfo.h ; then
   CPPFLAGS="$CPPFLAGS -I../topo"
   # FIXME!!
   LIBS="$LIBS -lbaseenv"
   LDFLAGS="$LDFLAGS -L.."
fi

# The following is old, to be deleted
if test "$with_topodesc" ; then
   if test -z "$with_topodesclib" ; then 
       with_topodesclib="$with_topodesc/lib"
   fi
   if test -z "$with_topodescinc" ; then 
       with_topodescinc="$with_topodesc/include"
   fi
fi
if test -n "$with_topodescinc" ; then 
    if test -s "$with_topodescinc/findcliques.h" ; then
        AC_DEFINE(HAVE_FINDCLIQUES_H,1,[Define if findcliques.h from topodesc available])
    else
        AC_MSG_ERROR([Could not find $with_topodesc/include/findcliques.h])
    fi
fi
if test -n "$with_topodesclib" ; then 
   LDFLAGS="$LDFLAGS -L$with_topodesclib"
   LIBS="$LIBS -ltopodisc"
   CPPFLAGS="$CPPFLAGS -I$with_topodescinc"
fi

have_cray_rca=no
AC_CHECK_LIB(rca,rca_get_meshcoord,have_cray_rca=yes)
if test "$have_cray_rca" = yes ; then
   AC_DEFINE(HAVE_CRAY_RCA,1,[Define if the Cray RCA library is available])
fi
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------

AC_OUTPUT(Makefile)
