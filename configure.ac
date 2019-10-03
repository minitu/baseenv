AC_PREREQ(2.63)
AC_INIT([baseenv],
        [0.94.1],
        [wgropp@illinois.edu],
        [baseenv],
        [http://wgropp.cs.illinois.edu/])

# Location of some extra macros
AC_CONFIG_MACRO_DIR([confdb])
# Uncomment the following to move the various autoconf required scripts
# into the confdb subdirectory
dnl AC_CONFIG_AUX_DIR([confdb])

AC_CONFIG_HEADER([baseenv.h])
AC_ARG_ENABLE([echo],
    [AS_HELP_STRING([--enable-echo],
    [Turn on verbose shell echo when running configure])],
    [set -x],
    [set +x])

dnl strict is only for GNU compilers and compilers that accept GNU warning
dnl options.  Note that we have to compensate for the idiotically named -Wall
dnl that is in no way even close to turning on all of the warnings.
AC_ARG_ENABLE([strict],
    [AS_HELP_STRING([--enable-strict],
    [])],
    [CFLAGS="$CFLAGS -Wall -Wextra -Wno-missing-field-initializers -Wstrict-prototypes -Wmissing-prototypes -Wno-unused-parameter -Wno-unused-label -Wshadow -Wmissing-declarations -Wno-long-long -Wundef -Wno-endif-labels -Wpointer-arith -Wbad-function-cast -Wwrite-strings -Wno-sign-compare -Wold-style-definition -Wno-multichar -Wno-deprecated-declarations -Wpacked -Wnested-externs -Winvalid-pch -Wno-pointer-sign -Wno-format-zero-length -Wtype-limits -Werror-implicit-function-declaration"],)
export CFLAGS
#
# Specify the location of MPI compilation scripts, or MPI install
# directory.  Not necessary in most cases (e.g., scripts in the
# search path)
AC_ARG_WITH([mpi],
   [AS_HELP_STRING([--with-mpi=dir],
                   [location of MPI compilation and run scripts])],
   [],
   [with_mpi=/usr/local/bin])
#
# Optional packages.  An location can be provided (directory)
# Not all of the packages may ship with the baseenv package
AC_ARG_WITH([fpmpi2],
   [AS_HELP_STRING([--with-fpmpi2=dir],
                   [location of FPMPI2 libraries])],
   [if test "x$with_fpmpi2" = "x" ; then with_fpmpi2=fpmpi2 ; fi],
   [with_fpmpi2=no])
AC_ARG_WITH([topo],
   [AS_HELP_STRING([--with-topo],
                   [Include the simple process topology routines])],
   [if test "x$with_topo" = "x" ; then with_topo=topo ; fi],
   [with_topo=topo])
AC_ARG_WITH([libtopodisc],
   [AS_HELP_STRING([--with-libtopodisc],
                   [Include the libtopodisc process topology routines])],
   [if test "x$with_libtopodisc" = "x" ; then with_libtopodisc=libtopodisc ; fi],
   [with_libtopodisc=no])
AC_ARG_WITH([hwloc],
   [AS_HELP_STRING([--with-hwloc],
                   [Include the hwloc hardware locality package])],
   [if test "x$with_hwloc" = "x" ; then with_hwloc=hwloc ; fi ],
   [with_hwloc=no])
AC_ARG_WITH([perftest],
   [AS_HELP_STRING([--with-perftest],
                   [Include the performance test package])],
   [if test "x$with_perftest" = "x" ; then with_perftest=perftest; fi ],
   [with_perftest=perftest])
AC_ARG_WITH([vectortest],
   [AS_HELP_STRING([--with-vectortest],
                   [Include the vector performance test package])],
   [if test "x$with_vectortest" = "x" ; then with_vectortest=TSVC; fi ],
   [with_vectortest=TSVC])
AC_ARG_WITH([haloperf],
   [AS_HELP_STRING([--with-haloperf],
                   [Include the performance test haloperf])],
   [if test "x$with_haloperf" = "x" ; then with_haloperf=haloperf; fi ],
   [with_haloperf=haloperf])

#
# Select the location and name for the base environment library
if test -z "$BASEENV_LIBNAME" ; then
   BASEENV_LIBNAME="baseenv"
   BASEENVF_LIBNAME="baseenvf"
fi
if test -z "$BASEENV_LIBDIR" ; then
   BASEENV_LIBDIR=`pwd`
fi
BASEENV_LIB="$BASEENV_LIBDIR/lib$BASEENV_LIBNAME.a"
BASEENVF_LIB="$BASEENV_LIBDIR/lib$BASEENVF_LIBNAME.a"
export BASEENV_LIB
export BASEENV_LIBNAME
export BASEENVF_LIBNAME
# Export that we are building the SEQ tools
BUILDING_SEQ=yes
export BUILDING_SEQ
AC_ARG_VAR([BUILDING_SEQ],[Set to yes if building the seq package])
AC_ARG_VAR([BASEENV_LIB],[Set to the library path for the basic environment])
AC_ARG_VAR([BASEENV_LIBNAME],[Set to the library name for the basic environment])
AC_ARG_VAR([BASEENVF_LIB],[Set to the library path for the Fortran basic environment])
AC_ARG_VAR([BASEENVF_LIBNAME],[Set to the library name for the Fortran basic environment])

# Query: Do we need to do the same for mputil?

# Set a default for the scratch dir
if test "X$MYSCRATCHDIR" = "X" ; then
    MYSCRATCHDIR="."
fi
AC_SUBST([MYSCRATCHDIR])


# ----------------------------------------------------------------------------
# This test is complicated by the fact that top_srcdir is not set until
# the very end of configure.  Instead, we get it ourselves
if test -z "$top_srcdir" ; then
   use_top_srcdir=$srcdir
else
   use_top_srcdir=$top_srcdir
fi
if test -z "$master_top_srcdir" ; then
    # This needs to be an absolute pathname
    case "$use_top_srcdir" in
    /*) ;;
    *)
        use_top_srcdir=`(cd $use_top_srcdir && pwd)`
	;;
    esac
    master_top_srcdir=$use_top_srcdir
fi
# Get the directory that we're running in...
if test -z "$master_top_builddir" ; then
   master_top_builddir="`pwd`"
fi
export master_top_builddir
export master_top_srcdir
AC_SUBST([master_top_builddir])
AC_SUBST([master_top_srcdir])
# ----------------------------------------------------------------------------
# Check for a file containing the git hash.  If none present, use unknown
# ----------------------------------------------------------------------------
githashvalue="Unknown"
if test -s $master_top_srcdir/githash ; then
    githashvalue=`cat $master_top_srcdir/githash`
fi
AC_DEFINE_UNQUOTED([BASEENV_GITHASH],"$githashvalue",[Git Hash value for baseenv for this release])
# ----------------------------------------------------------------------------
#
# Try to find the compiler
if test "x$with_mpi" != "x" ; then
    if test -d "$with_mpi" -a -x "$with_mpi/mpicc" ; then
        MPI_DIR=$with_mpi
    elif test -d "$with_mpi/bin" ; then
        MPI_DIR="$with_mpi/bin"
    fi
fi

# Set the C and Fortran FLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for this package
: ${CFLAGS=""}
: ${FCFLAGS=""}
# Optimization should be in COPT and FCOPT
AC_SUBST([COPT])
AC_SUBST([FCOPT])

PAC_PROG_MPICC

# Now, make sure that we have something that we can use
if test -n "$MPICC" ; then
    # Note that this may make CC a cross compiler
    CC=$MPICC
fi

PAC_PROG_MPIFORT
if test -n "$MPIFORT" ; then
    FC=$MPIFORT
fi

# Basic properties of the C compiler.  Small components that don't
# have separate builds can use these values
AC_PROG_CC
AC_C_VOLATILE
AC_C_RESTRICT
AC_C_INLINE

# Find how to build libraries
AC_PROG_RANLIB
# Find a working mkdir
AC_PROG_MKDIR_P
# Installation programs
AC_PROG_INSTALL
#
# Find the MPI location
PAC_MPI_CCWORKS

# Check for mpi.h (needed for some of the other tests)
AC_CHECK_HEADERS([mpi.h])

# Check for node/cpu information
# Linux getcpu and the more limited sched_getcpu.  Used for supporting
# node-awareness of socket and core.
AC_CHECK_HEADERS([utmpx.h])
AC_CHECK_FUNCS([getcpu sched_getcpu sched_getaffinity])
# sysctlbyname is a convenient interface to some configuration info
# on BSD. Linux, sadly, has nothing equivalent.
AC_CHECK_FUNCS([sysctlbyname])
# Look for hwloc, which provides extensive information on the node,
# particularly for locality (on systems that do not actively hide the
# information about processes and threads, such as OSX).
AC_CHECK_HEADERS([hwloc.h])
LIBHWLOC=false
HWLOCLIB=
dnl We use CHECK_LIB because the hwloc should be included only where needed
dnl using SEARCH_LIBS unconditionally adds -lhwloc to LIBS
# It only makes sense to use for the hwloc library if we found the hwloc
# header.
if test "$ac_cv_header_hwloc_h" = "yes" ; then
    AC_MSG_CHECKING([for hwloc routines in default libraries])
    AC_LINK_IFELSE([AC_LANG_CALL([], [hwloc_topology_init])],[LIBHWLOC=true],[])
    if test "$LIBHWLOC" = "false" ; then
        AC_MSG_RESULT([no])
        # Now look in -lhwloc
        AC_CHECK_LIB([hwloc],[hwloc_topology_init],[LIBHWLOC=true;HWLOCLIB=-lhwloc],[],[])
    else
        AC_MSG_RESULT([yes])
    fi
fi

if test -r /proc/cpuinfo ; then
    AC_DEFINE([HAVE_PROC_CPUINFO],1,[Define if /proc/cpuinfo is readable])
fi

# Check for the optional packages.  Make sure that that are present if selected
if test "$with_fpmpi2" != no ; then
    FPMPI2_DIR=$with_fpmpi2
    if test "$FPMPI2_DIR" = fpmpi2 ; then
        if test -d fpmpi2 ; then
            AC_CONFIG_SUBDIRS([fpmpi2])
        fi
    elif test -d "$FPMPI2_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $FPMPI2_DIR"
    else
        AC_MSG_ERROR([Directory $with_fpmpi2 for fpmpi2 does not exist])
    fi
fi

if test "$with_topo" != no ; then
    TOPO_DIR=$with_topo
    if test "$TOPO_DIR" = topo ; then
        if test -d topo ; then
            AC_CONFIG_SUBDIRS([topo])
	    TOPOINC="-I$master_top_srcdir/topo"
        fi
    elif test -d "$TOPO_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $TOPO_DIR"
	TOPOINC="-I$TOPO_DIR"
    else
        AC_MSG_ERROR([Directory $with_topo for topo does not exist])
    fi
    AC_DEFINE([HAVE_TOPOINIT],1,[Define when simple topology routines available])
fi
AC_SUBST([TOPOINC])

if test "$with_perftest" != no ; then
    PERFTEST_DIR=$with_perftest
    if test "$PERFTEST_DIR" = perftest ; then
        if test -d perftest ; then
            AC_CONFIG_SUBDIRS([perftest])
        fi
    elif test -d "$PERFTEST_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $PERFTEST_DIR"
    else
        AC_MSG_ERROR([Directory $with_perftest for perftest does not exist])
    fi
fi

if test "$with_vectortest" != no ; then
    VECTORTEST_DIR=$with_vectortest
    if test "$VECTORTEST_DIR" = TSVC ; then
        if test -d TSVC ; then
            AC_CONFIG_SUBDIRS([TSVC])
        fi
    elif test -d "$VECTORTEST_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $VECTORTEST_DIR"
    else
        AC_MSG_ERROR([Directory $with_vectortest for vectortest does not exist])
    fi
fi

if test "$with_haloperf" != no ; then
    HALOPERF_DIR=$with_haloperf
    if test "$HALOPERF_DIR" = haloperf ; then
        if test -d haloperf ; then
            AC_CONFIG_SUBDIRS([haloperf])
        fi
    elif test -d "$HALOPERF_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $HALOPERF_DIR"
    else
        AC_MSG_ERROR([Directory $with_haloperf for haloperf does not exist])
    fi
fi

if test "$with_libtopodisc" != no ; then
    LIBTOPODISC_DIR=$with_libtopodisc
    if test "$LIBTOPODISC_DIR" = libtopodisc ; then
        if test -d libtopodisc ; then
            AC_CONFIG_SUBDIRS([libtopodisc])
        fi
    elif test -d "$LIBTOPODISC_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $LIBTOPODISC_DIR"
    else
        AC_MSG_ERROR([Directory $with_libtopodisc for libtopodisc does not exist])
    fi
fi

if test "$with_hwloc" != no ; then
    HWLOC_DIR=$with_hwloc
    if test "$HWLOC_DIR" = hwloc ; then
        if test -d hwloc ; then
            AC_CONFIG_SUBDIRS([hwloc])
        fi
    elif test -d "$HWLOC_DIR" ; then
        # This is much harder.  Keep track of these; we'll handle them
	# later
	externalDirs="$externalDirs $HWLOC_DIR"
    else
        AC_MSG_ERROR([Directory $with_hwloc for hwloc does not exist])
    fi
fi

# Codes that use nodecart whold add HWLOCLIB to their lib list
# Most codes do not need hwloc, so we don't include it by default
if test "$LIBHWLOC" = "true" ; then
   AC_DEFINE([HAVE_HWLOC_LIB],[1],[Define if hwloc library is available])
fi
AC_SUBST([HWLOCLIB])

MPICC=$CC
# Handle the external packages.  We know that the directory exists
# We can't use the configure mechanism for this, so note that config.status
# won't know about these (we could add to the AC_OUTPUT actions, but
# that is more than I'm willing to do here
externsubdirs=""
AC_SUBST([externsubdirs])
for dir in $externalDirs ; do
    if test -x $dir/configure.gnu ; then
        externconfigure="./configure.gnu"
    elif test -x $dir/configure ; then
        externconfigure="./configure"
    else
        AC_MSG_ERROR([Cannot find configure for $dir])
    fi
    # Process args so that we can safely pass them to the configure
    # scripts.  A simplified version of what autoconf does
    pac_config_extern_args="--disable-option-checking "
    eval "set x $ac_configure_args"
    shift
    for ac_arg ; do
        case $ac_arg in
	MPICC=*) ;;
	CC=*) ;;
	*\'*) ac_arg=`$as_echo "$ac_arg" | sed "s/'/'\\\\\\\\''/g"` ;;
        esac
        as_fn_append pac_config_extern_args " '$ac_arg'"
    done
    as_fn_append pac_config_extern_args " 'MPICC=$MPICC'"
    as_fn_append pac_config_extern_args " 'CC=$CC'"
    ac_popdir=`pwd`
    AC_MSG_NOTICE([=== Configuring external package in $dir])
    cd "$dir"
    AC_MSG_NOTICE([Running $SHELL $externconfigure $pac_config_extern_args])
    eval "\$SHELL \"\$externconfigure\" $pac_config_extern_args \
        --cache-file=\"/dev/null\""
    # If successful, add to externsubdirs
    if test $? = 0 ; then externsubdirs="$externsubdirs $dir" ; fi
    cd $ac_popdir
done

# Check for the raw Fortran compiler
AC_PROG_FC
# There are some additional tests to consider for better support of
# Fortran.  This is complicated by the fact that we want to support both
# mpif.h (Fortran 77) and use mpi (Fortran 9x) codes, though note that many
# Fortran MPI applications and especially benchmarks use mpif.h
dnl AC_FC_SRCEXT(f90)
dnl AC_FC_FREEFORM
dnl AC_FC_IMPLICIT_NONE
dnl AC_FC_MODULE_FLAG
dnl AC_FC_MODULE_EXTENSION
dnl AC_FC_MODULE_OUTPUT_FLAG
# Enable finding the wrappers enable Fortran routines to be
# written in C
# The default autoconf implementation of this does not work correctly with
# the Cray compiler.  Neither Cray nor Gnu is willing to fix their respective
# codes (Cray's bizarre approach is to patch autoconf on Cray systems, forcing
# you to rebuild configure - and that doesn't even work.  GNU uses an
# implementation of this macro that is broken on Cray, and wants Cray to
# change how their Fortran compiler works, even though it works fine when
# used the way most Fortran programmers use it.
dnl AC_FC_WRAPPERS
# Instead, we use this macro
PAC_PROG_FC_NAME_MANGLE

# A special test - on some systems, the Fortran compiler can build MPI programs
# Test for that here
if test -z "$MPIFORT" ; then
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([whether $FC can link MPI programs])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],[
      integer ierror
      call mpi_init(ierror)])],
    [ftnMPIWorks=yes],[ftnMPIWorks=no])
    AC_MSG_RESULT([$ftnMPIWorks])
    if test "$ftnMPIWorks" = "yes" ; then
        MPIFORT=$FC
    fi
    AC_LANG_POP([Fortran])
fi
# In case the documentation targets are used, find doctext and attempt to
# find the source for the doctext LaTeX style files.  Use "false" if
# doctext is not found
AC_PATH_PROG(DOCTEXT,[doctext],[false])
export DOCTEXT
AC_SUBST(DOCTEXT)

# Subdirectory-specific tests
# These are used with the packages that do not have their own configure.
# If these become complex, consider moving them into a configure for that
# directory
# IO tests
AC_CHECK_HEADERS([lustre/lustreapi.h])
AC_CHECK_LIB([lustreapi],[llapi_file_open],[LLAPILIB=-llustreapi],[LLAPILIB=],[])
AC_SUBST([LLAPILIB])
if test "X$LLAPILIB" != "X" ; then
    AC_DEFINE([HAVE_LLAPI_FILE_OPEN],1,[Define if llapi_file_open in lustreapi])
fi

AC_SUBST(MPI_DIR)
AC_SUBST(FPMPI2_DIR)
AC_SUBST(MPICC)
AC_CONFIG_FILES([Makefile
		 seq/Makefile
		 util/Makefile
		 mputil/Makefile
                 io/Makefile
		 io/tests/Makefile
	       	 mpit/Makefile
		 nodeperf/Makefile
		 mpi-patterns/Makefile
		 mpi-patterns/stencil3/Makefile
                 nodecart/Makefile
                 nodecart/tests/Makefile
		 gpuvec/Makefile
])
AC_OUTPUT
