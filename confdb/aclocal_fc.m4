dnl
dnl/*D
dnl PAC_PROG_FC_NAME_MANGLE - Determine how the Fortran compiler mangles
dnl names
dnl
dnl Synopsis:
dnl PAC_PROG_FC_NAME_MANGLE([action])
dnl
dnl Output Effect:
dnl If no action is specified, one of the following names is defined:
dnl.vb
dnl If fortran names are mapped:
dnl   lower -> lower                  FC_NAME_LOWER
dnl   lower -> lower_                 FC_NAME_LOWER_USCORE
dnl   lower -> UPPER                  FC_NAME_UPPER
dnl   lower_lower -> lower__          FC_NAME_LOWER_2USCORE
dnl   mixed -> mixed                  FC_NAME_MIXED
dnl   mixed -> mixed_                 FC_NAME_MIXED_USCORE
dnl   mixed -> UPPER@STACK_SIZE       FC_NAME_UPPER_STDCALL
dnl.ve
dnl If an action is specified, it is executed instead.
dnl
dnl Notes:
dnl We assume that if lower -> lower (any underscore), upper -> upper with the
dnl same underscore behavior.  Previous versions did this by 
dnl compiling a Fortran program and running strings -a over it.  Depending on 
dnl strings is a bad idea, so instead we try compiling and linking with a 
dnl C program, since that is why we are doing this anyway.  A similar approach
dnl is used by FFTW, though without some of the cases we check (specifically, 
dnl mixed name mangling).  STD_CALL not only specifies a particular name
dnl mangling convention (adding the size of the calling stack into the function
dnl name, but also the stack management convention (callee cleans the stack,
dnl and arguments are pushed onto the stack from right to left)
dnl
dnl One additional problem is that some Fortran implementations include 
dnl references to the runtime (like pgf90_compiled for the pgf90 compiler
dnl used as the "Fortran 77" compiler).  This is not yet solved.
dnl
dnl D*/
dnl
AC_DEFUN([PAC_PROG_FC_NAME_MANGLE],[
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_CACHE_CHECK([for Fortran name mangling],
pac_cv_prog_fc_name_mangle,[
# Initialize pac_found to indicate if name mangling scheme has been found
pac_found=no
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([
        subroutine MY_name( ii )
        return
        end
    ])
],[
    PAC_RUNLOG([mv conftest.$OBJEXT fcconftest.$OBJEXT])
    saved_LIBS="$LIBS"
    dnl  FLIBS is set by AC_FC_LIBRARY_LDFLAGS
    LIBS="fcconftest.$OBJEXT $FLIBS $LIBS"
    AC_LANG_PUSH([C])
    for call in "" __stdcall ; do
        for sym in my_name_ my_name__ my_name MY_NAME MY_name MY_name_ NONE ; do
            AC_LINK_IFELSE([
                AC_LANG_PROGRAM([extern void ${call} ${sym}(int);],[${sym}(0);])
            ],[
                pac_found=yes
                break
            ])
        done
        test "$pac_found" = "yes" && break
    done
    AC_LANG_POP([C])
    LIBS="$saved_LIBS"
    rm -f fcconftest.$OBJEXT
])
AC_LANG_POP([Fortran])
dnl
# If we got to here and pac_cv_prog_fc_name_mangle is still NOT definable,
# it may be that the programs have to be linked with the Fortran compiler,
# not the C compiler.  Try reversing the language used for the test
if test  "$pac_found" != "yes" ; then
    AC_LANG_PUSH([C])
    for call in "" __stdcall ; do
        for sym in my_name_ my_name__ my_name MY_NAME MY_name MY_name_ NONE ; do
            AC_COMPILE_IFELSE([
                AC_LANG_SOURCE([void ${call} ${sym}(int a) {}])
            ],[
                PAC_RUNLOG([mv conftest.$OBJEXT cconftest.$OBJEXT])
                saved_LIBS="$LIBS"
                LIBS="cconftest.$OBJEXT $LIBS"
                AC_LANG_PUSH([Fortran])
                AC_LINK_IFELSE([
                    AC_LANG_PROGRAM([],[      call my_name(0)])
                ],[
                    pac_found=yes
                ])
                AC_LANG_POP([Fortran])
                LIBS="$saved_LIBS"
                rm -f cconftest.$OBJEXT
                test "$pac_found" = "yes" && break
            ])
        done
        test "$pac_found" = "yes" && break
    done
    AC_LANG_POP([C])
fi
if test "$pac_found" = "yes" ; then
    case ${sym} in
        my_name_)
            pac_cv_prog_fc_name_mangle="lower uscore" ;;
        my_name__)
            pac_cv_prog_fc_name_mangle="lower 2uscore" ;;
        my_name)
            pac_cv_prog_fc_name_mangle="lower" ;;
        MY_NAME)
            pac_cv_prog_fc_name_mangle="upper" ;;
        MY_name)
            pac_cv_prog_fc_name_mangle="mixed" ;;
        MY_name_)
            pac_cv_prog_fc_name_mangle="mixed uscore" ;;
        *)
            pac_cv_prog_fc_name_mangle=""
            pac_found=no;
            ;;
    esac
    if test "X$pac_cv_prog_fc_name_mangle" != "X" ; then
        if test "$call" = "__stdcall" ; then
            pac_cv_prog_fc_name_mangle="$pac_cv_prog_fc_name_mangle stdcall"
        fi
    fi
fi
])
dnl Endof ac_cache_check
case $pac_cv_prog_fc_name_mangle in
    *stdcall)
        FC_STDCALL="__stdcall" ;;
    *)
        FC_STDCALL="" ;;
esac
# Get the standard call definition
# FIXME: This should use FC_STDCALL, not STDCALL (non-conforming name)
FC_STDCALL="$call"
AC_DEFINE_UNQUOTED(STDCALL,[$FC_STDCALL],[Define calling convention])

# new_name="`echo $name | tr ' ' '_' | tr [a-z] [A-Z]`"
# We could have done the character conversion with 'tr'
# which may not be portable, e.g. solaris's /usr/ucb/bin/tr.
# So use a conservative approach.

# Replace blank with underscore
name_scheme="`echo $pac_cv_prog_fc_name_mangle | sed 's% %_%g'`"
# Turn lowercase into uppercase.
name_scheme="`echo $name_scheme | sed -e 'y%abcdefghijklmnopqrstuvwxyz%ABCDEFGHIJKLMNOPQRSTUVWXYZ%'`"
FC_NAME_MANGLE="FC_NAME_${name_scheme}"
AC_DEFINE_UNQUOTED([$FC_NAME_MANGLE])
AC_SUBST(FC_NAME_MANGLE)
if test "X$pac_cv_prog_fc_name_mangle" = "X" ; then
    AC_MSG_WARN([Unknown Fortran naming scheme])
fi
dnl
dnl Define the macros that is needed by AC_DEFINE_UNQUOTED([$FC_NAME_MANGLE])
AH_TEMPLATE([FC_NAME_LOWER],
    [Fortran names are lowercase with no trailing underscore])
AH_TEMPLATE([FC_NAME_LOWER_USCORE],
    [Fortran names are lowercase with one trailing underscore])
AH_TEMPLATE([FC_NAME_LOWER_2USCORE],
    [Fortran names are lowercase with two trailing underscores])
AH_TEMPLATE([FC_NAME_MIXED],
    [Fortran names preserve the original case])
AH_TEMPLATE([FC_NAME_MIXED_USCORE],
    [Fortran names preserve the original case with one trailing underscore])
AH_TEMPLATE([FC_NAME_UPPER],
    [Fortran names are uppercase])
AH_TEMPLATE([FC_NAME_LOWER_STDCALL],
    [Fortran names are lowercase with no trailing underscore in stdcall])
AH_TEMPLATE([FC_NAME_LOWER_USCORE_STDCALL],
    [Fortran names are lowercase with one trailing underscore in stdcall])
AH_TEMPLATE([FC_NAME_LOWER_2USCORE_STDCALL],
    [Fortran names are lowercase with two trailing underscores in stdcall])
AH_TEMPLATE([FC_NAME_MIXED_STDCALL],
    [Fortran names preserve the original case in stdcall])
AH_TEMPLATE([FC_NAME_MIXED_USCORE_STDCALL],
    [Fortran names preserve the original case with one trailing underscore in stdcall])
AH_TEMPLATE([FC_NAME_UPPER_STDCALL],
    [Fortran names are uppercase in stdcall])
# Define the function macros used by the GNU wrapper code
AH_TEMPLATE([FC_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran.])dnl
AH_TEMPLATE([FC_FUNC_],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), containing underscores, for linking with Fortran.])dnl
dnl
    case ${pac_cv_prog_fc_name_mangle} in
       "lower uscore")
       AC_DEFINE([FC_FUNC(name,NAME)],[name [##] _])
       AC_DEFINE([FC_FUNC_(name,NAME)],[name [##] _])
       ;;
       "lower 2uscore")
       AC_DEFINE([FC_FUNC(name,NAME)],[name [##] _])
       AC_DEFINE([FC_FUNC_(name,NAME)],[name [##] __])
       ;;
       "lower")
       AC_DEFINE([FC_FUNC(name,NAME)],[name])
       AC_DEFINE([FC_FUNC_(name,NAME)],[name])
       ;;
       "upper")
       AC_DEFINE([FC_FUNC(name,NAME)],[NAME])
       AC_DEFINE([FC_FUNC_(name,NAME)],[NAME])
       ;;
       "mixed")
       AC_MSG_WARN([Mixed Fortran wrapper name not supported by FC_FUNC])
       ;;
       "mixed uscore")
       AC_MSG_WARN([Mixed Fortran wrapper name not supported by FC_FUNC])
       ;;
       *)
       AC_MSG_WARN([Unknown Fortran wrapper name not supported by FC_FUNC])
       ;;
    esac
])
