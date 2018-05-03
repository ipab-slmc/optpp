#                        -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPTPP], [2.4], [http://software.sandia.gov/bugzilla])

AC_CONFIG_SRCDIR([src/globals.C])
AC_CONFIG_HEADER([include/OPT++_config.h])
AC_CONFIG_AUX_DIR([config])

# Ask Bill about the next line

#AC_CANONICAL_BUILD

AC_CANONICAL_HOST
AC_DISABLE_SHARED

AM_INIT_AUTOMAKE([1.9.2 foreign])

# Checks for programs.

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77

# Talk to Bill about the following block

#AC_PROG_CC([cc xlc gcc])
#if test X${build_os} = Xcygwin ; then
#  cpp_compilers="g++ gcc"
#else
#  cpp_compilers="CC cxx xlC c++ g++"
#fi
#AC_PROG_CXX($cpp_compilers)

AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Check whether we're building as part of DAKOTA.

AC_MSG_CHECKING([whether dakota_config.h.in exists])
if test -f "$srcdir/../Dakota/dakota_config.h.in" -o \
        -f "$srcdir/../../dakota_config.h.in"; then
  AC_MSG_RESULT([yes])
  AC_DEFINE([DAKOTA_OPTPP],[1],[Identifier for non-standalone builds.])
else
  AC_MSG_RESULT([no])
fi

AC_LANG([C++])

# Checks for libraries.

# Checks for header files.

AC_CXX_HAVE_STD
AC_LANG_PUSH([C])
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([fcntl.h float.h stddef.h stdlib.h string.h \	
		  sys/param.h sys/types.h sys/time.h sys/times.h \
		  sys/resource.h unistd.h values.h])
AC_LANG_POP([C])

# Checks for typedefs, structures, and compiler characteristics.

AC_LANG_PUSH([C])
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM
AC_LANG_POP([C])
AC_CXX_COMPILER_VENDOR([CXX_VENDOR])
AC_CXX_CPPFLAGS_STD_LANG([CXXFLAGS])
AC_CXX_NAMESPACES
AC_CXX_HAVE_SSTREAM

AC_ARG_ENABLE(warnings,
AS_HELP_STRING([--enable-warnings],[Turn on most compiler warnings]),
[
AX_CXXFLAGS_WARN_ALL_ANSI
AX_CFLAGS_WARN_ALL_ANSI
],
)

case "$ac_cv_cxx_compiler_vendor" in
    compaq)
        AX_CXXFLAGS_OSF_OPTION([-tlocal])
    ;;
    sgi)
        AX_CXXFLAGS_IRIX_OPTION([-ptused])
        AX_CXXFLAGS_IRIX_OPTION([-FE:eliminate_duplicate_inline_copies])
	AX_CXXFLAGS_IRIX_OPTION([-FE:template_in_elf_section])
    ;;
    *)
    ;;
esac

# Talk to Bill about the following block

#case "$build_os" in
#    *cygwin*)
#	CXXFLAGS="${CXXFLAGS} -DCYGWIN"
#	;;
#    *)
#	;;
#esac

# Checks for library functions.

AC_LANG_PUSH([C])
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MEMCMP
AC_CHECK_FUNCS([gettimeofday times strstr floor pow sqrt])
AC_CHECK_LIB([m], [floor])
AC_LANG_POP([C])

OPTPP_OPTIONS

# Output files.

AC_CONFIG_FILES([Makefile			docs/Makefile newmat11/Makefile
		 src/Makefile			src/Base/Makefile
                 src/Constraints/Makefile	src/GSS/Makefile
                 src/Newton/Makefile		src/PDS/Makefile
                 src/UserInterface/Makefile	src/Utils/Makefile
		 tests/Makefile			tests/constraints/Makefile
		 tests/hock/Makefile		tests/npsol/Makefile
		 tests/parallel/Makefile	tests/uncon/Makefile
		 tests/xml/Makefile])
AC_OUTPUT
