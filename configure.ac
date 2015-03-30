### configure.ac                                        -*- Autoconf -*-
###
### Process this file with autoconf to produce a configure script.
###
### ReliefSeq
### Bill White (c) 2005-2012
### http://insilico.utulsa.edu/reliefseq

AC_PREREQ(2.61)

AC_INIT([reliefseq], [1.0], [bill-white@utulsa.edu], [http://insilico.utulsa.edu/reliefseq])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([ReliefSeqCLI.cpp])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([reliefseq], [1.0])

AC_SUBST(SRCDIR, src)
AC_SUBST(LIBDIR, library)

AC_CANONICAL_BUILD
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
#AC_F77_LIBRARY_LDFLAGS
AC_PROG_LN_S
AC_PROG_LIBTOOL

# Use libtool for all builds
LT_INIT

# Set variables for basic compile debug/release
AC_ARG_ENABLE(debug, 
AC_HELP_STRING([--enable-debug],[Compile code for debugging.]), 
[debug_reliefseq=yes]
)

#CXXFLAGS="-O3 -Wall -I/usr/local/include -I/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin13.4.0/4.9.2/include $CXXFLAGS"
CXXFLAGS="-O3 -Wall -I/usr/local/include $CXXFLAGS"
if test "$debug_reliefseq" == 'yes' ; then
    CXXFLAGS="-O0 -Wall -DDEBUG"
fi

AC_ARG_ENABLE(static,
[  --enable-static    Turn on static builds],
[case "${enableval}" in
  yes) static=true ;;
  no)  static=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-static) ;;
esac],
[static=false]
)
AM_CONDITIONAL(STATIC, test x$static = xtrue)

# Checks for header files.
AC_CHECK_HEADERS([math.h stdlib.h string.h sys/time.h omp.h])
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_SIZEOF(void *)# AC_CHECK_FUNCS([round pow sqrt floor memset], [], AC_MSG_ERROR([required function missing]))

# makeinfo
AC_CHECK_PROG([MAKEINFO], [makeinfo], [makeinfo], AC_MSG_ERROR([could not find required program: makeinfo]))

AC_MSG_NOTICE([--------------------------------------------------])

LDFLAGS="-L/usr/local/lib -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9 $LDFLAGS"
AC_CHECK_LIB([gomp], [omp_get_num_threads],
    [], AC_MSG_ERROR(test for OpenMP failed))
LDFLAGS="$LDFLAGS -lgomp"

#AC_CHECK_LIB(pthread, pthread_create)
#LDFLAGS="$LDFLAGS -lpthread"

# boost.m4
BOOST_REQUIRE
BOOST_PROGRAM_OPTIONS
BOOST_FIND_HEADER(boost/program_options.hpp, AC_MSG_ERROR([could not find required headers of Boost c++. Install boost c++ libraries.]),)

# gsl.m4
AX_PATH_GSL(1.8,, AC_MSG_ERROR([could not find required version of GSL]),)

######################################################################

AC_MSG_NOTICE([--------------------------------------------------])

# Process makefiles
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
