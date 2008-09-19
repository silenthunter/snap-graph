# run autoconf to produce the configure script
AC_PREREQ(2.59)

AC_INIT(snap, 0.3, kmadduri@lbl.gov)
AC_CONFIG_SRCDIR(src/graph_kernels/breadth_first_search.c)
AM_CONFIG_HEADER(include/graph_defs.h)

AM_INIT_AUTOMAKE

# get platform information
AC_CANONICAL_HOST

AC_ARG_ENABLE(debug,
        [AC_HELP_STRING([--enable-debug],
            [build library and binaries with extra runtime checks for debugging])], 
        enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
    AC_DEFINE(SNAP_DEBUG,1,[Define to build snap with debug flags])
fi

# check for C compiler
AC_PROG_CC

# macro for determining compiler vendor
AX_COMPILER_VENDOR

# check for libtool
AC_PROG_LIBTOOL

# macro for determining optimization flags
AX_CC_MAXOPT

# override CFLAGS selection when debug option is selected
if test "${enable_debug}" = "yes"; then
	CFLAGS="-g"
fi

# add more warning flags for gcc, in debug/maintainer mode only
if test "$enable_debug" = "yes" || test "$USE_MAINTAINER_MODE" = "yes"; then
if test "${ac_test_CFLAGS}" != "set"; then
	if test "${ac_cv_prog_gcc}" = "yes"; then
		CFLAGS="$CFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wno-long-long -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs" # -Wundef -Wconversion -Wmissing-prototypes -Wmissing-declarations 
	fi
fi
fi

# check for the appropriate OpenMP flags and enable it.
AC_ARG_ENABLE(openmp,
        [AC_HELP_STRING([--enable-openmp],
            [use OpenMP directives for parallelism])],
        enable_openmp=$enableval, enable_openmp=no)
if test "$enable_openmp" = "yes"; then
    AC_DEFINE(HAVE_OPENMP,1,[Define to enable OpenMP])
fi

if test "$enable_openmp" = "yes"; then
   AX_OPENMP([THREADLIBS=" "
     AC_DEFINE(USING_OPENMP_THREADS, 1, [Define if we have and are using OpenMP multithreading directives])
     CFLAGS="$CFLAGS $OPENMP_CFLAGS"],
     [AC_MSG_ERROR([Unable to determine the appropriate OpenMP flag for your compiler.
         Please set the CFLAGS variable manually for OpenMP support.])])
fi

# checks for libraries
AC_CHECK_LIB(m, log2, [], [])

# header files
AC_CHECK_HEADERS([stdarg.h limits.h])

# typedefs, structures, and compiler characteristics
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# check functions
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor gettimeofday gethrtime memset memalign sqrt drand48 pow])

# list of config files
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/graph_generation/Makefile
                 src/graph_kernels/Makefile
                 src/graph_partitioning/Makefile
                 src/graph_metrics/Makefile
                 src/sprng2.0/Makefile
                 src/utils/Makefile
                 lib/Makefile
                 test/Makefile
                ])
AC_OUTPUT

# Print configuration output.
AC_MSG_RESULT([ ])
AC_MSG_RESULT([snap successfully configured! Please verify that])
AC_MSG_RESULT([this configuration matches with your expectations.])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([   User Option                         Value          ])
AC_MSG_RESULT([------------------------------------------------------])
AC_MSG_RESULT([   OpenMP support                      $enable_openmp])
AC_MSG_RESULT([   Enable debug                        $enable_debug])
AC_MSG_RESULT([])
AC_MSG_RESULT([   Compiler Configuration:])
AC_MSG_RESULT([------------------------------------------------------])
AC_MSG_RESULT([   C compiler                          $CC      ])
AC_MSG_RESULT([   Flags                               $CFLAGS  ])
AC_MSG_RESULT([])
