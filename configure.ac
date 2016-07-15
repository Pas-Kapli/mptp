#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([mptp], [0.1.1], [Tomas.Flouri@h-its.org])
AM_INIT_AUTOMAKE([subdir-objects])
AC_LANG([C])
AC_CONFIG_SRCDIR([src/mptp.c])
AC_CONFIG_HEADERS([config.h])
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_SED
AC_PROG_LEX
if test "x$LEX" != xflex; then
  AC_MSG_ERROR(could not find required installation of FLEX)
fi

AC_PROG_YACC
if test "x$YACC" != x"bison -y"; then
  AC_MSG_ERROR(could not find required installation of BISON)
fi

AC_PROG_INSTALL

# Checks for header files.
AC_CHECK_HEADERS([assert.h stdio.h stdarg.h string.h getopt.h stdlib.h regex.h ctype.h locale.h limits.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_FUNC_ALLOCA
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove memcpy gettimeofday memchr memset pow regcomp strcasecmp strchr strcspn sysinfo])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_cdf_chisq_P])

# Bash completions
AC_ARG_WITH([bash-completions],
  AC_HELP_STRING([--with-bash-completions=[DIR]], [Bash completions directory [default=no]]),
  [with_bash_completions="$withval"],
  [with_bash_completions="no"]
)
AS_CASE([$with_bash_completions],
#  [yes], [PKG_CHECK_VAR([bash_completions_dir], [bash-completion], [completionsdir], [], [AC_MSG_ERROR([bash completions not found])])],
  [no], [bash_completions_dir=],
  [bash_completions_dir="$with_bash_completions"]
)
AC_SUBST([bash_completions_dir])
AM_CONDITIONAL(HAVE_BASH_COMPLETIONS, test -n "$bash_completions_dir")
AS_IF([test -n "$bash_completions_dir"],
  [bash_completions_output="${bash_completions_dir}"],
  [bash_completions_output=no]
)

have_ps2pdf=no
AC_ARG_ENABLE(pdfman, AS_HELP_STRING([--disable-pdfman], [Disable PDF manual creation]))
AS_IF([test "x$enable_pdfman" != "xno"], [
  have_ps2pdf=yes
  AC_CHECK_PROG(HAVE_PS2PDF, ps2pdf, yes, no)
  if test "x$HAVE_PS2PDF" = "xno"; then
    AC_MSG_WARN([*** ps2pdf is required to build a PDF version of the manual])
    have_ps2pdf=no
  fi
])

AM_CONDITIONAL(HAVE_PS2PDF, test "x${have_ps2pdf}" = "xyes")
AM_PROG_CC_C_O

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile
                 completion/Makefile])

AC_OUTPUT

AC_MSG_RESULT([
  $PACKAGE $VERSION

    Target:    $host_os $host_cpu
    Compiler:  ${CC}
    CFLAGS:    ${CFLAGS} ${CPPFLAGS}
    LIBS:      ${LIBS} ${LDFLAGS}
  
  Continue  with 'make' command
])
