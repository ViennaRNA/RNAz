AC_INIT(RNAz, 0.1.1, wash@tbi.univie.ac.at, RNAz)
AM_INIT_AUTOMAKE

AC_PROG_CC

if test "$GCC" = yes; then
  AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
  AC_DEFINE(UNUSED,, only for gcc)
fi
      
AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CHECK_FUNCS(strdup strstr strchr erand48)

AC_PROG_CXX

AC_PROG_RANLIB

AC_PROG_INSTALL

AC_CONFIG_HEADERS(config.h) 

AC_OUTPUT(Makefile librna/Makefile libsvm/Makefile rnaz/Makefile models/Makefile examples/Makefile man/Makefile)


