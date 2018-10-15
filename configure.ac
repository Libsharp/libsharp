AC_INIT(config/config.auto.in)

AC_CHECK_PROG([uname_found],[uname],[1],[0])
if test $uname_found -eq 0 ; then
    echo "No uname found; setting system type to unknown."
    system="unknown"
else
    system=`uname -s`-`uname -r`
fi
AC_LANG([C])

AC_TRY_COMPILE([], [@%:@ifndef __INTEL_COMPILER
choke me
@%:@endif], [ICC=[yes]], [ICC=[no]])

if test $ICC = yes; then GCC=no; fi
CCTYPE=unknown
if test $GCC = yes; then CCTYPE=gcc; fi
if test $ICC = yes; then CCTYPE=icc; fi
AC_OPENMP

SILENT_RULE=".SILENT:"
AC_ARG_ENABLE(noisy-make,
  [  --enable-noisy-make     enable detailed make output],
  [if test "$enableval" = yes; then
     SILENT_RULE=""
   fi])

ENABLE_MPI=no
AC_ARG_ENABLE(mpi,
  [  --enable-mpi            enable generation of MPI-parallel code],
  [if test "$enableval" = yes; then
     ENABLE_MPI=yes
   fi])

ENABLE_DEBUG=no
AC_ARG_ENABLE(debug,
  [  --enable-debug          enable generation of debugging symbols],
  [if test "$enableval" = yes; then
     ENABLE_DEBUG=yes
   fi])

ENABLE_PIC=no
AC_ARG_ENABLE(pic,
  [  --enable-pic            enable generation of position independent code],
  [if test "$enableval" = yes; then
     ENABLE_PIC=yes
   fi])

case $CCTYPE in
  gcc)
    CCFLAGS="-O3 -fno-tree-vectorize -ffast-math -fomit-frame-pointer -std=c99 -pedantic -Wextra -Wall -Wno-unknown-pragmas -Wshadow -Wmissing-prototypes -Wfatal-errors -march=native"
    GCCVERSION="`$CC -dumpversion 2>&1`"
    echo "Using gcc version $GCCVERSION"
    AC_SUBST(GCCVERSION)
    changequote(,)
    gcc43=`echo $GCCVERSION | grep -c '^4\.[3456789]'`
    gcc44=`echo $GCCVERSION | grep -c '^4\.4'`
    changequote([,])
    if test $gcc44 -gt 0; then
      CCFLAGS="$CCFLAGS -fno-tree-fre"
    fi
    ;;
  icc)
    CCFLAGS="-O3 -xHOST -std=c99 -ip -Wbrief -Wall -vec-report0 -openmp-report0 -wd383,981,1419,1572"
    ;;
  *)
    CCFLAGS="-O2"
    # Don't do anything now
    ;;
esac

case $system in
  Darwin-*)
    ARCREATE="libtool -static -o"
    ;;
  *)
    ARCREATE="ar cr"
    ;;
esac

if test $ENABLE_DEBUG = yes; then
  DEBUG_CFLAGS="-g"
fi

if test $ENABLE_PIC = yes; then
  PIC_CFLAGS="-fPIC"
fi

if test $ENABLE_MPI = yes; then
  MPI_CFLAGS="-DUSE_MPI"
fi

CCFLAGS="$CCFLAGS $DEBUG_CFLAGS $OPENMP_CFLAGS $PIC_CFLAGS $MPI_CFLAGS"

CCFLAGS_NO_C="$CCFLAGS $CPPFLAGS"

LDCCFLAGS="$LDFLAGS $CCFLAGS"

AC_SUBST(SILENT_RULE)
AC_SUBST(CC)
AC_SUBST(CCFLAGS_NO_C)
AC_SUBST(LDCCFLAGS)
AC_SUBST(DEBUG_CFLAGS)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(PIC_CFLAGS)
AC_SUBST(ARCREATE)

AC_OUTPUT(config/config.auto)
