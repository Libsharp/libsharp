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

ENABLE_SSE2=yes
AC_ARG_ENABLE(sse2,
  [  --disable-sse2          disable generation of SSE2 instructions],
  [if test "$enableval" = no; then
     ENABLE_SSE2=no
   fi])

ENABLE_AVX=yes
AC_ARG_ENABLE(avx,
  [  --disable-avx           disable generation of AVX instructions],
  [if test "$enableval" = no; then
     ENABLE_AVX=no
   fi])

case $CCTYPE in
  gcc)
    CCFLAGS="-O3 -fno-tree-vectorize -ffast-math -fomit-frame-pointer -std=c99 -pedantic -Wextra -Wall -Wno-unknown-pragmas -Wshadow -Wmissing-prototypes -Wfatal-errors"
    GCCVERSION="`$CC -dumpversion 2>&1`"
    echo "Using gcc version $GCCVERSION"
    AC_SUBST(GCCVERSION)
    case $system in
      Darwin-*)
        ;;
      *)
        CCFLAGS="$CCFLAGS -ffunction-sections -fdata-sections"
        ;;
    esac
    changequote(,)
    gcc43=`echo $GCCVERSION | grep -c '4\.[3456789]'`
    changequote([,])
    if test $gcc43 -gt 0; then
      CCFLAGS="$CCFLAGS -march=native"
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

CCFLAGS="$CCFLAGS $OPENMP_CFLAGS"

if test $ENABLE_DEBUG = yes; then
  CCFLAGS="$CCFLAGS -g"
fi

if test $ENABLE_MPI = yes; then
  CCFLAGS="$CCFLAGS -DUSE_MPI"
fi

if test $ENABLE_SSE2 = no; then
  CCFLAGS="$CCFLAGS -DDISABLE_SSE2"
fi

if test $ENABLE_AVX = no; then
  CCFLAGS="$CCFLAGS -DDISABLE_AVX"
fi

CCFLAGS_NO_C="$CCFLAGS $CPPFLAGS"

LDCCFLAGS="$LDFLAGS $CCFLAGS"

AC_SUBST(SILENT_RULE)
AC_SUBST(CC)
AC_SUBST(CCFLAGS_NO_C)
AC_SUBST(LDCCFLAGS)
AC_SUBST(ARCREATE)

AC_OUTPUT(config/config.auto)
