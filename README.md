# Development has moved

This repository has been archived and is only kept so that packages depending on
it have a cacnonical place to download the last version.

For new projects using spherical harmonic transforms we recommend to use
- `libsharp2` (https://gitlab.mpcdf.mpg.de/mtr/libsharp) if you need MPI
  functionality, or
- the `sht` component of https://gitlab.mpcdf.mpg.de/mtr/ducc if your project
  is written in modern C++ or Python, and you don't need MPI support within
  the transforms.

Both of the above libraries are successors to `libsharp` which are
significantly faster than the last `libsharp` version available here. Since the
switch to the new algorithms required slight changes to the API, we decided to
develop them in separate repositories and keep the original `libsharp`
repository unchanged.


# Libsharp

*IMPORTANT NOTE*: It appears that the default branch upon cloning from
github.com/dagss/libsharp was an outdated 'dagss' branch instead of
the 'master' branch. To get the latest copy,
please do `git checkout master; git pull`. New clones are no longer affected.

## Paper

https://arxiv.org/abs/1303.4945

## Compilation

GNU make is required for compilation.

Libsharp compilation has been successfully tested with GNU and Intel compilers.
When using gcc, version 4.x is required [1].
Since libsharp was written in standard C99, other compilers should work fine,
but SSE2/AVX support will most likely be deactivated.

If you obtained libsharp directly from the git repository, you will also
need a copy of the GNU autotools. In this case, run "autoconf" in libsharp's
main directory before any other steps.
For libsharp releases distributed as a .tar.gz file, this step is not necessary.

Afterwards, simply run "./configure"; if this fails, please refer to the output
of "./configure --help" for additional hints and, if necessary, provide
additional flags to the configure script.
Once the script finishes successfully, run "make"
(or "gmake"). This should install the compilation products in the
subdirectory "auto/".

Documentation can be created by the command "(g)make doc".
However this requires the doxygen application to be installed
on your system.
The documentation will be created in the subdirectory doc/.


[1] Some versions of the gcc 4.4.x release series contain a bug which causes
the compiler to crash during libsharp compilation. This appears to be fixed
in the gcc 4.4.7 release. It is possible to work around this problem by adding
the compiler flag "-fno-tree-fre" after the other optimization flags - the
configure script should do this automatically.
