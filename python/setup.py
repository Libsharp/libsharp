#! /usr/bin/env python

descr   = """Spherical Harmionic transforms package

Python API for the libsharp spherical harmonic transforms library
"""

import os
import sys

DISTNAME            = 'libsharp'
DESCRIPTION         = 'libsharp library for fast Spherical Harmonic Transforms'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'Dag Sverre Seljebotn',
MAINTAINER_EMAIL    = 'd.s.seljebotn@astro.uio.no',
URL                 = 'http://sourceforge.net/projects/libsharp/'
LICENSE             = 'GPL'
DOWNLOAD_URL        = "http://sourceforge.net/projects/libsharp/"
VERSION             = '0.1'

# Add our fake Pyrex at the end of the Python search path
# in order to fool setuptools into allowing compilation of
# pyx files to C files. Importing Cython.Distutils then
# makes Cython the tool of choice for this rather than
# (the possibly nonexisting) Pyrex.
project_path = os.path.split(__file__)[0]
sys.path.append(os.path.join(project_path, 'fake_pyrex'))

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import numpy as np

libsharp = os.environ.get('LIBSHARP', None)
libsharp_include = os.environ.get('LIBSHARP_INCLUDE', libsharp and os.path.join(libsharp, 'include'))
libsharp_lib = os.environ.get('LIBSHARP_LIB', libsharp and os.path.join(libsharp, 'lib'))

if libsharp_include is None or libsharp_lib is None:
    sys.stderr.write('Please set LIBSHARP environment variable to the install directly of libsharp, '
                     'this script will refer to the lib and include sub-directories. Alternatively '
                     'set LIBSHARP_INCLUDE and LIBSHARP_LIB\n')
    sys.exit(1)

if __name__ == "__main__":
    setup(install_requires = ['numpy'],
          packages = find_packages(),
          test_suite="nose.collector",
          # Well, technically zipping the package will work, but since it's
          # all compiled code it'll just get unzipped again at runtime, which
          # is pointless:
          zip_safe = False,
          name = DISTNAME,
          version = VERSION,
          maintainer = MAINTAINER,
          maintainer_email = MAINTAINER_EMAIL,
          description = DESCRIPTION,
          license = LICENSE,
          url = URL,
          download_url = DOWNLOAD_URL,
          long_description = LONG_DESCRIPTION,
          classifiers =
            [ 'Development Status :: 3 - Alpha',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Topic :: Scientific/Engineering'],
          cmdclass = {"build_ext": build_ext},
          ext_modules = [
              Extension("libsharp.libsharp",
                        ["libsharp/libsharp.pyx"],
                        libraries=["sharp", "fftpack", "c_utils"],
                        include_dirs=[libsharp_include],
                        library_dirs=[libsharp_lib],
                        extra_link_args=["-fopenmp"],
              ),
              Extension("libsharp.libsharp_mpi",
                        ["libsharp/libsharp_mpi.pyx"],
                        libraries=["sharp", "fftpack", "c_utils"],
                        include_dirs=[libsharp_include],
                        library_dirs=[libsharp_lib],
                        extra_link_args=["-fopenmp"],
              ),
              ],
          )
