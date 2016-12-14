/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Functionality for reading wall clock time
 *
 *  Copyright (C) 2010-2016 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#if defined (_OPENMP)
#include <omp.h>
#elif defined (USE_MPI)
#include "mpi.h"
#elif defined (_WIN32)
#include <Windows.h>
#else
#include <sys/time.h>
#include <stdlib.h>
#endif

#include "walltime_c.h"

double wallTime(void)
  {
#if defined (_OPENMP)
  return omp_get_wtime();
#elif defined (USE_MPI)
  return MPI_Wtime();
#elif defined (_WIN32)
  static double inv_freq = -1.;
  if (inv_freq<0)
    {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    inv_freq = 1. / double(freq.QuadPart);
    }
  LARGE_INTEGER count;
  QueryPerformanceCounter(&count);
  return count.QuadPart*inv_freq;
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1e-6*t.tv_usec;
#endif
  }
