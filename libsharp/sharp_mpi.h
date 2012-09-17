/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file sharp_mpi.h
 *  Interface for the spherical transform library with MPI support.
 *
 *  Copyright (C) 2011,2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SHARP_MPI_H
#define PLANCK_SHARP_MPI_H

#include <mpi.h>
#include "sharp.h"

#ifdef __cplusplus
extern "C" {
#endif

void sharp_execute_mpi (MPI_Comm comm, sharp_jobtype type, int spin,
  int add_output, void **alm, void **map, const sharp_geom_info *geom_info,
  const sharp_alm_info *alm_info, int ntrans, int dp, int nv, double *time,
  unsigned long long *opcnt);

#ifdef __cplusplus
}
#endif

#endif
