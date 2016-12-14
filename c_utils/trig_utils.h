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

/*! \file trig_utils.h
 *
 *  Copyright (C) 2016 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_TRIGHELPER_H
#define PLANCK_TRIGHELPER_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! Computes sine and cosine of \a i*alpha+beta for \a i=[0;n[. Stores the sines
    in \a s[i*stride] and the cosines in c[i*stride]. */
void sincos_multi (size_t n, double alpha, double beta, double *s, double *c,
  int stride);

/*! Computes sine and cosine of \a i*2pi/n for \a i=[0;nang[. Stores the sines
    in \a s[i*stride] and the cosines in c[i*stride]. */
void sincos_2pibyn (size_t n, size_t nang, double *s, double *c, int stride);

#ifdef __cplusplus
}
#endif

#endif
