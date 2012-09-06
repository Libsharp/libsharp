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

/*! \file sharp_lowlevel.h
 *  Low-level, portable interface for the spherical transform library.
 *
 *  Copyright (C) 2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SHARP_LOWLEVEL_H
#define PLANCK_SHARP_LOWLEVEL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \internal
    Helper type containing information about a single ring. */
typedef struct
  {
  double theta, phi0, weight, cth, sth;
  ptrdiff_t ofs;
  int nph, stride;
  } sharp_ringinfo;

/*! \internal
    Helper type containing information about a pair of rings with colatitudes
    symmetric around the equator. */
typedef struct
  {
  sharp_ringinfo r1,r2;
  } sharp_ringpair;

/*! \internal
    Type holding all required information about a map geometry. */
typedef struct
  {
  sharp_ringpair *pair;
  int npairs;
  } sharp_geom_info;

/*! \defgroup almgroup Helpers for dealing with a_lm */
/*! \{ */

/*! \internal
    Helper type for index calculation in a_lm arrays. */
typedef struct
  {
  /*! Maximum \a l index of the array */
  int lmax;
  /*! Number of different \a m values in this object */
  int nm;
  /*! Array with \a nm entries containing the individual m values */
  int *mval;
  /*! Array with \a nm entries containing the (hypothetical) indices of
      the coefficients with quantum numbers 0,\a mval[i] */
  ptrdiff_t *mvstart;
  /*! Stride between a_lm and a_(l+1),m */
  ptrdiff_t stride;
  } sharp_alm_info;

/*! Creates an Alm data structure information from the following parameters:
    \param lmax maximum \a l quantum number (>=0)
    \param mmax maximum \a m quantum number (0<= \a mmax <= \a lmax)
    \param stride the stride between consecutive a_lm entries
    \param mstart the index of the (hypothetical) coefficient with the
      quantum numbers 0,\a m. Must have \a mmax+1 entries.
    \param alm_info will hold a pointer to the newly created data structure
 */
void sharp_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, sharp_alm_info **alm_info);
/*! Creates an Alm data structure information from the following parameters:
    \param lmax maximum \a l quantum number (>=0)
    \param nm number of different \a m (<=\a lmax+1)
    \param stride the stride between consecutive a_lm entries
    \param mval array with \a nm entries containing the individual m values
    \param mvstart array with \a nm entries containing the (hypothetical)
      indices of the coefficients with the quantum numbers 0,\a mval[i]
    \param alm_info will hold a pointer to the newly created data structure
 */
void sharp_make_general_alm_info (int lmax, int nm, int stride, const int *mval,
  const ptrdiff_t *mvstart, sharp_alm_info **alm_info);
/*! Returns the index of the coefficient with quantum numbers \a l,
    \a mval[mi]. */
ptrdiff_t sharp_alm_index (const sharp_alm_info *self, int l, int mi);
/*! Deallocates the a_lm info object. */
void sharp_destroy_alm_info (sharp_alm_info *info);

/*! \} */

/*! \defgroup geominfogroup Functions for dealing with geometry information */
/*! \{ */

/*! Creates a geometry information from a set of ring descriptions.
    All arrays passed to this function must have \a nrings elements.
    \param nrings the number of rings in the map
    \param nph the number of pixels in each ring
    \param ofs the index of the first pixel in each ring in the map array
    \param stride the stride between consecutive pixels
    \param phi0 the azimuth (in radians) of the first pixel in each ring
    \param theta the colatitude (in radians) of each ring
    \param weight the pixel weight to be used for the ring
    \param geom_info will hold a pointer to the newly created data structure
 */
void sharp_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, sharp_geom_info **geom_info);

/*! Deallocates the geometry information in \a info. */
void sharp_destroy_geom_info (sharp_geom_info *info);

/*! \} */

/*! \defgroup jobgroup Functionality for defining and executing SHTs */
/*! \{ */

/*! Enumeration of SHARP job types. */
typedef enum { SHARP_MAP2ALM,       /*!< analysis */
               SHARP_ALM2MAP,       /*!< synthesis */
               SHARP_ALM2MAP_DERIV1 /*!< synthesis of first derivatives */
             } sharp_jobtype;

void sharp_execute_ll (sharp_jobtype type, int spin, int add_output, void **alm,
  void **map, const sharp_geom_info *geom_info, const sharp_alm_info *alm_info,
  int ntrans, int dp);

#ifdef __cplusplus
}
#endif

#endif
