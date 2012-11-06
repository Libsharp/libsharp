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
    \param weight the pixel weight to be used for the ring. Pass NULL to use
      1.0 as weight for all rings. By default weights are used for alm2map
      but not map2alm, but the execution flags can override this.
    \param geom_info will hold a pointer to the newly created data structure
 */
void sharp_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, sharp_geom_info **geom_info);

/*! Deallocates the geometry information in \a info. */
void sharp_destroy_geom_info (sharp_geom_info *info);

/*! \} */

/*! \defgroup lowlevelgroup Low-level libsharp SHT interface */
/*! \{ */

/*! Enumeration of SHARP job types. */
typedef enum { SHARP_MAP2ALM,       /*!< analysis */
               SHARP_ALM2MAP,       /*!< synthesis */
               SHARP_ALM2MAP_DERIV1 /*!< synthesis of first derivatives */
             } sharp_jobtype;

/*! Job flags */
typedef enum { SHARP_SP = 0, /*!< map and alm is in single precision */
               SHARP_DP = 1 << 1, /*!< map and alm is in double precision */

               SHARP_ALM2MAP_USE_WEIGHTS = 1 << 2, /*!< apply ring weights for alm2map */
               SHARP_MAP2ALM_IGNORE_WEIGHTS = 1 << 3, /*!< do not use ring weights for map2alm */

               /* convenience flag combinations (stable API even if the default changes) */
               SHARP_USE_WEIGHTS = SHARP_ALM2MAP_USE_WEIGHTS, /*!< use ring weights for both map2alm and alm2map */
               SHARP_IGNORE_WEIGHTS = SHARP_MAP2ALM_IGNORE_WEIGHTS /*!< do not use ring weights for either map2alm or map2alm */
             } sharp_jobflags;

/*! Performs a libsharp SHT job. The interface deliberately does not use
  the C99 "complex" data type, in order to be callable from C. 
  \param type the type of SHT
  \param spin the spin of the quantities to be transformed
  \param add_output if 0, the output arrays will be overwritten,
    else the result will be added to the output arrays.
  \param alm contains pointers to the a_lm coefficients. If \a spin==0,
    alm[0] points to the a_lm of the first SHT, alm[1] to those of the second
    etc. If \a spin>0, alm[0] and alm[1] point to the a_lm of the first SHT,
    alm[2] and alm[3] to those of the second, etc. The exact data type of \a alm
    depends on the \a dp parameter.
  \param map contains pointers to the maps. If \a spin==0,
    map[0] points to the map of the first SHT, map[1] to that of the second
    etc. If \a spin>0, map[0] and map[1] point to the maps of the first SHT,
    map[2] and map[3] to those of the second, etc. The exact data type of \a map
    depends on the \a dp parameter.
  \param geom_info A \c sharp_geom_info object compatible with the provided
    \a map arrays.
  \param alm_info A \c sharp_alm_info object compatible with the provided
    \a alm arrays. All \c m values from 0 to some \c mmax<=lmax must be present
    exactly once.
  \param ntrans the number of simultaneous SHTs
  \param flags See sharp_jobflags. In particular, if SHARP_SP is set, the \a alm is
    expected to have the type "complex float **"
    and \a map is expected to have the type "float **"; if SHARP_DP is set, the expected
    types are "complex double **" and "double **", respectively.
  \param nv Internally used SHT parameter. Set to 0 unless you know what you are
    doing.
  \param time If not NULL, the wall clock time required for this SHT
    (in seconds)will be written here.
  \param opcnt If not NULL, a conservative estimate of the total floating point
    operation count for this SHT will be written here. */
void sharp_execute (sharp_jobtype type, int spin, int add_output, void *alm,
  void *map, const sharp_geom_info *geom_info, const sharp_alm_info *alm_info,
  int ntrans, int flags, int nv, double *time, unsigned long long *opcnt);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
