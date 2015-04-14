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
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  \author Martin Reinecke \author Dag Sverre Seljebotn
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
  int npairs, nphmax;
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
  /*! Combination of flags from sharp_almflags */
  int flags;
  /*! Array with \a nm entries containing the (hypothetical) indices of
      the coefficients with quantum numbers 0,\a mval[i] */
  ptrdiff_t *mvstart;
  /*! Stride between a_lm and a_(l+1),m */
  ptrdiff_t stride;
  } sharp_alm_info;

/*! alm_info flags */
typedef enum { SHARP_PACKED = 1,
               /*!< m=0-coefficients are packed so that the (zero) imaginary part is
                    not present. mvstart is in units of *real* float/double for all
                    m; stride is in units of reals for m=0 and complex for m!=0 */
               SHARP_REAL_HARMONICS  = 1<<6
               /*!< Use the real spherical harmonic convention. For
                    m==0, the alm are treated exactly the same as in
                    the complex case.  For m!=0, alm[i] represent a
                    pair (+abs(m), -abs(m)) instead of (real, imag),
                    and the coefficients are scaled by a factor of
                    sqrt(2) relative to the complex case.  In other
                    words, (sqrt(.5) * alm[i]) recovers the
                    corresponding complex coefficient (when accessed
                    as complex).
                */
             } sharp_almflags;



/*! Creates an a_lm data structure from the following parameters:
    \param lmax maximum \a l quantum number (>=0)
    \param mmax maximum \a m quantum number (0<= \a mmax <= \a lmax)
    \param stride the stride between entries with identical \a m, and \a l
      differing by 1.
    \param mstart the index of the (hypothetical) coefficient with the
      quantum numbers 0,\a m. Must have \a mmax+1 entries.
    \param alm_info will hold a pointer to the newly created data structure
 */
void sharp_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, sharp_alm_info **alm_info);
/*! Creates an a_lm data structure which from the following parameters:
    \param lmax maximum \a l quantum number (\a >=0)
    \param nm number of different \a m (\a 0<=nm<=lmax+1)
    \param stride the stride between entries with identical \a m, and \a l
      differing by 1.
    \param mval array with \a nm entries containing the individual m values
    \param mvstart array with \a nm entries containing the (hypothetical)
      indices of the coefficients with the quantum numbers 0,\a mval[i]
    \param flags a combination of sharp_almflags (pass 0 unless you know you need this)
    \param alm_info will hold a pointer to the newly created data structure
 */
void sharp_make_general_alm_info (int lmax, int nm, int stride, const int *mval,
  const ptrdiff_t *mvstart, int flags, sharp_alm_info **alm_info);
/*! Returns the index of the coefficient with quantum numbers \a l,
    \a mval[mi].
    \note for a \a sharp_alm_info generated by sharp_make_alm_info() this is
    the index for the coefficient with the quantum numbers \a l, \a mi. */
ptrdiff_t sharp_alm_index (const sharp_alm_info *self, int l, int mi);
/*! Returns the number of alm coefficients described by \a self. If the SHARP_PACKED
    flag is set, this is number of "real" coeffecients (for m < 0 and m >= 0),
    otherwise it is the number of complex coefficients (with m>=0). */
ptrdiff_t sharp_alm_count(const sharp_alm_info *self);
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
    \param wgt the pixel weight to be used for the ring in map2alm
      and adjoint map2alm transforms.
      Pass NULL to use 1.0 as weight for all rings.
    \param geom_info will hold a pointer to the newly created data structure
 */
void sharp_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *wgt, sharp_geom_info **geom_info);

/*! Counts the number of grid points needed for (the local part of) a map described
    by \a info.
 */
ptrdiff_t sharp_map_size(const sharp_geom_info *info);

/*! Deallocates the geometry information in \a info. */
void sharp_destroy_geom_info (sharp_geom_info *info);

/*! \} */

/*! \defgroup lowlevelgroup Low-level libsharp SHT interface */
/*! \{ */

/*! Enumeration of SHARP job types. */
typedef enum { SHARP_YtW=0,               /*!< analysis */
               SHARP_MAP2ALM=SHARP_YtW,   /*!< analysis */
               SHARP_Y=1,                 /*!< synthesis */
               SHARP_ALM2MAP=SHARP_Y,     /*!< synthesis */
               SHARP_Yt=2,                /*!< adjoint synthesis */
               SHARP_WY=3,                /*!< adjoint analysis */
               SHARP_ALM2MAP_DERIV1=4     /*!< synthesis of first derivatives */
             } sharp_jobtype;

/*! Job flags */
typedef enum { SHARP_DP              = 1<<4,
               /*!< map and a_lm are in double precision */
               SHARP_ADD             = 1<<5,
               /*!< results are added to the output arrays, instead of
                    overwriting them */

               /* NOTE: SHARP_REAL_HARMONICS, 1<<6, is also available in sharp_jobflags,
                  but its use here is deprecated in favor of having it in the sharp_alm_info */

               SHARP_NO_FFT          = 1<<7,

               SHARP_USE_WEIGHTS     = 1<<20,    /* internal use only */
               SHARP_NO_OPENMP       = 1<<21,    /* internal use only */
               SHARP_NVMAX           = (1<<4)-1 /* internal use only */
             } sharp_jobflags;

/*! Performs a libsharp SHT job. The interface deliberately does not use
  the C99 "complex" data type, in order to be callable from C89 and C++.
  \param type the type of SHT
  \param spin the spin of the quantities to be transformed
  \param alm contains pointers to the a_lm coefficients. If \a spin==0,
    alm[0] points to the a_lm of the first SHT, alm[1] to those of the second
    etc. If \a spin>0, alm[0] and alm[1] point to the a_lm of the first SHT,
    alm[2] and alm[3] to those of the second, etc. The exact data type of \a alm
    depends on whether the SHARP_DP flag is set.
  \param map contains pointers to the maps. If \a spin==0,
    map[0] points to the map of the first SHT, map[1] to that of the second
    etc. If \a spin>0, or \a type is SHARP_ALM2MAP_DERIV1, map[0] and map[1]
    point to the maps of the first SHT, map[2] and map[3] to those of the
    second, etc. The exact data type of \a map depends on whether the SHARP_DP
    flag is set.
  \param geom_info A \c sharp_geom_info object compatible with the provided
    \a map arrays.
  \param alm_info A \c sharp_alm_info object compatible with the provided
    \a alm arrays. All \c m values from 0 to some \c mmax<=lmax must be present
    exactly once.
  \param ntrans the number of simultaneous SHTs
  \param flags See sharp_jobflags. In particular, if SHARP_DP is set, then
    \a alm is expected to have the type "complex double **" and \a map is
    expected to have the type "double **"; otherwise, the expected
    types are "complex float **" and "float **", respectively.
  \param time If not NULL, the wall clock time required for this SHT
    (in seconds) will be written here.
  \param opcnt If not NULL, a conservative estimate of the total floating point
    operation count for this SHT will be written here. */
void sharp_execute (sharp_jobtype type, int spin, void *alm, void *map,
  const sharp_geom_info *geom_info, const sharp_alm_info *alm_info, int ntrans,
  int flags, double *time, unsigned long long *opcnt);

void sharp_set_chunksize_min(int new_chunksize_min);
void sharp_set_nchunks_max(int new_nchunks_max);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
