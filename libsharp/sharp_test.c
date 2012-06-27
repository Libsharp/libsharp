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

/*! \file sharp_test.c
    Accuracy test for libsharp's map analysis.

    This program first generates a_lm coefficients up to
    a user-specified lmax (with mmax=lmax); where applicable, the
    real and imaginary parts of the coefficients are uniform
    random numbers of the interval [-1;1[.
    Afterwards, the random a_lm are converted to a map.
    This map is analyzed (optionally using an iterative scheme
    with a user-supplied number of steps).
    After every iteration, the code then outputs the RMS of the residual a_lm
    (i.e. the difference between the current and original a_lm), divided by
    the RMS of the original a_lm, as well as the maximum absolute change of any
    real or imaginary part between the current and original a_lm.

    This operation can be performed for several different pixelisations:
      - a Gaussian with the minimal number of rings for exact analysis
        and a user-defined ring resolution
      - an ECP grid with the minimal number of rings for exact analysis
        and a user-defined ring resolution
      - a Healpix grid with a user-defined Nside parameter.

    The user can specify the spin of the desired transform.

    Copyright (C) 2006-2012 Max-Planck-Society
    \author Martin Reinecke
*/

#include <stdio.h>
#include <string.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "sharp.h"
#include "sharp_geomhelpers.h"
#include "sharp_almhelpers.h"
#include "c_utils.h"
#include "sharp_core.h"

typedef complex double dcmplx;

static double drand (double min, double max)
  { return min + (max-min)*rand()/(RAND_MAX+1.0); }

static void random_alm (dcmplx *alm, sharp_alm_info *helper, int spin)
  {
  for (int mi=0;mi<helper->nm; ++mi)
    {
    int m=helper->mval[mi];
    for (int l=m;l<=helper->lmax; ++l)
      {
      if ((l<spin)&&(m<spin))
        alm[sharp_alm_index(helper,l,mi)] = 0.;
      else
        {
        double rv = drand(-1,1);
        double iv = (m==0) ? 0 : drand(-1,1);
        alm[sharp_alm_index(helper,l,mi)] = rv+_Complex_I*iv;
        }
      }
    }
  }

static void measure_errors (dcmplx **alm, dcmplx **alm2,
  ptrdiff_t nalms, int ncomp)
  {
  for (int i=0; i<ncomp; ++i)
    {
    double sum=0, sum2=0, maxdiff=0;
    for (ptrdiff_t m=0; m<nalms; ++m)
      {
      double x=creal(alm[i][m])-creal(alm2[i][m]),
             y=cimag(alm[i][m])-cimag(alm2[i][m]);
      sum+=x*x+y*y;
      sum2+=creal(alm[i][m])*creal(alm[i][m])+cimag(alm[i][m])*cimag(alm[i][m]);
      if (fabs(x)>maxdiff) maxdiff=fabs(x);
      if (fabs(y)>maxdiff) maxdiff=fabs(y);
      }
    sum=sqrt(sum/nalms);
    sum2=sqrt(sum2/nalms);
    printf("component %i: rms %e, maxerr %e\n",i, sum/sum2, maxdiff);
    }
  }

static void map2alm_iter (sharp_geom_info *tinfo, double **map,
  dcmplx **alm_orig, dcmplx **alm, int lmax, int mmax,
  ptrdiff_t npix, ptrdiff_t nalms, int spin, int ntrans, int niter)
  {
  int ncomp = ntrans*((spin==0) ? 1 : 2);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);

  sharp_job job;
  sharpd_build_job(&job,MAP2ALM,spin,0,&alm[0],&map[0],tinfo,alms,ntrans);
  sharp_execute_job(&job);
  printf("wall time for map2alm: %fs\n",job.time);
  printf("Performance: %fGFLOPs/s\n",1e-9*job.opcnt/job.time);
  measure_errors(alm_orig,alm,nalms,ncomp);

  for (int iter=0; iter<niter; ++iter)
    {
    double **map2;
    ALLOC2D(map2,double,ncomp,npix);
    printf ("\niteration %i:\n", iter+1);
    sharpd_build_job(&job,ALM2MAP,spin,0,&alm[0],&map2[0],tinfo,alms,ntrans);
    sharp_execute_job(&job);
    printf("wall time for alm2map: %fs\n",job.time);
    printf("Performance: %fGFLOPs/s\n",1e-9*job.opcnt/job.time);
    for (int i=0; i<ncomp; ++i)
      for (ptrdiff_t m=0; m<npix; ++m)
        map2[i][m] = map[i][m]-map2[i][m];

    sharpd_build_job(&job,MAP2ALM,spin,1,&alm[0],&map2[0],tinfo,alms,ntrans);
    sharp_execute_job(&job);
    printf("wall time for map2alm: %fs\n",job.time);
    printf("Performance: %fGFLOPs/s\n",1e-9*job.opcnt/job.time);
    DEALLOC2D(map2);
    measure_errors(alm_orig,alm,nalms,ncomp);
    }

  sharp_destroy_alm_info(alms);
  }

static void check_accuracy (sharp_geom_info *tinfo, ptrdiff_t lmax,
  ptrdiff_t mmax, ptrdiff_t npix, int spin, int ntrans, int niter)
  {
  ptrdiff_t nalms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);
  int ncomp = ntrans*((spin==0) ? 1 : 2);

  double **map;
  ALLOC2D(map,double,ncomp,npix);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);

  srand(4);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,ncomp,nalms);
  for (int i=0; i<ncomp; ++i)
    random_alm(alm[i],alms,spin);

  dcmplx **alm2;
  ALLOC2D(alm2,dcmplx,ncomp,nalms);

  sharp_job job;
  printf ("\niteration 0:\n");
  sharpd_build_job(&job,ALM2MAP,spin,0,&alm[0],&map[0],tinfo,alms,ntrans);
  sharp_execute_job(&job);
  printf("wall time for alm2map: %fs\n",job.time);
  printf("Performance: %fGFLOPs/s\n",1e-9*job.opcnt/job.time);

  map2alm_iter(tinfo,map,alm,alm2,lmax,mmax,npix,nalms,spin,ntrans,niter);

  DEALLOC2D(map);
  DEALLOC2D(alm);
  DEALLOC2D(alm2);

  sharp_destroy_alm_info(alms);
  }

int main(int argc, char **argv)
  {
#ifdef USE_MPI
  MPI_Init(NULL,NULL);
#endif
  module_startup_c("sharp_test",argc,7,
    "<healpix|ecp|gauss> <lmax> <nside|nphi> <niter> <spin> <ntrans>",1);

  int lmax=atoi(argv[2]);
  int niter=atoi(argv[4]);
  int spin=atoi(argv[5]);
  int ntrans=atoi(argv[6]);

  printf("Testing map analysis accuracy.\n");
  printf("lmax=%d, %d iterations, spin=%d\n", lmax, niter, spin);

  sharp_geom_info *tinfo;
  if (strcmp(argv[1],"gauss")==0)
    {
    int nrings=lmax+1;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting Gaussian grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"ecp")==0)
    {
    int nrings=2*lmax+2;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting ECP grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    sharp_make_ecp_geom_info (nrings, ppring, 0., 1, ppring, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"healpix")==0)
    {
    int nside=atoi(argv[3]);
    if (nside<1) nside=1;
    ptrdiff_t npix=12*(ptrdiff_t)nside*nside;
    printf("\nTesting Healpix grid (nside=%d, %ld pixels)\n",
          nside,(long)npix);
    sharp_make_healpix_geom_info (nside, 1, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else
    UTIL_FAIL("unknown grid geometry");

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
  }
