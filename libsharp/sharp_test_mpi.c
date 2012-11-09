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

/*! \file sharp_test_mpi.c
    Accuracy test for libsharp's map analysis with MPI support.

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

#ifdef USE_MPI

#include <stdio.h>
#include <string.h>
#include "sharp_mpi.h"
#include "sharp_geomhelpers.h"
#include "sharp_almhelpers.h"
#include "c_utils.h"
#include "walltime_c.h"
#include "sharp_announce.h"
#include "sharp_core.h"

typedef complex double dcmplx;

int ntasks, mytask;

static unsigned long long totalops (unsigned long long val)
  {
  unsigned long long tmp;
  MPI_Allreduce (&val, &tmp,1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  return tmp;
  }

static double maxTime (double val)
  {
  double tmp;
  MPI_Allreduce (&val, &tmp,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return tmp;
  }

static double drand (double min, double max)
  { return min + (max-min)*rand()/(RAND_MAX+1.0); }

static ptrdiff_t get_nalms(const sharp_alm_info *ainfo)
  {
  ptrdiff_t res=0;
  for (int i=0; i<ainfo->nm; ++i)
    res += ainfo->lmax-ainfo->mval[i]+1;
  return res;
  }

static ptrdiff_t get_npix(const sharp_geom_info *ginfo)
  {
  ptrdiff_t res=0;
  for (int i=0; i<ginfo->npairs; ++i)
    {
    res += ginfo->pair[i].r1.nph;
    if (ginfo->pair[i].r2.nph>0) res += ginfo->pair[i].r2.nph;
    }
  return res;
  }

static void reduce_alm_info(sharp_alm_info *ainfo)
  {
  int nmnew=0;
  ptrdiff_t ofs = 0;
  for (int i=mytask; i<ainfo->nm; i+=ntasks,++nmnew)
    {
    ainfo->mval[nmnew]=ainfo->mval[i];
    ainfo->mvstart[nmnew]=ofs-ainfo->mval[nmnew];
    ofs+=ainfo->lmax-ainfo->mval[nmnew]+1;
    }
  ainfo->nm=nmnew;
  }

static void reduce_geom_info(sharp_geom_info *ginfo)
  {
  int npairsnew=0;
  ptrdiff_t ofs = 0;
  for (int i=mytask; i<ginfo->npairs; i+=ntasks,++npairsnew)
    {
    ginfo->pair[npairsnew]=ginfo->pair[i];
    ginfo->pair[npairsnew].r1.ofs=ofs;
    ofs+=ginfo->pair[npairsnew].r1.nph;
    ginfo->pair[npairsnew].r2.ofs=ofs;
    if (ginfo->pair[npairsnew].r2.nph>0) ofs+=ginfo->pair[npairsnew].r2.nph;
    }
  ginfo->npairs=npairsnew;
  }

static void random_alm (dcmplx *alm, sharp_alm_info *helper, int spin)
  {
  static int cnt=0;
  ++cnt;
  for (int mi=0;mi<helper->nm; ++mi)
    {
    int m=helper->mval[mi];
    srand(1234567*cnt+8912*m);
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
  const sharp_alm_info *ainfo, int ncomp)
  {
  long nalms=get_nalms(ainfo), nalms_tot;
  MPI_Allreduce(&nalms,&nalms_tot,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  for (int i=0; i<ncomp; ++i)
    {
    double sum=0, sum2=0, maxdiff=0, sumtot, sum2tot, maxdifftot;
    for (int mi=0; mi<ainfo->nm; ++mi)
      {
      int m=ainfo->mval[mi];
      for (int l=m; l<=ainfo->lmax; ++l)
        {
        ptrdiff_t idx=sharp_alm_index(ainfo,l,mi);
        double x=creal(alm[i][idx])-creal(alm2[i][idx]),
               y=cimag(alm[i][idx])-cimag(alm2[i][idx]);
        sum+=x*x+y*y;
        sum2+=creal(alm[i][idx])*creal(alm[i][idx])
             +cimag(alm[i][idx])*cimag(alm[i][idx]);
        if (fabs(x)>maxdiff) maxdiff=fabs(x);
        if (fabs(y)>maxdiff) maxdiff=fabs(y);
        }
      }

    MPI_Allreduce(&sum,&sumtot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sum2,&sum2tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&maxdiff,&maxdifftot,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    sumtot=sqrt(sumtot/nalms_tot);
    sum2tot=sqrt(sum2tot/nalms_tot);
    if (mytask==0)
      printf("component %i: rms %e, maxerr %e\n",i, sumtot/sum2tot, maxdifftot);
    }
  }

static void map2alm_iter (sharp_geom_info *tinfo, double **map,
  dcmplx **alm_orig, dcmplx **alm, int lmax, int mmax,
  ptrdiff_t npix, int spin, int ntrans, int niter)
  {
  int ncomp = ntrans*((spin==0) ? 1 : 2);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);
  reduce_alm_info(alms);

  double jtime;
  unsigned long long jopcnt;

  sharp_execute_mpi(MPI_COMM_WORLD,SHARP_MAP2ALM,spin,&alm[0],&map[0],
    tinfo,alms,ntrans,SHARP_DP,&jtime,&jopcnt);
  unsigned long long opcnt=totalops(jopcnt);
  double timer=maxTime(jtime);
  if (mytask==0) printf("wall time for map2alm: %fs\n",timer);
  if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/timer);
  measure_errors(alm_orig,alm,alms,ncomp);

  for (int iter=0; iter<niter; ++iter)
    {
    double **map2;
    ALLOC2D(map2,double,ncomp,npix);
    if (mytask==0) printf ("\niteration %i:\n", iter+1);
    sharp_execute_mpi(MPI_COMM_WORLD,SHARP_ALM2MAP,spin,&alm[0],&map2[0],
      tinfo,alms,ntrans,SHARP_DP,&jtime,&jopcnt);
    opcnt=totalops(jopcnt);
    timer=maxTime(jtime);
    if (mytask==0) printf("wall time for alm2map: %fs\n",timer);
    if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/timer);
    for (int i=0; i<ncomp; ++i)
      for (ptrdiff_t m=0; m<npix; ++m)
        map2[i][m] = map[i][m]-map2[i][m];

    sharp_execute_mpi(MPI_COMM_WORLD,SHARP_MAP2ALM,spin,&alm[0],&map2[0],
      tinfo,alms,ntrans,SHARP_DP|SHARP_ADD,&jtime,&jopcnt);
    opcnt=totalops(jopcnt);
    timer=maxTime(jtime);
    if (mytask==0) printf("wall time for map2alm: %fs\n",wallTime()-timer);
    if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/timer);
    DEALLOC2D(map2);
    measure_errors(alm_orig,alm,alms,ncomp);
    }

  sharp_destroy_alm_info(alms);
  }

static void check_accuracy (sharp_geom_info *tinfo, ptrdiff_t lmax,
  ptrdiff_t mmax, ptrdiff_t npix, int spin, int ntrans, int niter)
  {
  int ncomp = ntrans*((spin==0) ? 1 : 2);

  double **map;
  ALLOC2D(map,double,ncomp,npix);

  double jtime;
  unsigned long long jopcnt;

  sharp_alm_info *alms;
  ptrdiff_t nalms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);
  reduce_alm_info(alms);
  nalms=get_nalms(alms);

  dcmplx **alm;
  ALLOC2D(alm,dcmplx,ncomp,nalms);
  srand(4);
  for (int i=0; i<ncomp; ++i)
    random_alm(alm[i],alms,spin);

  dcmplx **alm2;
  ALLOC2D(alm2,dcmplx,ncomp,nalms);

  if (mytask==0) printf ("\niteration 0:\n");
  sharp_execute_mpi(MPI_COMM_WORLD,SHARP_ALM2MAP,spin,&alm[0],&map[0],
    tinfo,alms,ntrans,SHARP_DP,&jtime,&jopcnt);
  unsigned long long opcnt=totalops(jopcnt);
  double timer=maxTime(jtime);
  if (mytask==0) printf("wall time for alm2map: %fs\n",timer);
  if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/timer);

  map2alm_iter(tinfo, map, alm, alm2, lmax, mmax, npix, spin, ntrans, niter);

  DEALLOC2D(map);
  DEALLOC2D(alm);
  DEALLOC2D(alm2);

  sharp_destroy_alm_info(alms);
  }

int main(int argc, char **argv)
  {
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mytask);

  sharp_module_startup("sharp_test_mpi",argc,7,
    "<healpix|ecp|gauss> <lmax> <nside|nphi> <niter> <spin> <ntrans>",
    mytask==0);
  int lmax=atoi(argv[2]);
  int niter=atoi(argv[4]);
  int spin=atoi(argv[5]);
  int ntrans=atoi(argv[6]);

  if (mytask==0)
    {
    printf("Testing map analysis accuracy.\n");
    printf("lmax=%d, %d iterations, spin=%d\n", lmax, niter, spin);
    }

  sharp_geom_info *tinfo;
  if (strcmp(argv[1],"gauss")==0)
    {
    int nrings=lmax+1;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    if (mytask==0)
      printf("\nTesting Gaussian grid (%d rings, %d pixels/ring, %ld pixels)\n",
             nrings,ppring,(long)npix);
    sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);
    reduce_geom_info(tinfo);
    npix=get_npix(tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"ecp")==0)
    {
    int nrings=2*lmax+2;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    if (mytask==0)
      printf("\nTesting ECP grid (%d rings, %d pixels/ring, %ld pixels)\n",
             nrings,ppring,(long)npix);
    sharp_make_ecp_geom_info (nrings, ppring, 0., 1, ppring, &tinfo);
    reduce_geom_info(tinfo);
    npix=get_npix(tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"healpix")==0)
    {
    int nside=atoi(argv[3]);
    if (nside<1) nside=1;
    ptrdiff_t npix=12*(ptrdiff_t)nside*nside;
    if (mytask==0)
      printf("\nTesting Healpix grid (nside=%d, %ld pixels)\n",
             nside,(long)npix);
    sharp_make_healpix_geom_info (nside, 1, &tinfo);
    reduce_geom_info(tinfo);
    npix=get_npix(tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,ntrans,niter);
    sharp_destroy_geom_info(tinfo);
    }
  else
    UTIL_FAIL("unknown grid geometry");

  MPI_Finalize();
  return 0;
  }

#else

#include "c_utils.h"

int main(void)
  { UTIL_FAIL("MPI support not enabled."); return 1; }

#endif
