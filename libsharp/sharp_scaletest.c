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

/*! \file sharp_scaletest.c
    Copyright (C) 2012 Max-Planck-Society
    \author Martin Reinecke
*/

#include <stdio.h>

#if (defined(_OPENMP) && defined(USE_MPI))

#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include "sharp_mpi.h"
#include "sharp.h"
#include "sharp_vecutil.h"
#include "sharp_geomhelpers.h"
#include "sharp_almhelpers.h"
#include "c_utils.h"
#include "sharp_announce.h"
#include "sharp_core.h"
#include "memusage.h"

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

static double totalMem (double val)
  {
  double tmp;
  MPI_Allreduce (&val, &tmp,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return tmp;
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

static double drand (double min, double max)
  { return min + (max-min)*rand()/(RAND_MAX+1.0); }

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

static void measure_errors (dcmplx **alm, double *sqsum,
  const sharp_alm_info *ainfo, int ncomp)
  {
  long nalms=get_nalms(ainfo), nalms_tot;
  MPI_Allreduce(&nalms,&nalms_tot,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  for (int i=0; i<ncomp; ++i)
    {
    double sum=0, maxdiff=0, sumtot, sqsumtot, maxdifftot;
    for (int mi=0; mi<ainfo->nm; ++mi)
      {
      int m=ainfo->mval[mi];
      for (int l=m; l<=ainfo->lmax; ++l)
        {
        ptrdiff_t idx=sharp_alm_index(ainfo,l,mi);
        sum+=creal(alm[i][idx])*creal(alm[i][idx])
            +cimag(alm[i][idx])*cimag(alm[i][idx]);
        if (fabs(creal(alm[i][idx]))>maxdiff) maxdiff=fabs(creal(alm[i][idx]));
        if (fabs(cimag(alm[i][idx]))>maxdiff) maxdiff=fabs(cimag(alm[i][idx]));
        }
      }

    MPI_Allreduce(&sum,&sumtot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sqsum[i],&sqsumtot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&maxdiff,&maxdifftot,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    sumtot=sqrt(sumtot/nalms_tot);
    sqsumtot=sqrt(sqsumtot/nalms_tot);
    if (mytask==0)
      printf("component %i: rms %e, maxerr %e\n",i,sumtot/sqsumtot, maxdifftot);
    }
  }

int main(int argc, char **argv)
  {
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mytask);

  sharp_module_startup("sharp_scaletest",argc,4,"<lmax> <nphi> <spin>",
    mytask==0);

  int lmax=atoi(argv[1]);
  int spin=atoi(argv[3]);

  sharp_geom_info *tinfo;
  ptrdiff_t npix=0;
  int nrings=lmax+1;
  int ppring=atoi(argv[2]);
  sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);

  reduce_geom_info(tinfo);
  npix=get_npix(tinfo);

  int mmax=lmax;
  int ncomp = (spin==0) ? 1 : 2;

  double **map;
  ALLOC2D(map,double,ncomp,npix);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);

  reduce_alm_info(alms);
  ptrdiff_t nalms=get_nalms(alms);

  dcmplx **alm;
  ALLOC2D(alm,dcmplx,ncomp,nalms);

  if (mytask==0) printf("Testing map analysis accuracy.\n");
  if (mytask==0) printf("lmax=%d, spin=%d\n", lmax, spin);
  if (mytask==0)
    printf("\nTesting Gaussian grid (%d rings, %d pixels/ring, %ld pixels)\n",
      nrings,ppring,(long)npix);

  for (int n=0; n<ncomp; ++n)
    random_alm(alm[n],alms,spin);

  double time=0;
  unsigned long long opcnt=0;
  sharp_execute_mpi(MPI_COMM_WORLD,SHARP_ALM2MAP,spin,&alm[0],&map[0],
    tinfo,alms,1,SHARP_DP,&time,&opcnt);

  time=maxTime(time);
  opcnt=totalops(opcnt);
  if (mytask==0) printf("wall time for alm2map: %fs\n",time);
  if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/time);

  double *sqsum=RALLOC(double,ncomp);
  for (int i=0; i<ncomp; ++i)
    {
    sqsum[i]=0;
    for (ptrdiff_t j=0; j<nalms; ++j)
      {
      sqsum[i]+=creal(alm[i][j])*creal(alm[i][j])
               +cimag(alm[i][j])*cimag(alm[i][j]);
      alm[i][j]=-alm[i][j];
      }
    }

  sharp_execute_mpi(MPI_COMM_WORLD,SHARP_MAP2ALM,spin,&alm[0],&map[0],
    tinfo,alms,1,SHARP_DP|SHARP_ADD,&time,&opcnt);

  time=maxTime(time);
  opcnt=totalops(opcnt);
  if (mytask==0) printf("wall time for map2alm: %fs\n",time);
  if (mytask==0) printf("Performance: %fGFLOPs/s\n",1e-9*opcnt/time);

  measure_errors(alm,sqsum,alms,ncomp);

  DEALLOC2D(map);
  DEALLOC2D(alm);

  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  double mHWM=totalMem(VmHWM());
  if (mytask==0) printf("\nMemory high water mark: %.2f MB\n",mHWM/(1<<20));

  MPI_Finalize();
  return 0;
  }

#else

#include "c_utils.h"

int main(void)
  { UTIL_FAIL("Need OpenMP and MPI"); return 1; }

#endif
