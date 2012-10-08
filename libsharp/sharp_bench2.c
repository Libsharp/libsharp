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

/*! \file sharp_bench2.c
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

int main(int argc, char **argv)
  {
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mytask);
  int master=(mytask==0);

  sharp_module_startup("sharp_bench2",argc,7,
    "<healpix|ecp|gauss> <lmax> <nside|nphi> <a2m/m2a> <spin> <ntrans>",0);

  int lmax=atoi(argv[2]);
  sharp_jobtype jtype = (strcmp(argv[4],"a2m")==0) ?
    SHARP_ALM2MAP : SHARP_MAP2ALM;
  int spin=atoi(argv[5]);
  int ntrans=atoi(argv[6]);

  sharp_geom_info *tinfo;
  ptrdiff_t npix=0;
  int geom2=0;
  if (strcmp(argv[1],"gauss")==0)
    {
    int nrings=lmax+1;
    int ppring=geom2=atoi(argv[3]);
    sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);
    }
  else if (strcmp(argv[1],"ecp")==0)
    {
    int nrings=2*lmax+2;
    int ppring=geom2=atoi(argv[3]);
    sharp_make_ecp_geom_info (nrings, ppring, 0., 1, ppring, &tinfo);
    }
  else if (strcmp(argv[1],"healpix")==0)
    {
    int nside=atoi(argv[3]);
    if (nside<1) nside=1;
    geom2=nside;
    sharp_make_healpix_geom_info (nside, 1, &tinfo);
    }
  else
    UTIL_FAIL("unknown grid geometry");

  reduce_geom_info(tinfo);
  npix=get_npix(tinfo);

  int mmax=lmax;
  int ncomp = ntrans*((spin==0) ? 1 : 2);

  double **map;
  ALLOC2D(map,double,ncomp,npix);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);

  reduce_alm_info(alms);
  ptrdiff_t nalms=get_nalms(alms);

  dcmplx **alm;
  ALLOC2D(alm,dcmplx,ncomp,nalms);

  for (int n=0; n<ncomp; ++n)
    {
    for (int i=0; i<npix; ++i) map[n][i]=1;
    for (int i=0; i<nalms; ++i) alm[n][i]=1;
    }

  double time=1e20;
  unsigned long long opcnt=0;
  for (int ntries=0; (ntries<2)||(ntries*time<5); ++ntries)
    {
    double ltime;
    unsigned long long lopcnt;
    sharp_execute_mpi(MPI_COMM_WORLD,jtype,spin,0,&alm[0],&map[0],
      tinfo,alms,ntrans,1,0,&ltime,&lopcnt);

    ltime=maxTime(ltime);
    if (ltime<time) { time=ltime; opcnt=totalops(lopcnt); }
    }
  DEALLOC2D(map);
  DEALLOC2D(alm);

  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  double mHWM=totalMem(VmHWM());

  int nomp=omp_get_max_threads();

  if (master)
    printf("%-12s %-7s %-3s %2d %d %2d %3d %5d %5d %1d %.2e %7.2f %9.2f\n",
      getenv("HOST"),argv[1],argv[4],spin,VLEN,nomp,ntasks,lmax,geom2,ntrans,
      time,opcnt/(time*1e9),mHWM/(1<<20));

  MPI_Finalize();
  return 0;
  }

#else

#include "c_utils.h"

int main(void)
  { UTIL_FAIL("Need OpenMP and MPI"); return 1; }

#endif
