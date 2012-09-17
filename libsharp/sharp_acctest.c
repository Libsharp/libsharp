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

/*! \file sharp_acctest.c
    Systematic accuracy test for libsharp.

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
#include "sharp_announce.h"
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
    UTIL_ASSERT((maxdiff<1e-10)&&(sum/sum2<1e-10),"error");
    }
  }

static void check_sign_scale(void)
  {
  int lmax=50;
  int mmax=lmax;
  sharp_geom_info *tinfo;
  int nrings=lmax+1;
  int ppring=2*lmax+2;
  ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
  sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);

  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(lmax,mmax,1,&alms);
  ptrdiff_t nalms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);

  for (int ntrans=1; ntrans<10; ++ntrans)
    {
    double **map;
    ALLOC2D(map,double,2*ntrans,npix);

    dcmplx **alm;
    ALLOC2D(alm,dcmplx,2*ntrans,nalms);
    for (int i=0; i<2*ntrans; ++i)
      for (int j=0; j<nalms; ++j)
        alm[i][j]=1.+_Complex_I;

    sharp_execute(SHARP_ALM2MAP,0,0,(void **)&alm[0],(void **)&map[0],tinfo,
      alms,ntrans,1,0,NULL,NULL);
    for (int it=0; it<ntrans; ++it)
      {
      UTIL_ASSERT(FAPPROX(map[it][0     ], 3.588246976618616912e+00,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[it][npix/2], 4.042209792157496651e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[it][npix-1],-1.234675107554816442e+01,1e-12),
        "error");
      }
    sharp_execute(SHARP_ALM2MAP,1,0,(void **)&alm[0],(void **)&map[0],tinfo,
      alms,ntrans,1,0,NULL,NULL);
    for (int it=0; it<ntrans; ++it)
      {
      UTIL_ASSERT(FAPPROX(map[2*it  ][0     ], 2.750897760535633285e+00,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix/2], 3.137704477368562905e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix-1],-8.405730859837063917e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][0     ],-2.398026536095463346e+00,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix/2],-4.961140548331700728e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix-1],-1.412765834230440021e+01,1e-12),
        "error");
      }

    sharp_execute(SHARP_ALM2MAP,2,0,(void **)&alm[0],(void **)&map[0],tinfo,
      alms,ntrans,1,0,NULL,NULL);
    for (int it=0; it<ntrans; ++it)
      {
      UTIL_ASSERT(FAPPROX(map[2*it  ][0     ],-1.398186224727334448e+00,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix/2],-2.456676000884031197e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix-1],-1.516249174408820863e+02,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][0     ],-3.173406200299964119e+00,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix/2],-5.831327404513146462e+01,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix-1],-1.863257892248353897e+01,1e-12),
        "error");
      }

    sharp_execute(SHARP_ALM2MAP_DERIV1,1,0,(void **)&alm[0],(void **)&map[0],
      tinfo,alms,ntrans,1,0,NULL,NULL);
    for (int it=0; it<ntrans; ++it)
      {
      UTIL_ASSERT(FAPPROX(map[2*it  ][0     ],-6.859393905369091105e-01,1e-11),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix/2],-2.103947835973212364e+02,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it  ][npix-1],-1.092463246472086439e+03,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][0     ],-1.411433220713928165e+02,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix/2],-1.146122859381925082e+03,1e-12),
        "error");
      UTIL_ASSERT(FAPPROX(map[2*it+1][npix-1], 7.821618677689795049e+02,1e-12),
        "error");
      }

    DEALLOC2D(map);
    DEALLOC2D(alm);
    }

  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);
  }

static void check_accuracy (sharp_geom_info *tinfo, ptrdiff_t lmax,
  ptrdiff_t mmax, ptrdiff_t npix, int spin, int ntrans, int nv)
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

  sharp_execute(SHARP_ALM2MAP,spin,0,(void **)(&alm[0]),(void **)(&map[0]),
    tinfo,alms,ntrans,1,nv,NULL,NULL);
  sharp_execute(SHARP_MAP2ALM,spin,0,(void **)(&alm2[0]),(void **)(&map[0]),
    tinfo,alms,ntrans,1,nv,NULL,NULL);
  measure_errors(alm,alm2,nalms,ncomp);

  DEALLOC2D(map);
  DEALLOC2D(alm);
  DEALLOC2D(alm2);

  sharp_destroy_alm_info(alms);
  }

int main(void)
  {
#ifdef USE_MPI
  MPI_Init(NULL,NULL);
#endif
  sharp_module_startup("sharp_acctest",1,1,"",1);

  int lmax=127;

  printf("Checking signs and scales.\n");
  check_sign_scale();
  printf("Passed.\n\n");

  printf("Testing map analysis accuracy.\n");

  sharp_geom_info *tinfo;
  int nrings=lmax+1;
  int ppring=2*lmax+2;
  ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
  sharp_make_gauss_geom_info (nrings, ppring, 1, ppring, &tinfo);
  for (int nv=1; nv<=6; ++nv)
    for (int ntrans=1; ntrans<=6; ++ntrans)
      {
      check_accuracy(tinfo,lmax,lmax,npix,0,ntrans,nv);
      check_accuracy(tinfo,lmax,lmax,npix,1,ntrans,nv);
      check_accuracy(tinfo,lmax,lmax,npix,2,ntrans,nv);
      check_accuracy(tinfo,lmax,lmax,npix,3,ntrans,nv);
      check_accuracy(tinfo,lmax,lmax,npix,30,ntrans,nv);
      }
  sharp_destroy_geom_info(tinfo);
  printf("Passed.\n\n");

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
  }
