/*
 *  This file is part of libfftpack.
 *
 *  libfftpack is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libfftpack is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libfftpack; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libfftpack is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2005-2016 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "trig_utils.h"
#include "fftpack.h"
#include "bluestein.h"

/* returns the sum of all prime factors of n */
size_t prime_factor_sum (size_t n)
  {
  size_t result=0,tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  size_t limit=(size_t)sqrt(n+0.01);
  for (size_t x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+=x;
    n=tmp;
    limit=(size_t)sqrt(n+0.01);
    }
  if (n>1) result+=n;

  return result;
  }

/* returns the smallest composite of 2, 3 and 5 which is >= n */
static size_t good_size(size_t n)
  {
  if (n<=6) return n;

  size_t bestfac=2*n;
  for (size_t f2=1; f2<bestfac; f2*=2)
    for (size_t f23=f2; f23<bestfac; f23*=3)
      for (size_t f235=f23; f235<bestfac; f235*=5)
        if (f235>=n) bestfac=f235;
  return bestfac;
  }

void bluestein_i (size_t n, double **tstorage, size_t *worksize)
  {
  size_t n2=good_size(n*2-1);
  *worksize=2+2*n+8*n2+16;
  *tstorage = RALLOC(double,2+2*n+8*n2+16);
  ((size_t *)(*tstorage))[0]=n2;
  double *bk  = *tstorage+2;
  double *bkf = *tstorage+2+2*n;
  double *work= *tstorage+2+2*(n+n2);

/* initialize b_k */
  double *tmp=RALLOC(double,4*n);
  sincos_2pibyn(2*n,2*n,&tmp[1],&tmp[0],2);
  bk[0] = 1;
  bk[1] = 0;

  size_t coeff=0;
  for (size_t m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    bk[2*m  ] = tmp[2*coeff  ];
    bk[2*m+1] = tmp[2*coeff+1];
    }
  DEALLOC(tmp);

/* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
  double xn2 = 1./n2;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (size_t m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
    }
  for (size_t m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;
  cffti (n2,work);
  cfftf (n2,bkf,work);
  }

void bluestein (size_t n, double *data, double *tstorage, int isign)
  {
  size_t n2=*((size_t *)tstorage);
  double *bk  = tstorage+2;
  double *bkf = tstorage+2+2*n;
  double *work= tstorage+2+2*(n+n2);
  double *akf = tstorage+2+2*n+6*n2+16;

/* initialize a_k and FFT it */
  if (isign>0)
    for (size_t m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   - data[m+1]*bk[m+1];
      akf[m+1] = data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  else
    for (size_t m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   + data[m+1]*bk[m+1];
      akf[m+1] =-data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  for (size_t m=2*n; m<2*n2; ++m)
    akf[m]=0;

  cfftf (n2,akf,work);

/* do the convolution */
  if (isign>0)
    for (size_t m=0; m<2*n2; m+=2)
      {
      double im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (size_t m=0; m<2*n2; m+=2)
      {
      double im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }


/* inverse FFT */
  cfftb (n2,akf,work);

/* multiply by b_k* */
  if (isign>0)
    for (size_t m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] - bk[m+1]*akf[m+1];
      data[m+1] = bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  else
    for (size_t m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] + bk[m+1]*akf[m+1];
      data[m+1] =-bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  }
