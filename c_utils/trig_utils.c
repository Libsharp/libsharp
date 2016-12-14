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

/*! \file trig_utils.c
 *
 *  Copyright (C) 2016 Max-Planck-Society
 *  \author Martin Reinecke
 *
 *  Many inspirations for this code come from Tasche and Zeuner: "Improved
 *  Roundoff Error Analysis for Precomputed Twiddle Factors", Journal of
 *  Computational Analysis and Applications, 4, 2002.
 */

#include <math.h>
#include "c_utils.h"
#include "trig_utils.h"

void sincos_multi (size_t n, double alpha, double beta, double *s, double *c,
  int stride)
  {
  if (beta!=0.)
    { s[0]=sin(beta); c[0]=cos(beta); }
  else
    { s[0] = 0.; c[0]=1.; }
  for (size_t vl=1;;vl<<=1)
    {
    double sa=sin(alpha), ca=cos(alpha);
    size_t ilim = (2*vl>n) ? n : 2*vl;
    int d=vl*stride;
    for (int i=vl*stride; i<(int)(ilim*stride); i+=stride)
      { s[i] = s[i-d]*ca + c[i-d]*sa; c[i] = c[i-d]*ca - s[i-d]*sa; }
    if (ilim==n) return;
    alpha*=2;
    }
  }

void sincos_2pibyn (size_t n, size_t nang, double *s, double *c, int stride)
  {
  static const double twopi=6.28318530717958647692;
  // nmax: number of sin/cos pairs that must be genuinely computed; the rest
  // can be obtained via symmetries
  size_t nmax = ((n&3)==0) ? n/8+1 : ( ((n&1)==0) ? n/4+1 : n/2+1 );
  size_t ngoal = (nang<nmax) ? nang : nmax;
  sincos_multi (ngoal, twopi/n, 0., s, c, stride);
  size_t ndone=ngoal;
  if (ndone==nang) return;
  if ((n&3)==0)
    {
    ngoal=n/4+1;
    if (nang<ngoal) ngoal=nang;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      s[i*stride]=c[(n/4-i)*stride];
      c[i*stride]=s[(n/4-i)*stride];
      }
    ndone=ngoal;
    if (ngoal==nang) return;
    }
  if ((n&1)==0)
    {
    ngoal=n/2+1;
    if (nang<ngoal) ngoal=nang;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      c[i*stride]=-c[(n/2-i)*stride];
      s[i*stride]= s[(n/2-i)*stride];
      }
    ndone=ngoal;
    if (ngoal==nang) return;
    }
  ngoal=n;
  if (nang<ngoal) ngoal=nang;
  for (size_t i=ndone; i<ngoal; ++i)
    {
    c[i*stride]= c[(n-i)*stride];
    s[i*stride]=-s[(n-i)*stride];
    }
  ndone=ngoal;
  if (ngoal==nang) return;
  for (size_t i=ndone; i<nang; ++i)
    {
    c[i*stride]= c[(i-n)*stride];
    s[i*stride]= s[(i-n)*stride];
    }
  }


/* Code for sin/cos(2*pi*m/n). Taken from FFTW. */
static void mysincos(int m, int n, double *s, double *c)
  {
  static const double twopi=6.28318530717958647692;
  double theta, t;
  unsigned octant = 0;
  int quarter_n = n;

  n += n; n += n;
  m += m; m += m;

  if (m < 0) m += n;
  if (m > n - m) { m = n - m; octant |= 4; }
  if (m - quarter_n > 0) { m = m - quarter_n; octant |= 2; }
  if (m > quarter_n - m) { m = quarter_n - m; octant |= 1; }

  theta = (twopi*m)/n;
  *c = cos(theta); *s = sin(theta);

  if (octant & 1) { t = *c; *c = *s; *s = t; }
  if (octant & 2) { t = *c; *c = -(*s); *s = t; }
  if (octant & 4) { *s = -(*s); }
  }

void trigtest(void);
void trigtest(void)
  {
#define LENGTH 1234
  static const double twopi=6.28318530717958647692;
  double sc[2*LENGTH+34];
  for (int i=1; i<LENGTH; ++i)
    {
    sc[0]=sc[1]=sc[2*i+30+2]=sc[2*i+30+3]=10;
    sincos_2pibyn(i,i+15,&sc[2],&sc[3],2);
    UTIL_ASSERT(fabs(sc[0]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[1]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[2*i+30+2]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[2*i+30+3]-10.)<1e-16,"bad memory access");
    for (int j=0; j<i; ++j)
      {
      double c, s;
      mysincos(j,i,&s,&c);
      UTIL_ASSERT(fabs(sc[2*j+2]-s)<8e-16,"bad sin");
      UTIL_ASSERT(fabs(sc[2*j+3]-c)<8e-16,"bad cos");
      }
    sc[0]=sc[1]=sc[2*i+2]=sc[2*i+3]=10;
    sincos_multi(i,twopi*1.1/i,0.1,&sc[2],&sc[3],2);
    UTIL_ASSERT(fabs(sc[0]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[1]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[2*i+2]-10.)<1e-16,"bad memory access");
    UTIL_ASSERT(fabs(sc[2*i+3]-10.)<1e-16,"bad memory access");
    for (int j=0; j<i; ++j)
      {
      double ang=twopi*1.1/i*j+0.1;
      UTIL_ASSERT(fabs(sc[2*j+2]-sin(ang))<1e-15,"bad sin");
      UTIL_ASSERT(fabs(sc[2*j+3]-cos(ang))<1e-15,"bad cos");
      }
    }
  }
