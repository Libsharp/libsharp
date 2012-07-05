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

/*! \file sharp.c
 *  Spherical transform library
 *
 *  Copyright (C) 2006-2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include "ls_fft.h"
#include "sharp_ylmgen_c.h"
#include "sharp.h"
#include "c_utils.h"
#include "sharp_core.h"
#include "sharp_vecutil.h"
#include "walltime_c.h"

typedef complex double dcmplx;
typedef complex float  fcmplx;

static void get_chunk_info (int ndata, int nmult, int *nchunks, int *chunksize)
  {
  static const int chunksize_min=500, nchunks_max=10;
  *chunksize = IMAX(chunksize_min,(ndata+nchunks_max-1)/nchunks_max);
  *chunksize = ((*chunksize+nmult-1)/nmult)*nmult;
  *nchunks = (ndata+*chunksize-1) / *chunksize;
  }

typedef struct
  {
  double s;
  int i;
  } idxhelper;

static int idx_compare (const void *xa, const void *xb)
  {
  const idxhelper *a=xa, *b=xb;
  return (a->s > b->s) ? -1 : (a->s < b->s) ? 1 : 0;
  }

typedef struct
  {
  double phi0_;
  dcmplx *shiftarr, *work;
  int s_shift, s_work;
  real_plan plan;
  int norot;
  } ringhelper;

static void ringhelper_init (ringhelper *self)
  {
  static ringhelper rh_null = { 0, NULL, NULL, 0, 0, NULL, 0 };
  *self = rh_null;
  }

static void ringhelper_destroy (ringhelper *self)
  {
  if (self->plan) kill_real_plan(self->plan);
  DEALLOC(self->shiftarr);
  DEALLOC(self->work);
  ringhelper_init(self);
  }

static void ringhelper_update (ringhelper *self, int nph, int mmax, double phi0)
  {
  self->norot = (fabs(phi0)<1e-14);
  if (!(self->norot))
    if ((mmax!=self->s_shift-1) || (!FAPPROX(phi0,self->phi0_,1e-12)))
      {
      RESIZE (self->shiftarr,dcmplx,mmax+1);
      self->s_shift = mmax+1;
      self->phi0_ = phi0;
      for (int m=0; m<=mmax; ++m)
        self->shiftarr[m] = cos(m*phi0) + _Complex_I*sin(m*phi0);
      }
  if (!self->plan) self->plan=make_real_plan(nph);
  if (nph!=(int)self->plan->length)
    {
    kill_real_plan(self->plan);
    self->plan=make_real_plan(nph);
    }
  GROW(self->work,dcmplx,self->s_work,nph);
  }

static int ringinfo_compare (const void *xa, const void *xb)
  {
  const sharp_ringinfo *a=xa, *b=xb;
  return (a->sth < b->sth) ? -1 : (a->sth > b->sth) ? 1 : 0;
  }
static int ringpair_compare (const void *xa, const void *xb)
  {
  const sharp_ringpair *a=xa, *b=xb;
  if (a->r1.nph==b->r1.nph)
    return (a->r1.phi0 < b->r1.phi0) ? -1 : (a->r1.phi0 > b->r1.phi0) ? 1 : 0;
  return (a->r1.nph<b->r1.nph) ? -1 : 1;
  }

void sharp_make_general_alm_info (int lmax, int nm, int stride, const int *mval,
  const ptrdiff_t *mstart, sharp_alm_info **alm_info)
  {
  sharp_alm_info *info = RALLOC(sharp_alm_info,1);
  info->lmax = lmax;
  info->nm = nm;
  info->mval = RALLOC(int,nm);
  info->mvstart = RALLOC(ptrdiff_t,nm);
  info->stride = stride;
  for (int mi=0; mi<nm; ++mi)
    {
    info->mval[mi] = mval[mi];
    info->mvstart[mi] = mstart[mi];
    }
  *alm_info = info;
  }

void sharp_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, sharp_alm_info **alm_info)
  {
  int *mval=RALLOC(int,mmax+1);
  for (int i=0; i<=mmax; ++i)
    mval[i]=i;
  sharp_make_general_alm_info (lmax, mmax+1, stride, mval, mstart, alm_info);
  DEALLOC(mval);
  }

ptrdiff_t sharp_alm_index (const sharp_alm_info *self, int l, int mi)
  { return self->mvstart[mi]+self->stride*l; }

void sharp_destroy_alm_info (sharp_alm_info *info)
  {
  DEALLOC (info->mval);
  DEALLOC (info->mvstart);
  DEALLOC (info);
  }

void sharp_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, sharp_geom_info **geom_info)
  {
  sharp_geom_info *info = RALLOC(sharp_geom_info,1);
  sharp_ringinfo *infos = RALLOC(sharp_ringinfo,nrings);

  int pos=0;
  info->pair=RALLOC(sharp_ringpair,nrings);
  info->npairs=0;
  *geom_info = info;

  for (int m=0; m<nrings; ++m)
    {
    infos[m].theta = theta[m];
    infos[m].cth = cos(theta[m]);
    infos[m].sth = sin(theta[m]);
    infos[m].weight = weight[m];
    infos[m].phi0 = phi0[m];
    infos[m].ofs = ofs[m];
    infos[m].stride = stride[m];
    infos[m].nph = nph[m];
    }
  qsort(infos,nrings,sizeof(sharp_ringinfo),ringinfo_compare);
  while (pos<nrings)
    {
    info->pair[info->npairs].r1=infos[pos];
    if ((pos<nrings-1) && FAPPROX(infos[pos].cth,-infos[pos+1].cth,1e-12))
      {
      info->pair[info->npairs].r2=infos[pos+1];
      ++pos;
      }
    else
      info->pair[info->npairs].r2.nph=-1;
    ++pos;
    ++info->npairs;
    }
  DEALLOC(infos);

  qsort(info->pair,info->npairs,sizeof(sharp_ringpair),ringpair_compare);
  }

void sharp_destroy_geom_info (sharp_geom_info *geom_info)
  {
  DEALLOC (geom_info->pair);
  DEALLOC (geom_info);
  }

static int sharp_get_mmax (int *mval, int nm)
  {
  int *mcheck=RALLOC(int,nm);
  SET_ARRAY(mcheck,0,nm,0);
  for (int i=0; i<nm; ++i)
    {
    int m_cur=mval[i];
    UTIL_ASSERT((m_cur>=0) && (m_cur<nm), "m out of range");
    UTIL_ASSERT(mcheck[m_cur]==0, "duplicate m value");
    mcheck[m_cur]=1;
    }
  DEALLOC(mcheck);
  return nm-1; // FIXME: this looks wrong
  }

static void ringhelper_phase2ring (ringhelper *self,
  const sharp_ringinfo *info, void *data, int mmax, const dcmplx *phase,
  int pstride, sharp_fde fde)
  {
  int nph = info->nph;
  int stride = info->stride;

  ringhelper_update (self, nph, mmax, info->phi0);
  self->work[0]=phase[0];
  SET_ARRAY(self->work,1,nph,0.);

#if 0
  if (self->norot)
    for (int m=1; m<=mmax; ++m)
      {
      int idx1 = m%nph;
      int idx2 = nph-1-((m-1)%nph);
      self->work[idx1]+=phase[m*pstride];
      self->work[idx2]+=conj(phase[m*pstride]);
      }
  else
    for (int m=1; m<=mmax; ++m)
      {
      int idx1 = m%nph;
      int idx2 = nph-1-((m-1)%nph);
      dcmplx tmp = phase[m*pstride]*self->shiftarr[m];
      self->work[idx1]+=tmp;
      self->work[idx2]+=conj(tmp);
      }
#else
  int idx1=1, idx2=nph-1;
  for (int m=1; m<=mmax; ++m)
    {
    dcmplx tmp = phase[m*pstride];
    if(!self->norot) tmp*=self->shiftarr[m];
    self->work[idx1]+=tmp;
    self->work[idx2]+=conj(tmp);
    if (++idx1>=nph) idx1=0;
    if (--idx2<0) idx2=nph-1;
    }
#endif
  real_plan_backward_c (self->plan, (double *)(self->work));
  if (fde==DOUBLE)
    for (int m=0; m<nph; ++m)
      ((double *)data)[m*stride+info->ofs] += creal(self->work[m]);
  else
    for (int m=0; m<nph; ++m)
      ((float *)data)[m*stride+info->ofs] += (float)creal(self->work[m]);
  }

static void ringhelper_ring2phase (ringhelper *self,
  const sharp_ringinfo *info, const void *data, int mmax, dcmplx *phase,
  int pstride, sharp_fde fde)
  {
  int nph = info->nph;
#if 1
  int maxidx = mmax; /* Enable this for traditional Healpix compatibility */
#else
  int maxidx = IMIN(nph-1,mmax);
#endif

  ringhelper_update (self, nph, mmax, -info->phi0);
  if (fde==DOUBLE)
    for (int m=0; m<nph; ++m)
      self->work[m] = ((double *)data)[info->ofs+m*info->stride]*info->weight;
  else
    for (int m=0; m<nph; ++m)
      self->work[m] = ((float *)data)[info->ofs+m*info->stride]*info->weight;

  real_plan_forward_c (self->plan, (double *)self->work);

  if (self->norot)
    for (int m=0; m<=maxidx; ++m)
      phase[m*pstride] = self->work[m%nph];
  else
    for (int m=0; m<=maxidx; ++m)
      phase[m*pstride]=self->work[m%nph]*self->shiftarr[m];

  for (int m=maxidx+1;m<=mmax; ++m)
    phase[m*pstride]=0.;
  }

static void ringhelper_pair2phase (ringhelper *self, int mmax,
  const sharp_ringpair *pair, const void *data, dcmplx *phase1, dcmplx *phase2,
  int pstride, sharp_fde fde)
  {
  ringhelper_ring2phase (self, &(pair->r1), data, mmax, phase1, pstride, fde);
  if (pair->r2.nph>0)
    ringhelper_ring2phase (self, &(pair->r2), data, mmax, phase2, pstride, fde);
  }

static void ringhelper_phase2pair (ringhelper *self, int mmax,
  const dcmplx *phase1, const dcmplx *phase2, int pstride,
  const sharp_ringpair *pair, void *data, sharp_fde fde)
  {
  ringhelper_phase2ring (self, &(pair->r1), data, mmax, phase1, pstride, fde);
  if (pair->r2.nph>0)
    ringhelper_phase2ring (self, &(pair->r2), data, mmax, phase2, pstride, fde);
  }

static void fill_map (const sharp_geom_info *ginfo, void *map, double value,
  sharp_fde fde)
  {
  for (int j=0;j<ginfo->npairs;++j)
    {
    if (fde==DOUBLE)
      {
      for (int i=0;i<ginfo->pair[j].r1.nph;++i)
        ((double *)map)[ginfo->pair[j].r1.ofs+i*ginfo->pair[j].r1.stride]=value;
      for (int i=0;i<ginfo->pair[j].r2.nph;++i)
        ((double *)map)[ginfo->pair[j].r2.ofs+i*ginfo->pair[j].r2.stride]=value;
      }
    else
      {
      for (int i=0;i<ginfo->pair[j].r1.nph;++i)
        ((float *)map)[ginfo->pair[j].r1.ofs+i*ginfo->pair[j].r1.stride]
          =(float)value;
      for (int i=0;i<ginfo->pair[j].r2.nph;++i)
        ((float *)map)[ginfo->pair[j].r2.ofs+i*ginfo->pair[j].r2.stride]
          =(float)value;
      }
    }
  }

static void fill_alm (const sharp_alm_info *ainfo, void *alm, dcmplx value,
  sharp_fde fde)
  {
  if (fde==DOUBLE)
    for (int mi=0;mi<ainfo->nm;++mi)
      for (int l=ainfo->mval[mi];l<=ainfo->lmax;++l)
        ((dcmplx *)alm)[sharp_alm_index(ainfo,l,mi)] = value;
  else
    for (int mi=0;mi<ainfo->nm;++mi)
      for (int l=ainfo->mval[mi];l<=ainfo->lmax;++l)
        ((fcmplx *)alm)[sharp_alm_index(ainfo,l,mi)] = (fcmplx)value;
  }

static void init_output (sharp_job *job)
  {
  if (job->add_output) return;
  if (job->type == MAP2ALM)
    for (int i=0; i<job->ntrans*job->nalm; ++i)
      fill_alm (job->ainfo,job->alm[i],0.,job->fde);
  else
    for (int i=0; i<job->ntrans*job->nmaps; ++i)
      fill_map (job->ginfo,job->map[i],0.,job->fde);
  }

static void alloc_phase (sharp_job *job, int nm, int ntheta)
  { job->phase=RALLOC(dcmplx,2*job->ntrans*job->nmaps*nm*ntheta); }

static void dealloc_phase (sharp_job *job)
  { DEALLOC(job->phase); }

//FIXME: set phase to zero if not MAP2ALM?
static void map2phase (sharp_job *job, int mmax, int llim, int ulim)
  {
  if (job->type != MAP2ALM) return;
  int pstride = 2*job->ntrans*job->nmaps;
#pragma omp parallel
{
  ringhelper helper;
  ringhelper_init(&helper);
#pragma omp for schedule(dynamic,1)
  for (int ith=llim; ith<ulim; ++ith)
    {
    int dim2 = pstride*(ith-llim)*(mmax+1);
    for (int i=0; i<job->ntrans*job->nmaps; ++i)
      ringhelper_pair2phase(&helper,mmax,&job->ginfo->pair[ith], job->map[i],
        &job->phase[dim2+2*i], &job->phase[dim2+2*i+1], pstride, job->fde);
    }
  ringhelper_destroy(&helper);
} /* end of parallel region */
  }

static void alloc_almtmp (sharp_job *job, int lmax)
  { job->almtmp=RALLOC(dcmplx,job->ntrans*job->nalm*(lmax+1)); }

static void dealloc_almtmp (sharp_job *job)
  { DEALLOC(job->almtmp); }

static void alm2almtmp (sharp_job *job, int lmax, int mi)
  {
  if (job->type!=MAP2ALM)
    for (int l=job->ainfo->mval[mi]; l<=lmax; ++l)
      {
      ptrdiff_t aidx = sharp_alm_index(job->ainfo,l,mi);
      double fct = (job->type==ALM2MAP) ? job->norm_l[l] :
                    fabs(job->norm_l[l])*sqrt(l*(l+1.));
      for (int i=0; i<job->ntrans*job->nalm; ++i)
        if (job->fde==DOUBLE)
          job->almtmp[job->ntrans*job->nalm*l+i]
            = ((dcmplx *)job->alm[i])[aidx]*fct;
        else
          job->almtmp[job->ntrans*job->nalm*l+i]
            = ((fcmplx *)job->alm[i])[aidx]*fct;
      }
  else
    SET_ARRAY(job->almtmp,job->ntrans*job->nalm*job->ainfo->mval[mi],
              job->ntrans*job->nalm*(lmax+1),0.);
  }

static void almtmp2alm (sharp_job *job, int lmax, int mi)
  {
  if (job->type != MAP2ALM) return;
  for (int l=job->ainfo->mval[mi]; l<=lmax; ++l)
    {
    ptrdiff_t aidx = sharp_alm_index(job->ainfo,l,mi);
    for (int i=0;i<job->ntrans*job->nalm;++i)
      if (job->fde==DOUBLE)
        ((dcmplx *)job->alm[i])[aidx] +=
          job->almtmp[job->ntrans*job->nalm*l+i]*job->norm_l[l];
      else
        ((fcmplx *)job->alm[i])[aidx] +=
          (fcmplx)(job->almtmp[job->ntrans*job->nalm*l+i]*job->norm_l[l]);
    }
  }

static void phase2map (sharp_job *job, int mmax, int llim, int ulim)
  {
  if (job->type == MAP2ALM) return;
  int pstride = 2*job->ntrans*job->nmaps;
#pragma omp parallel
{
  ringhelper helper;
  ringhelper_init(&helper);
#pragma omp for schedule(dynamic,1)
  for (int ith=llim; ith<ulim; ++ith)
    {
    int dim2 = pstride*(ith-llim)*(mmax+1);
    for (int i=0; i<job->ntrans*job->nmaps; ++i)
      ringhelper_phase2pair(&helper,mmax,&job->phase[dim2+2*i],
        &job->phase[dim2+2*i+1],pstride,&job->ginfo->pair[ith],job->map[i],
        job->fde);
    }
  ringhelper_destroy(&helper);
} /* end of parallel region */
  }

void sharp_execute_job (sharp_job *job)
  {
  double timer=wallTime();
  job->opcnt=0;
  int lmax = job->ainfo->lmax,
      mmax=sharp_get_mmax(job->ainfo->mval, job->ainfo->nm);

  job->norm_l = Ylmgen_get_norm (lmax, job->spin);

/* clear output arrays if requested */
  init_output (job);

  int nchunks, chunksize;
  get_chunk_info(job->ginfo->npairs,job->nv*VLEN,&nchunks,&chunksize);
  alloc_phase (job,mmax+1,chunksize);

/* chunk loop */
  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=IMIN(llim+chunksize,job->ginfo->npairs);
    int *ispair = RALLOC(int,ulim-llim);
    double *cth = RALLOC(double,ulim-llim), *sth = RALLOC(double,ulim-llim);
    idxhelper *stmp = RALLOC(idxhelper,ulim-llim);
    for (int i=0; i<ulim-llim; ++i)
      {
      ispair[i] = job->ginfo->pair[i+llim].r2.nph>0;
      cth[i] = job->ginfo->pair[i+llim].r1.cth;
      sth[i] = job->ginfo->pair[i+llim].r1.sth;
      stmp[i].s=sth[i];
      stmp[i].i=i;
      }
    qsort (stmp,ulim-llim,sizeof(idxhelper),idx_compare);
    int *idx = RALLOC(int,ulim-llim);
    for (int i=0; i<ulim-llim; ++i)
      idx[i]=stmp[i].i;
    DEALLOC(stmp);

/* map->phase where necessary */
    map2phase (job, mmax, llim, ulim);

#pragma omp parallel
{
    sharp_job ljob = *job;
    ljob.opcnt=0;
    Ylmgen_C generator;
    Ylmgen_init (&generator,lmax,mmax,ljob.spin);
    alloc_almtmp(&ljob,lmax);

#pragma omp for schedule(dynamic,1)
    for (int mi=0; mi<job->ainfo->nm; ++mi)
      {
/* alm->alm_tmp where necessary */
      alm2almtmp (&ljob, lmax, mi);

      inner_loop (&ljob, ispair, cth, sth, llim, ulim, &generator, mi, idx);

/* alm_tmp->alm where necessary */
      almtmp2alm (&ljob, lmax, mi);
      }

    Ylmgen_destroy(&generator);
    dealloc_almtmp(&ljob);

#pragma omp critical
    job->opcnt+=ljob.opcnt;
} /* end of parallel region */

/* phase->map where necessary */
    phase2map (job, mmax, llim, ulim);

    DEALLOC(ispair);
    DEALLOC(cth);
    DEALLOC(sth);
    DEALLOC(idx);
    } /* end of chunk loop */

  DEALLOC(job->norm_l);
  dealloc_phase (job);
  job->time=wallTime()-timer;
  }

static void sharp_build_job_common (sharp_job *job, sharp_jobtype type, int spin,
  int add_output, const sharp_geom_info *geom_info,
  const sharp_alm_info *alm_info, int ntrans)
  {
  UTIL_ASSERT((ntrans>0),"bad number of simultaneous transforms");
  if (type==ALM2MAP_DERIV1) spin=1;
  UTIL_ASSERT((spin>=0)&&(spin<=30), "bad spin");
  job->type = type;
  job->spin = spin;
  job->norm_l = NULL;
  job->add_output = add_output;
  job->nmaps = (type==ALM2MAP_DERIV1) ? 2 : ((spin>0) ? 2 : 1);
  job->nalm = (type==ALM2MAP_DERIV1) ? 1 : ((spin>0) ? 2 : 1);
  job->ginfo = geom_info;
  job->ainfo = alm_info;
  job->nv = sharp_nv_oracle (type, spin, ntrans);
  job->time = 0.;
  job->opcnt = 0;
  job->ntrans = ntrans;
  }

void sharpd_build_job (sharp_job *job, sharp_jobtype type, int spin,
  int add_output, dcmplx **alm, double **map, const sharp_geom_info *geom_info,
  const sharp_alm_info *alm_info, int ntrans)
  {
  sharp_build_job_common (job, type, spin, add_output, geom_info, alm_info,
    ntrans);
  job->alm=(void **)alm;
  job->map=(void **)map;
  job->fde=DOUBLE;
  }

void sharps_build_job (sharp_job *job, sharp_jobtype type, int spin,
  int add_output, fcmplx **alm, float **map, const sharp_geom_info *geom_info,
  const sharp_alm_info *alm_info, int ntrans)
  {
  sharp_build_job_common (job, type, spin, add_output, geom_info, alm_info,
    ntrans);
  job->alm=(void **)alm;
  job->map=(void **)map;
  job->fde=FLOAT;
  }

int sharp_get_nv_max (void)
{ return 6; }

int sharp_nv_oracle (sharp_jobtype type, int spin, int ntrans)
  {
  if (type==ALM2MAP_DERIV1) spin=1;
  UTIL_ASSERT((ntrans>0),"bad number of simultaneous transforms");
  UTIL_ASSERT((spin>=0)&&(spin<=30), "bad spin");
#include "oracle.inc"

  return nv_opt[IMIN(ntrans,maxtr)-1][spin!=0][type];
  }

#include "sharp_mpi.c"
