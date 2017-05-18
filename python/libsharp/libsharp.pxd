cdef extern from "sharp.h":
    ctypedef long ptrdiff_t

    void sharp_legendre_transform_s(float *bl, float *recfac, ptrdiff_t lmax, float *x,
                                    float *out, ptrdiff_t nx)
    void sharp_legendre_transform(double *bl, double *recfac, ptrdiff_t lmax, double *x,
                                  double *out, ptrdiff_t nx)
    void sharp_legendre_transform_recfac(double *r, ptrdiff_t lmax)
    void sharp_legendre_transform_recfac_s(float *r, ptrdiff_t lmax)
    void sharp_legendre_roots(int n, double *x, double *w)

    # sharp_lowlevel.h
    ctypedef struct sharp_alm_info:
        pass

    ctypedef struct sharp_geom_info:
        pass

    void sharp_make_alm_info (int lmax, int mmax, int stride,
                             ptrdiff_t *mvstart, sharp_alm_info **alm_info)

    void sharp_make_geom_info (int nrings, int *nph, ptrdiff_t *ofs,
                               int *stride, double *phi0, double *theta,
                               double *wgt, sharp_geom_info **geom_info)

    void sharp_destroy_alm_info(sharp_alm_info *info)
    void sharp_destroy_geom_info(sharp_geom_info *info)

    ptrdiff_t sharp_map_size(sharp_geom_info *info)
    ptrdiff_t sharp_alm_count(sharp_alm_info *self)


    ctypedef enum sharp_jobtype:
        SHARP_YtW
        SHARP_Yt
        SHARP_WY
        SHARP_Y

    ctypedef enum:
        SHARP_DP
        SHARP_ADD

    void sharp_execute(sharp_jobtype type_,
                       int spin,
                       void *alm,
                       void *map,
                       sharp_geom_info *geom_info,
                       sharp_alm_info *alm_info,
                       int ntrans,
                       int flags,
                       double *time,
                       unsigned long long *opcnt) nogil

    ctypedef enum:
        SHARP_ERROR_NO_MPI

    int sharp_execute_mpi_maybe (void *pcomm, sharp_jobtype type, int spin,
        void *alm, void *map, sharp_geom_info *geom_info,
        sharp_alm_info *alm_info, int ntrans, int flags, double *time,
        unsigned long long *opcnt) nogil

    void sharp_normalized_associated_legendre_table(int m, int lmax, int ntheta,
        double *theta, int ncols, double *out) nogil


cdef extern from "sharp_geomhelpers.h":
    void sharp_make_subset_healpix_geom_info(
        int nside, int stride, int nrings,
        int *rings, double *weight, sharp_geom_info **geom_info)
    void sharp_make_gauss_geom_info(
        int nrings, int nphi, double phi0,
        int stride_lon, int stride_lat, sharp_geom_info **geom_info)

cdef extern from "sharp_almhelpers.h":
    void sharp_make_triangular_alm_info (int lmax, int mmax, int stride,
        sharp_alm_info **alm_info)
    void sharp_make_rectangular_alm_info (int lmax, int mmax, int stride,
        sharp_alm_info **alm_info)
    void sharp_make_mmajor_real_packed_alm_info (int lmax, int stride,
        int nm, const int *ms, sharp_alm_info **alm_info)

