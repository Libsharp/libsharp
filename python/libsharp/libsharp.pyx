import numpy as np
cimport numpy as np
cimport cython

__all__ = ['legendre_transform', 'legendre_roots', 'sht', 'synthesis', 'adjoint_synthesis',
           'analysis', 'adjoint_analysis', 'healpix_grid', 'triangular_order', 'rectangular_order',
           'packed_real_order', 'normalized_associated_legendre_table']


def legendre_transform(x, bl, out=None):
    if out is None:
        out = np.empty_like(x)
    if out.shape[0] == 0:
        return out
    elif x.dtype == np.float64:
        if bl.dtype != np.float64:
            bl = bl.astype(np.float64)
        return _legendre_transform(x, bl, out=out)
    elif x.dtype == np.float32:
        if bl.dtype != np.float32:
            bl = bl.astype(np.float32)
        return _legendre_transform_s(x, bl, out=out)
    else:
        raise ValueError("unsupported dtype")


def _legendre_transform(double[::1] x, double[::1] bl, double[::1] out):
    if out.shape[0] != x.shape[0]:
        raise ValueError('x and out must have same shape')
    sharp_legendre_transform(&bl[0], NULL, bl.shape[0] - 1, &x[0], &out[0], x.shape[0])
    return np.asarray(out)


def _legendre_transform_s(float[::1] x, float[::1] bl, float[::1] out):
    if out.shape[0] != x.shape[0]:
        raise ValueError('x and out must have same shape')
    sharp_legendre_transform_s(&bl[0], NULL, bl.shape[0] - 1, &x[0], &out[0], x.shape[0])
    return np.asarray(out)


def legendre_roots(n):
    x = np.empty(n, np.double)
    w = np.empty(n, np.double)
    cdef double[::1] x_buf = x, w_buf = w
    if not (x_buf.shape[0] == w_buf.shape[0] == n):
        raise AssertionError()
    if n > 0:
        sharp_legendre_roots(n, &x_buf[0], &w_buf[0])
    return x, w


JOBTYPE_TO_CONST = {
    'Y': SHARP_Y,
    'Yt': SHARP_Yt,
    'WY': SHARP_WY,
    'YtW': SHARP_YtW
}

def sht(jobtype, geom_info ginfo, alm_info ainfo, double[:, :, ::1] input,
        int spin=0, comm=None, add=False):
    cdef void *comm_ptr
    cdef int flags = SHARP_DP | (SHARP_ADD if add else 0)
    cdef int r
    cdef sharp_jobtype jobtype_i
    cdef double[:, :, ::1] output_buf
    cdef int ntrans = input.shape[0]
    cdef int ntotcomp = ntrans * input.shape[1]
    cdef int i, j

    if spin == 0 and input.shape[1] != 1:
        raise ValueError('For spin == 0, we need input.shape[1] == 1')
    elif spin != 0 and input.shape[1] != 2:
        raise ValueError('For spin != 0, we need input.shape[1] == 2')


    cdef size_t[::1] ptrbuf = np.empty(2 * ntotcomp, dtype=np.uintp)
    cdef double **alm_ptrs = <double**>&ptrbuf[0]
    cdef double **map_ptrs = <double**>&ptrbuf[ntotcomp]

    try:
        jobtype_i = JOBTYPE_TO_CONST[jobtype]
    except KeyError:
        raise ValueError('jobtype must be one of: %s' % ', '.join(sorted(JOBTYPE_TO_CONST.keys())))

    if jobtype_i == SHARP_Y or jobtype_i == SHARP_WY:
        output = np.empty((input.shape[0], input.shape[1], ginfo.local_size()), dtype=np.float64)
        output_buf = output
        for i in range(input.shape[0]):
            for j in range(input.shape[1]):
                alm_ptrs[i * input.shape[1] + j] = &input[i, j, 0]
                map_ptrs[i * input.shape[1] + j] = &output_buf[i, j, 0]
    else:
        output = np.empty((input.shape[0], input.shape[1], ainfo.local_size()), dtype=np.float64)
        output_buf = output
        for i in range(input.shape[0]):
            for j in range(input.shape[1]):
                alm_ptrs[i * input.shape[1] + j] = &output_buf[i, j, 0]
                map_ptrs[i * input.shape[1] + j] = &input[i, j, 0]

    if comm is None:
        with nogil:
            sharp_execute (
                jobtype_i,
                geom_info=ginfo.ginfo, alm_info=ainfo.ainfo,
                spin=spin, alm=alm_ptrs, map=map_ptrs,
                ntrans=ntrans, flags=flags, time=NULL, opcnt=NULL)
    else:
        from mpi4py import MPI
        if not isinstance(comm, MPI.Comm):
            raise TypeError('comm must be an mpi4py communicator')
        from .libsharp_mpi import _addressof
        comm_ptr = <void*><size_t>_addressof(comm)
        with nogil:
            r = sharp_execute_mpi_maybe (
                comm_ptr, jobtype_i,
                geom_info=ginfo.ginfo, alm_info=ainfo.ainfo,
                spin=spin, alm=alm_ptrs, map=map_ptrs,
                ntrans=ntrans, flags=flags, time=NULL, opcnt=NULL)
        if r == SHARP_ERROR_NO_MPI:
            raise Exception('MPI requested, but not available')

    return output


def synthesis(*args, **kw):
    return sht('Y', *args, **kw)

def adjoint_synthesis(*args, **kw):
    return sht('Yt', *args, **kw)

def analysis(*args, **kw):
    return sht('YtW', *args, **kw)

def adjoint_analysis(*args, **kw):
    return sht('WY', *args, **kw)


#
# geom_info
#
class NotInitializedError(Exception):
    pass


cdef class geom_info:
    cdef sharp_geom_info *ginfo

    def __cinit__(self, *args, **kw):
        self.ginfo = NULL

    def local_size(self):
        if self.ginfo == NULL:
            raise NotInitializedError()
        return sharp_map_size(self.ginfo)

    def __dealloc__(self):
        if self.ginfo != NULL:
            sharp_destroy_geom_info(self.ginfo)
        self.ginfo = NULL


cdef class healpix_grid(geom_info):

    _weight_cache = {}  # { (nside, 'T'/'Q'/'U') -> numpy array of ring weights cached from file }

    def __init__(self, int nside, stride=1, int[::1] rings=None, double[::1] weights=None):
        if weights is not None and weights.shape[0] != 2 * nside:
            raise ValueError('weights must have length 2 * nside')
        sharp_make_subset_healpix_geom_info(nside, stride,
                                            nrings=4 * nside - 1 if rings is None else rings.shape[0],
                                            rings=NULL if rings is None else &rings[0],
                                            weight=NULL if weights is None else &weights[0],
                                            geom_info=&self.ginfo)

    @classmethod
    def load_ring_weights(cls, nside, fields):
        """
        Loads HEALPix ring weights from file. The environment variable
        HEALPIX should be set, and this routine will look in the `data`
        subdirectory.

        Parameters
        ----------

        nside: int
            HEALPix nside parameter

        fields: tuple of str
            Which weights to extract; pass ('T',) to only get scalar
            weights back, or ('T', 'Q', 'U') to get all the weights

        Returns
        -------

        List of NumPy arrays, according to fields parameter.

        """
        import os
        from astropy.io import fits
        data_path = os.path.join(os.environ['HEALPIX'], 'data')
        fits_field_names = {
            'T': 'TEMPERATURE WEIGHTS',
            'Q': 'Q-POLARISATION WEIGHTS',
            'U': 'U-POLARISATION WEIGHTS'}

        must_load = [field for field in fields if (nside, field) not in cls._weight_cache]

        if must_load:
            hdulist = fits.open(os.path.join(data_path, 'weight_ring_n%05d.fits' % nside))
            try:
                for field in must_load:
                    w = hdulist[1].data.field(fits_field_names[field]).ravel().astype(np.double)
                    w += 1
                    cls._weight_cache[nside, field] = w
            finally:
                hdulist.close()
        return [cls._weight_cache[(nside, field)].copy() for field in fields]

#
# alm_info
#


cdef class alm_info:
    cdef sharp_alm_info *ainfo

    def __cinit__(self, *args, **kw):
        self.ainfo = NULL

    def local_size(self):
        if self.ainfo == NULL:
            raise NotInitializedError()
        return sharp_alm_count(self.ainfo)

    def mval(self):
        if self.ainfo == NULL:
            raise NotInitializedError()
        return np.asarray(<int[:self.ainfo.nm]> self.ainfo.mval)

    def mvstart(self):
        if self.ainfo == NULL:
            raise NotInitializedError()
        return np.asarray(<long[:self.ainfo.nm]> self.ainfo.mvstart)

    def __dealloc__(self):
        if self.ainfo != NULL:
            sharp_destroy_alm_info(self.ainfo)
        self.ainfo = NULL

    @cython.boundscheck(False)
    def almxfl(self, np.ndarray[double, ndim=3, mode='c'] alm, np.ndarray[double, ndim=2, mode='c'] fl):
        """Multiply Alm by a Ell based array


        Parameters
        ----------
        alm : np.ndarray
            input alm, 3 dimensions = (different signal x polarizations x lm-ordering)
        fl : np.ndarray
            either 1 dimension, e.g. gaussian beam, or 2 dimensions e.g. a polarized beam

        Returns
        -------
        None, it modifies alms in-place

        """
        cdef int mvstart = 0
        cdef bint has_multiple_beams = alm.shape[2] > 1 and fl.shape[1] > 1
        cdef int f, i_m, m, num_ells, i_l, i_signal, i_pol, i_mv

        for i_m in range(self.ainfo.nm):
            m = self.ainfo.mval[i_m]
            f = 1 if (m==0) else 2
            num_ells = self.ainfo.lmax + 1 - m

            if not has_multiple_beams:
                for i_signal in range(alm.shape[0]):
                    for i_pol in range(alm.shape[1]):
                        for i_l in range(num_ells):
                            l = m + i_l
                            for i_mv in range(mvstart + f*i_l, mvstart + f*i_l +f):
                                alm[i_signal, i_pol, i_mv] *= fl[l, 0]
            else:
                for i_signal in range(alm.shape[0]):
                    for i_pol in range(alm.shape[1]):
                        for i_l in range(num_ells):
                            l = m + i_l
                            for i_mv in range(mvstart + f*i_l, mvstart + f*i_l +f):
                                alm[i_signal, i_pol, i_mv] *= fl[l, i_pol]
            mvstart += f * num_ells

cdef class triangular_order(alm_info):
    def __init__(self, int lmax, mmax=None, stride=1):
        mmax = mmax if mmax is not None else lmax
        sharp_make_triangular_alm_info(lmax, mmax, stride, &self.ainfo)


cdef class rectangular_order(alm_info):
    def __init__(self, int lmax, mmax=None, stride=1):
        mmax = mmax if mmax is not None else lmax
        sharp_make_rectangular_alm_info(lmax, mmax, stride, &self.ainfo)


cdef class packed_real_order(alm_info):
    def __init__(self, int lmax, stride=1, int[::1] ms=None):
        sharp_make_mmajor_real_packed_alm_info(lmax=lmax, stride=stride,
                                               nm=lmax + 1 if ms is None else ms.shape[0],
                                               ms=NULL if ms is None else &ms[0],
                                               alm_info=&self.ainfo)

#
# 
#

@cython.boundscheck(False)
def normalized_associated_legendre_table(int lmax, int m, theta):
    cdef double[::1] theta_ = np.ascontiguousarray(theta, dtype=np.double)
    out = np.zeros((theta_.shape[0], lmax - m + 1), np.double)
    cdef double[:, ::1] out_ = out
    if lmax < m:
        raise ValueError("lmax < m")
    with nogil:
        sharp_normalized_associated_legendre_table(m, 0, lmax, theta_.shape[0], &theta_[0], lmax - m + 1, 1, 1, &out_[0,0])
    return out
