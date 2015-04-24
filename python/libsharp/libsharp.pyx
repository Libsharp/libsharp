import numpy as np

cdef extern from "sharp.h":
    ctypedef long ptrdiff_t

    void sharp_legendre_transform_s(float *bl, float *recfac, ptrdiff_t lmax, float *x,
                                    float *out, ptrdiff_t nx)
    void sharp_legendre_transform(double *bl, double *recfac, ptrdiff_t lmax, double *x,
                                  double *out, ptrdiff_t nx)
    void sharp_legendre_transform_recfac(double *r, ptrdiff_t lmax)
    void sharp_legendre_transform_recfac_s(float *r, ptrdiff_t lmax)


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

