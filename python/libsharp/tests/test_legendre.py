import numpy as np
from scipy.special import legendre
from scipy.special import p_roots
import libsharp

from numpy.testing import assert_allclose


def check_legendre_transform(lmax, ntheta):
    l = np.arange(lmax + 1)
    if lmax >= 1:
        sigma = -np.log(1e-3) / lmax / (lmax + 1)
        bl = np.exp(-sigma*l*(l+1))
        bl *= (2 * l + 1)
    else:
        bl = np.asarray([1], dtype=np.double)

    theta = np.linspace(0, np.pi, ntheta, endpoint=True)
    x = np.cos(theta)

    # Compute truth using scipy.special.legendre
    P = np.zeros((ntheta, lmax + 1))
    for l in range(lmax + 1):
        P[:, l] = legendre(l)(x)
    y0 = np.dot(P, bl)


    # double-precision
    y = libsharp.legendre_transform(x, bl)

    assert_allclose(y, y0, rtol=1e-12, atol=1e-12)

    # single-precision
    y32 = libsharp.legendre_transform(x.astype(np.float32), bl)
    assert_allclose(y, y0, rtol=1e-5, atol=1e-5)


def test_legendre_transform():
    nthetas_to_try = [0, 9, 17, 19] + list(np.random.randint(500, size=20))
    for ntheta in nthetas_to_try:
        for lmax in [0, 1, 2, 3, 20] + list(np.random.randint(50, size=4)):
            yield check_legendre_transform, lmax, ntheta

def check_legendre_roots(n):
    xs, ws = ([], []) if n == 0 else p_roots(n) # from SciPy
    xl, wl = libsharp.legendre_roots(n)
    assert_allclose(xs, xl, rtol=1e-14, atol=1e-14)
    assert_allclose(ws, wl, rtol=1e-14, atol=1e-14)

def test_legendre_roots():
    """
    Test the Legendre root-finding algorithm from libsharp by comparing it with
    the SciPy version.
    """
    yield check_legendre_roots, 0
    yield check_legendre_roots, 1
    yield check_legendre_roots, 32
    yield check_legendre_roots, 33
    yield check_legendre_roots, 128
