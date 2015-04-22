import numpy as np
from scipy.special import legendre
import libsharp


def test_legendre_transform():
    lmax = 20
    ntheta = 19

    l = np.arange(lmax + 1)
    bl = np.exp(-l*(l+1))
    bl *= (2 * l + 1)

    theta = np.linspace(0, np.pi, ntheta, endpoint=True)
    x = np.cos(theta)

    # Compute truth using scipy.special.legendre
    P = np.zeros((ntheta, lmax + 1))
    for l in range(lmax + 1):
        P[:, l] = legendre(l)(x)
    y0 = np.dot(P, bl)

    # double-precision
    y = libsharp.legendre_transform(x, bl)
    assert np.max(np.abs(y - y) / np.abs(y)) < 1e-12

    # single-precision
    y32 = libsharp.legendre_transform(x.astype(np.float32), bl)
    assert np.max(np.abs(y32 - y) / np.abs(y32)) < 1e-4
