from __future__ import print_function
import numpy as np

from numpy.testing import assert_almost_equal
from nose.tools import eq_, ok_

from libsharp import normalized_associated_legendre_table
from scipy.special import sph_harm, p_roots

def test_compare_legendre_table_with_scipy():
    def test(theta, m, lmax):
        Plm = normalized_associated_legendre_table(lmax, m, theta)

        Plm_p = sph_harm(m, np.arange(m, lmax + 1), 0, theta)[None, :]
        if not np.allclose(Plm_p, Plm):
            print(Plm_p)
            print(Plm)
        return ok_, np.allclose(Plm_p, Plm)

    yield test(np.pi/2, 0, 10)
    yield test(np.pi/4, 0, 10)
    yield test(3 * np.pi/4, 0, 10)
    yield test(np.pi/4, 1, 4)
    yield test(np.pi/4, 2, 4)
    yield test(np.pi/4, 50, 50)
    yield test(np.pi/2, 49, 50)


def test_legendre_table_wrapper_logic():
    # tests the SSE 2 logic in the high-level wrapper by using an odd number of thetas
    theta = np.asarray([np.pi/2, np.pi/4, 3 * np.pi / 4])
    m = 3
    lmax = 10
    Plm = normalized_associated_legendre_table(lmax, m, theta)
    assert np.allclose(Plm[1, :], normalized_associated_legendre_table(lmax, m, np.pi/4)[0, :])
    assert np.allclose(Plm[2, :], normalized_associated_legendre_table(lmax, m, 3 * np.pi/4)[0, :])
