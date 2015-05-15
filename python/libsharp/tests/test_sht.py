import numpy as np
import healpy
from scipy.special import legendre
from scipy.special import p_roots
from numpy.testing import assert_allclose
import libsharp

from mpi4py import MPI


def test_basic():
    lmax = 10
    nside = 8
    rank = MPI.COMM_WORLD.Get_rank()
    ms = np.arange(rank, lmax + 1, MPI.COMM_WORLD.Get_size(), dtype=np.int32)
    
    order = libsharp.packed_real_order(lmax, ms=ms)
    grid = libsharp.healpix_grid(nside)

    
    alm = np.zeros(order.local_size())
    if rank == 0:
        alm[0] = 1
    elif rank == 1:
        alm[0] = 1


    
    map = libsharp.synthesis(grid, order, alm, comm=MPI.COMM_WORLD)
    if rank == 0:
        healpy.mollzoom(map)
        from matplotlib.pyplot import show
        show()
