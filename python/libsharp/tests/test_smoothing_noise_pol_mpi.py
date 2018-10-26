# This test needs to be run with:

# mpirun -np X python test_smoothing_noise_pol_mpi.py

from mpi4py import MPI

import numpy as np

import healpy as hp

import libsharp

mpi = True
rank = MPI.COMM_WORLD.Get_rank()

nside = 256
npix = hp.nside2npix(nside)

np.random.seed(100)
input_map = np.random.normal(size=(3, npix))
fwhm_deg = 10
lmax = 512

nrings = 4 * nside - 1  # four missing pixels

if rank == 0:
    print("total rings", nrings)

n_mpi_processes = MPI.COMM_WORLD.Get_size()
rings_per_process = nrings // n_mpi_processes + 1
# ring indices are 1-based

ring_indices_emisphere = np.arange(2*nside, dtype=np.int32) + 1
local_ring_indices = ring_indices_emisphere[rank::n_mpi_processes]

# to improve performance, simmetric rings north/south need to be in the same rank
# therefore we use symmetry to create the full ring indexing

if local_ring_indices[-1] == 2 * nside:
    # has equator ring
    local_ring_indices = np.concatenate(
      [local_ring_indices[:-1],
       nrings - local_ring_indices[::-1] + 1]
    )
else:
    # does not have equator ring
    local_ring_indices = np.concatenate(
      [local_ring_indices,
       nrings - local_ring_indices[::-1] + 1]
    )

print("rank", rank, "n_rings", len(local_ring_indices))

if not mpi:
    local_ring_indices = None
grid = libsharp.healpix_grid(nside, rings=local_ring_indices)

# returns start index of the ring and number of pixels
startpix, ringpix, _, _, _ = hp.ringinfo(nside, local_ring_indices.astype(np.int64))

local_npix = grid.local_size()

def expand_pix(startpix, ringpix, local_npix):
    """Turn first pixel index and number of pixel in full array of pixels

    to be optimized with cython or numba
    """
    local_pix = np.empty(local_npix, dtype=np.int64)
    i = 0
    for start, num in zip(startpix, ringpix):
        local_pix[i:i+num] = np.arange(start, start+num)
        i += num
    return local_pix

local_pix = expand_pix(startpix, ringpix, local_npix)

local_map = input_map[:, local_pix]

local_hitmap = np.zeros(npix)
local_hitmap[local_pix] = 1
hp.write_map("hitmap_{}.fits".format(rank), local_hitmap, overwrite=True)

print("rank", rank, "npix", npix, "local_npix", local_npix, "local_map len", len(local_map), "unique pix", len(np.unique(local_pix)))

local_m_indices = np.arange(rank, lmax + 1, MPI.COMM_WORLD.Get_size(), dtype=np.int32)
if not mpi:
    local_m_indices = None

order = libsharp.packed_real_order(lmax, ms=local_m_indices) 
local_nl = order.local_size()
print("rank", rank, "local_nl", local_nl, "mval", order.mval())

mpi_comm = MPI.COMM_WORLD if mpi else None

# map2alm
# maps in libsharp are 3D, 2nd dimension is IQU, 3rd is pixel

alm_sharp_I = libsharp.analysis(grid, order,
                                np.ascontiguousarray(local_map[0].reshape((1, 1, -1))),
                                spin=0, comm=mpi_comm)
alm_sharp_P = libsharp.analysis(grid, order,
                                np.ascontiguousarray(local_map[1:].reshape((1, 2, -1))),
                                spin=2, comm=mpi_comm)

beam = hp.gauss_beam(fwhm=np.radians(fwhm_deg), lmax=lmax, pol=True)

print("Smooth")
# smooth in place (zonca implemented this function)
order.almxfl(alm_sharp_I, np.ascontiguousarray(beam[:, 0:1]))
order.almxfl(alm_sharp_P, np.ascontiguousarray(beam[:, (1, 2)]))

# alm2map

new_local_map_I = libsharp.synthesis(grid, order, alm_sharp_I, spin=0, comm=mpi_comm)
new_local_map_P = libsharp.synthesis(grid, order, alm_sharp_P, spin=2, comm=mpi_comm)

# Transfer map to first process for writing

local_full_map = np.zeros(input_map.shape, dtype=np.float64)
local_full_map[0, local_pix] = new_local_map_I
local_full_map[1:, local_pix] = new_local_map_P

output_map = np.zeros(input_map.shape, dtype=np.float64) if rank == 0 else None
mpi_comm.Reduce(local_full_map, output_map, root=0, op=MPI.SUM)

if rank == 0:
    # hp.write_map("sharp_smoothed_map.fits", output_map, overwrite=True)
    # hp_smoothed = hp.alm2map(hp.map2alm(input_map, lmax=lmax), nside=nside) # transform only
    hp_smoothed = hp.smoothing(input_map, fwhm=np.radians(fwhm_deg), lmax=lmax)
    std_diff = (hp_smoothed-output_map).std()
    print("Std of difference between libsharp and healpy", std_diff)
    # hp.write_map(
    #     "healpy_smoothed_map.fits",
    #     hp_smoothed,
    #     overwrite=True
    # )
    assert std_diff < 1e-5
