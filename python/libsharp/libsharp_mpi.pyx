"""
We keep MPI support in a seperate module because we need to import mpi4py,
and that will fire up the linked MPI library, which is something we only
want to do if a comm was passed.
"""

from mpi4py.MPI cimport Comm
from mpi4py.mpi_c cimport MPI_Comm

cimport libsharp

cdef extern from "sharp_mpi.h":
    void sharp_execute_mpi (MPI_Comm comm, sharp_jobtype type, int spin,
       void *alm, void *map, sharp_geom_info *geom_info,
       sharp_alm_info *alm_info, int ntrans, int flags, double *time,
       unsigned long long *opcnt) nogil

def _sht_mpi(MPI_comm comm, int jobtype, int spin, int flags, int ntrans, grid_info ginfo, double[::1] map,
             alm_info ainfo, double[::1] alm):
    with nogil:
        with nogil:
            sharp_execute_mpi(comm=comm.ob_mpi, jobtype=jobtype, spin=spin,
                              alm=&alm[0], map=&map[0], geom_info=ginfo.ginfo,
                              alm_info=ainfo.ainfo, ntrans=ntrans, flags=flags,
                              time=NULL, opcnt=NULL)
    
def sht_mpi(object comm, int jobtype, grid g, double[::1] map, alm_order, double[::1] alm, spin=0, comm=None, add=False):
    
def mpi4py_Comm_to_handle(Comm mpi4py_comm, capsule_to_c_comm):
    cdef MPI_Comm *out = <MPI_Comm*>PyCapsule_GetPointer(capsule_to_c_comm,
                                                         "_mpi4py_bridge_MPI_Comm")
    out[0] = mpi4py_comm.ob_mpi
