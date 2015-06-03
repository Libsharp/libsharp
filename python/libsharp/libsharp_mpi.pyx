cdef extern from "mpi.h":
    ctypedef void *MPI_Comm

cdef extern from "Python.h":
    object PyLong_FromVoidPtr(void*)

cdef extern:
    ctypedef class mpi4py.MPI.Comm [object PyMPICommObject]:
        cdef MPI_Comm ob_mpi
        cdef unsigned flags

# For compatibility with mpi4py <= 1.3.1
# Newer versions could use the MPI._addressof function
def _addressof(Comm comm):
    cdef void *ptr = NULL
    ptr = <void*>&comm.ob_mpi
    return PyLong_FromVoidPtr(ptr)
