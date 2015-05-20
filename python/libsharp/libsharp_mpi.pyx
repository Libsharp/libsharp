from mpi4py.MPI cimport Comm
cdef extern from "Python.h":
    object PyLong_FromVoidPtr(void*)

# For compatibility with mpi4py <= 1.3.1
# Newer versions could use the MPI._addressof function
def _addressof(comm):
    cdef void *ptr = NULL
    ptr = <void*>&(<Comm>comm).ob_mpi
    return PyLong_FromVoidPtr(ptr)
