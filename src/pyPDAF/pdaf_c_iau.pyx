import sys
import numpy as np
cimport numpy as cnp
from . cimport pdaf_c_cb_interface as pdaf_cb
from .cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from .cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from .cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

try:
    import mpi4py
    mpi4py.rc.initialize = False
except ImportError:
    pass

# Global error handler
def global_except_hook(exctype, value, traceback):
    from traceback import print_exception
    try:
        import mpi4py.MPI

        if mpi4py.MPI.Is_initialized():
            try:
                sys.stderr.write('Uncaught exception was ''detected on rank {}.\n'.format(
                    mpi4py.MPI.COMM_WORLD.Get_rank()))
                print_exception(exctype, value, traceback)
                sys.stderr.write("\n")
                sys.stderr.flush()
            finally:
                try:
                    mpi4py.MPI.COMM_WORLD.Abort(1)
                except Exception as e:
                    sys.stderr.write('MPI Abort failed, this process will hang.\n')
                    sys.stderr.flush()
                    raise e
        else:
            sys.__excepthook__(exctype, value, traceback)
    except ImportError:
        sys.__excepthook__(exctype, value, traceback)

sys.excepthook = global_except_hook

def iau_init(int  type_iau_in, int  nsteps_iau_in):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    type_iau_in : int 
        Type of IAU, (0) no IAU
    nsteps_iau_in : int 
        number of time steps in IAU

    Returns
    -------
    flag : int 
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_init(&type_iau_in, &nsteps_iau_in, &flag)

    return flag


def iau_reset(int  type_iau_in, int  nsteps_iau_in):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    type_iau_in : int 
        Type of IAU, (0) no IAU
    nsteps_iau_in : int 
        number of time steps in IAU

    Returns
    -------
    flag : int 
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_reset(&type_iau_in, &nsteps_iau_in, &flag)

    return flag


def iau_set_weights(int  iweights, double [::1] weights):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    iweights : int 
        Length of weights input vector
    weights : ndarray[np.float64, ndim=1]
        Input weight vector
        Array shape: (iweights)

    Returns
    -------
    """
    with nogil:
        c__pdaf_iau_set_weights(&iweights, &weights[0])



def iau_set_pointer():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    cdef CFI_cdesc_rank2 iau_ptr_cfi
    cdef CFI_cdesc_t *iau_ptr_ptr = <CFI_cdesc_t *> &iau_ptr_cfi
    cdef int  flag
    with nogil:
        c__pdaf_iau_set_pointer(iau_ptr_ptr, &flag)

    cdef CFI_index_t iau_ptr_subscripts[2]
    iau_ptr_subscripts[0] = 0
    iau_ptr_subscripts[1] = 0
    cdef double * iau_ptr_ptr_np
    iau_ptr_ptr_np = <double *>CFI_address(iau_ptr_ptr, iau_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] iau_ptr_np = np.asarray(<double [:iau_ptr_ptr.dim[0].extent:1,:iau_ptr_ptr.dim[1].extent]> iau_ptr_ptr_np, order="F")
    return iau_ptr_np, flag


def iau_init_inc(int  dim_p, int  dim_ens_l, double [::1,:] ens_inc):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int 
        PE-local dimension of model state
    dim_ens_l : int 
        Task-local size of ensemble
    ens_inc : ndarray[np.float64, ndim=2]
        PE-local increment ensemble
        Array shape: (dim_p, dim_ens_l)

    Returns
    -------
    flag : int 
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_init_inc(&dim_p, &dim_ens_l, &ens_inc[0,0], &flag)

    return flag


def iau_add_inc(py__collect_state_pdaf, py__distribute_state_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

        Callback Parameters
        -------------------
        dim_p : int 
                pe-local state dimension
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

        Callback Parameters
        -------------------
        dim_p : int 
                PE-local state dimension
        state_p : ndarray[np.float64, ndim=1]
                PE-local state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local state vector
                Array shape: (dim_p)


    Returns
    -------
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    with nogil:
        c__pdaf_iau_add_inc(pdaf_cb.c__collect_state_pdaf, 
                            pdaf_cb.c__distribute_state_pdaf)



