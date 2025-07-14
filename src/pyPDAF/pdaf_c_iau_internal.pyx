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

def _iau_init_weights(int  type_iau, int  nsteps_iau):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    type_iau : int 
        Type of IAU, (0) no IAU
    nsteps_iau : int 
        number of time steps in IAU

    Returns
    -------
    """
    with nogil:
        c__pdaf_iau_init_weights(&type_iau, &nsteps_iau)



def _iau_update_inc(double [::1,:] ens_ana):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    ens_ana : ndarray[np.float64, ndim=2]
        PE-local analysis ensemble
        Array shape: (:, :)

    Returns
    -------
    ens_ana : ndarray[np.float64, ndim=2]
        PE-local analysis ensemble
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 ens_ana_cfi
    cdef CFI_cdesc_t *ens_ana_ptr = <CFI_cdesc_t *> &ens_ana_cfi
    cdef size_t ens_ana_nbytes = ens_ana.nbytes
    cdef CFI_index_t ens_ana_extent[2]
    ens_ana_extent[0] = ens_ana.shape[0]
    ens_ana_extent[1] = ens_ana.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_ana_np = np.asarray(ens_ana, dtype=np.float64, order="F")
    with nogil:
        CFI_establish(ens_ana_ptr, &ens_ana[0,0], CFI_attribute_other,
                      CFI_type_double , ens_ana_nbytes, 2, ens_ana_extent)

        c__pdaf_iau_update_inc(ens_ana_ptr)

    return ens_ana_np


def _iau_add_inc_ens(int  step, int  dim_p, int  dim_ens_task, 
    double [::1,:] ens, py__collect_state_pdaf, py__distribute_state_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int 
        Time step
    dim_p : int 
        PE-local dimension of model state
    dim_ens_task : int 
        Ensemble size of model task
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
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
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 ens_cfi
    cdef CFI_cdesc_t *ens_ptr = <CFI_cdesc_t *> &ens_cfi
    cdef size_t ens_nbytes = ens.nbytes
    cdef CFI_index_t ens_extent[2]
    ens_extent[0] = ens.shape[0]
    ens_extent[1] = ens.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.asarray(ens, dtype=np.float64, order="F")
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    with nogil:
        CFI_establish(ens_ptr, &ens[0,0], CFI_attribute_other,
                      CFI_type_double , ens_nbytes, 2, ens_extent)

        c__pdaf_iau_add_inc_ens(&step, &dim_p, &dim_ens_task, ens_ptr, 
                                pdaf_cb.c__collect_state_pdaf, 
                                pdaf_cb.c__distribute_state_pdaf)

    return ens_np


def _iau_update_ens(double [::1,:] ens):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)

    Returns
    -------
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 ens_cfi
    cdef CFI_cdesc_t *ens_ptr = <CFI_cdesc_t *> &ens_cfi
    cdef size_t ens_nbytes = ens.nbytes
    cdef CFI_index_t ens_extent[2]
    ens_extent[0] = ens.shape[0]
    ens_extent[1] = ens.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.asarray(ens, dtype=np.float64, order="F")
    with nogil:
        CFI_establish(ens_ptr, &ens[0,0], CFI_attribute_other,
                      CFI_type_double , ens_nbytes, 2, ens_extent)

        c__pdaf_iau_update_ens(ens_ptr)

    return ens_np


def _iau_dealloc():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_iau_dealloc()



