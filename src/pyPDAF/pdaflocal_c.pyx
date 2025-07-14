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

def set_indices(int  dim_l, int [::1] map):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_l : int 
        Dimension of local state vector
    map : ndarray[np.intc, ndim=1]
        Index array for mapping between local and global state vector
        Array shape: (dim_l)

    Returns
    -------
    """
    with nogil:
        c__pdaflocal_set_indices(&dim_l, &map[0])



def set_increment_weights(int  dim_l, double [::1] weights):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_l : int 
        Dimension of local state vector
    weights : ndarray[np.float64, ndim=1]
        Weights array
        Array shape: (dim_l)

    Returns
    -------
    """
    with nogil:
        c__pdaflocal_set_increment_weights(&dim_l, &weights[0])



def clear_increment_weights():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaflocal_clear_increment_weights()



