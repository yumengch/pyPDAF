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

def get_assim_flag():
    """Return the flag that
    indicates if the DA is performed in the last time step.
    It only works for online DA systems. 
    """
    cdef int  did_assim
    with nogil:
        c__pdaf_get_assim_flag(&did_assim)

    return did_assim


def get_localfilter():
    """Return whether a local filter is used. 
    """
    cdef int  localfilter_out
    with nogil:
        c__pdaf_get_localfilter(&localfilter_out)

    return localfilter_out


def get_local_type():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    cdef int  localtype
    with nogil:
        c__pdaf_get_local_type(&localtype)

    return localtype


def get_memberid(int  memberid):
    """Return the ensemble member id on the current process.

    For example, it can be called during the ensemble
    integration if ensemble-specific forcing is read.
    It can also be used in the user-supplied functions
    such as :func:`py__collect_state_pdaf` and
    :func:`py__distribute_state_pdaf`.

    Parameters
    ----------
    memberid : int 
        Index in the local ensemble

    Returns
    -------
    memberid : int 
        Index in the local ensemble
    """
    with nogil:
        c__pdaf_get_memberid(&memberid)

    return memberid


def get_obsmemberid(int  memberid):
    """Return the ensemble member id
    when observation operator is being applied.

    This function is used specifically for
    user-supplied function :func:`py__obs_op_pdaf`.

    Parameters
    ----------
    memberid : int 
        Index in the local ensemble

    Returns
    -------
    memberid : int 
        Index in the local ensemble
    """
    with nogil:
        c__pdaf_get_obsmemberid(&memberid)

    return memberid


def get_smootherens():
    """Return the smoothed ensemble in earlier time steps.

    It is only used when the smoother options is used .
    """
    cdef CFI_cdesc_rank3 sens_point_cfi
    cdef CFI_cdesc_t *sens_point_ptr = <CFI_cdesc_t *> &sens_point_cfi
    cdef int  maxlag
    cdef int  status
    with nogil:
        c__pdaf_get_smootherens(sens_point_ptr, &maxlag, &status)

    cdef CFI_index_t sens_point_subscripts[3]
    sens_point_subscripts[0] = 0
    sens_point_subscripts[1] = 0
    sens_point_subscripts[2] = 0
    cdef double * sens_point_ptr_np
    sens_point_ptr_np = <double *>CFI_address(sens_point_ptr, sens_point_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_point_np = np.asarray(<double [:sens_point_ptr.dim[0].extent:1,:sens_point_ptr.dim[1].extent,:sens_point_ptr.dim[2].extent]> sens_point_ptr_np, order="F")
    return sens_point_np, maxlag, status


