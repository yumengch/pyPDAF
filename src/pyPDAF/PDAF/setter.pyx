import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

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

def set_iparam_filters(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_set_iparam_filters(&id, &value, &flag)

    return flag


def set_rparam_filters(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_set_rparam_filters(&id, &value, &flag)

    return flag


def set_comm_pdaf(int  in_comm_pdaf):
    """Set the MPI communicator used by PDAF.

    By default, PDAF assumes it can use all available
    processes, i.e., `MPI_COMM_WORLD`.
    By using this function, we limit the number of processes
    that can be used by PDAF to given MPI communicator.

    Parameters
    ----------
    in_comm_pdaf : int
        MPI communicator for PDAF

    Returns
    -------
    """
    with nogil:
        c__pdaf_set_comm_pdaf(&in_comm_pdaf)



def set_debug_flag(int  debugval):
    """Activate the debug output of the PDAF.

    Starting from the use of this function,
    the debug infomation is sent to screen output.
    The screen output end when the debug flag is
    set to 0 by this function.

    For the sake of simplicity,
    we recommend using debugging output for
    a single local domain, e.g.
    `if domain_p == 1: pyPDAF.PDAF.set_debug_flag(1)`

    Parameters
    ----------
    debugval : int
        Value for debugging flag

    Returns
    -------
    """
    with nogil:
        c__pdaf_set_debug_flag(&debugval)



def set_ens_pointer():
    """Return the ensemble in a numpy array.

    Here the internal array data has the same memoery address
    as PDAF ensemble array allowing for manual ensemble modification.
    """
    cdef CFI_cdesc_rank2 ens_ptr_cfi
    cdef CFI_cdesc_t *ens_ptr_ptr = <CFI_cdesc_t *> &ens_ptr_cfi
    cdef int  status
    with nogil:
        c__pdaf_set_ens_pointer(ens_ptr_ptr, &status)

    cdef CFI_index_t ens_ptr_subscripts[2]
    ens_ptr_subscripts[0] = 0
    ens_ptr_subscripts[1] = 0
    cdef double * ens_ptr_ptr_np
    ens_ptr_ptr_np = <double *>CFI_address(ens_ptr_ptr, ens_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_ptr_np = np.asarray(<double [:ens_ptr_ptr.dim[0].extent:1,:ens_ptr_ptr.dim[1].extent]> ens_ptr_ptr_np, order="F")
    return ens_ptr_np, status


def set_iparam(int  id, int  value, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value
    flag : int
        Status flag: 0 for no error

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    with nogil:
        c__pdaf_set_iparam(&id, &value, &flag)

    return flag


def set_memberid(int  memberid):
    """Set the ensemble member index to given value.

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
        c__pdaf_set_memberid(&memberid)

    return memberid


def set_offline_mode(int  screen):
    """Activate offline mode of PDAF.

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    """
    with nogil:
        c__pdaf_set_offline_mode(&screen)



def set_rparam(int  id, double  value, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value
    flag : int
        Status flag: 0 for no error

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    with nogil:
        c__pdaf_set_rparam(&id, &value, &flag)

    return flag


def set_seedset(int  seedset_in):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    seedset_in : int
        Seedset index (1-20)

    Returns
    -------
    """
    with nogil:
        c__pdaf_set_seedset(&seedset_in)



def set_smootherens(int  maxlag):
    """Get a pointer to smoother ensemble.

    When smoother is used, the smoothed ensemble states
    at earlier times are stored in an internal array of PDAF.
    To be able to smooth post times,
    the smoother algorithm must have access to the past ensembles.

    In this function, the user can obtain a numpy array of
    smoother ensemble. This array has the same memory address
    as the internal PDAF smoother ensemble array.
    This allows for manual modification of the smoother ensemble.

    In the offline mode the user has to manually
    fill the smoother ensemble array
    from ensembles read in from files.
    This function is typically called in
    :func:`py__init_ens_pdaf` in the call to
    :func:`pyPDAF.PDAF.PDAF_init`.

    In the online mode, the smoother array is filled
    automatically during the cycles of forecast phases and analysis steps.

    Parameters
    ----------
    maxlag : int
        Number of past timesteps in sens

    Returns
    -------
    sens_point : ndarray[np.float64, ndim=3]
        Pointer to smoother array
        Array shape: (:,:,:)
    status : int
        Status flag,
    """
    cdef CFI_cdesc_rank3 sens_point_cfi
    cdef CFI_cdesc_t *sens_point_ptr = <CFI_cdesc_t *> &sens_point_cfi
    cdef int  status
    with nogil:
        c__pdaf_set_smootherens(sens_point_ptr, &maxlag, &status)

    cdef CFI_index_t sens_point_subscripts[3]
    sens_point_subscripts[0] = 0
    sens_point_subscripts[1] = 0
    sens_point_subscripts[2] = 0
    cdef double * sens_point_ptr_np
    sens_point_ptr_np = <double *>CFI_address(sens_point_ptr, sens_point_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_point_np = np.asarray(<double [:sens_point_ptr.dim[0].extent:1,:sens_point_ptr.dim[1].extent,:sens_point_ptr.dim[2].extent]> sens_point_ptr_np, order="F")
    return sens_point_np, status


def set_forget(int  step, int  localfilter, int  dim_obs_p, int  dim_ens,
    double [::1,:] mens_p, double [::1] mstate_p, double [::1] obs_p,
    py__init_obsvar_pdaf, double  forget_in, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    localfilter : int
        Whether filter is domain-local
    dim_obs_p : int
        Dimension of observation vector
    dim_ens : int
        Ensemble size
    mens_p : ndarray[np.float64, ndim=2]
        Observed PE-local ensemble
        Array shape: (dim_obs_p, dim_ens)
    mstate_p : ndarray[np.float64, ndim=1]
        Observed PE-local mean state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        Observation vector
        Array shape: (dim_obs_p)
    py__init_obsvar_pdaf : Callable
        Initialize mean obs. error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    forget_in : double
        Prescribed forgetting factor
    screen : int
        Verbosity flag

    Returns
    -------
    forget_out : double
        Adaptively estimated forgetting factor
    """
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef double  forget_out
    with nogil:
        c__pdaf_set_forget(&step, &localfilter, &dim_obs_p, &dim_ens,
                           &mens_p[0,0], &mstate_p[0], &obs_p[0],
                           pdaf_cb.c__init_obsvar_pdaf, &forget_in,
                           &forget_out, &screen)

    return forget_out


