import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

def set_comm_pdaf(int  in_comm_pdaf):
    """set_comm_pdaf(in_comm_pdaf:int) -> None

    Set the MPI communicator used by PDAF.

    By default, PDAF assumes it can use all available
    processes, i.e., `MPI_COMM_WORLD`.
    By using this function, we limit the number of processes
    that can be used by PDAF to given MPI communicator.

    Parameters
    ----------
    in_comm_pdaf : int
        MPI communicator for PDAF
    """
    c__pdaf_set_comm_pdaf(&in_comm_pdaf)



def set_debug_flag(int  debugval):
    """set_debug_flag(debugval:int) -> None

    Activate the debug output of the PDAF.

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

    """
    c__pdaf_set_debug_flag(&debugval)



def set_ens_pointer():
    """set_ens_pointer() -> Tuple[np.ndarray, int]

    Return the ensemble in a numpy array.

    Here the internal array data has the same memoery address
    as PDAF ensemble array allowing for manual ensemble modification.

    Returns
    -------
    ens_ptr_np : np.ndarray
        Numpy array view of the ensemble
    status : int
        Status flag
    """
    cdef CFI_cdesc_rank2 ens_ptr_cfi
    cdef int  status
    c__pdaf_set_ens_pointer(<CFI_cdesc_t *> &ens_ptr_cfi, &status)

    cdef CFI_index_t ens_ptr_subscripts[2]
    ens_ptr_subscripts[0] = ens_ptr_cfi.dim[0].lower_bound
    ens_ptr_subscripts[1] = ens_ptr_cfi.dim[1].lower_bound
    cdef double * ens_ptr_ptr_np
    ens_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &ens_ptr_cfi, ens_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_ptr_np = np.asarray(<double [:ens_ptr_cfi.dim[0].extent:1,:ens_ptr_cfi.dim[1].extent]> ens_ptr_ptr_np, order="F")
    return ens_ptr_np, status


def set_iparam(int  idval, int  value, int  flag):
    """set_iparam(idval: int, value: int, flag: int) -> int

    Set integer parameters for PDAF.

    The integer parameters specific to a DA method can be set in the array
    `filter_param_i` that is an argument of :func:`pyPDAF.PDAF.init`
    (see the page on `initializing PDAF <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_).

    This function provides an alternative way.
    Instead of providing all parameters in the call to :func:`pyPDAF.PDAF.init`,
    one can provide only the required minimum options for this call.
    Afterwards, one can then call this function for each integer parameter
    that one intends to specify differently from the default value.
    An advantage of using this function is that one only needs to call it
    for parameters that one intends to change, while in the call to
    :func:`pyPDAF.PDAF.init` all parameters up to the index one intends to change
    have to be specified, even if one does not want to change a parameter value.

    The routine is usually called by all processes after the call to
    :func:`pyPDAF.PDAF.init` in init_pdaf. One can also call the routine at a
    later time during an assimilation process to change parameters.
    The parameter will be set for the DA method that was specified
    in the call to :func:`pyPDAF.PDAF.init`.

    Parameters
    ----------
    idval : int
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
    c__pdaf_set_iparam(&idval, &value, &flag)

    return flag


def set_memberid(int  memberid):
    """set_memberid(int  memberid) -> int

    Set the ensemble member index to given value.

    Parameters
    ----------
    memberid : int
        Index in the local ensemble

    Returns
    -------
    memberid : int
        Index in the local ensemble
    """
    c__pdaf_set_memberid(&memberid)

    return memberid


def set_offline_mode(int  screen):
    """set_offline_mode(screen: int) -> None

    Activate offline mode of PDAF.

    Parameters
    ----------
    screen : int
        Verbosity flag

    """
    c__pdaf_set_offline_mode(&screen)



def set_rparam(int  idval, double  value, int  flag):
    """set_rparam(idval: int, value: float, flag: int) -> int

    Set floating-point parameters for PDAF.

    The floating-point parameters specific to a DA method can be set in the array
    `filter_param_r` that is an argument of :func:`pyPDAF.PDAF.init`
    (see the page on `initializing PDAF <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_).

    This function provides an alternative way.
    Instead of providing all parameters in the call to :func:`pyPDAF.PDAF.init`,
    one can provide only the required minimum options for this call.
    Afterwards, one can then call this function for each floating-point parameter
    that one intends to specify differently from the default value.
    An advantage of using this function is that one only needs to call it
    for parameters that one intends to change, while in the call to
    :func:`pyPDAF.PDAF.init` all parameters up to the index one intends to change
    have to be specified, even if one does not want to change a parameter value.

    The routine is usually called by all processes after the call to
    :func:`pyPDAF.PDAF.init` in init_pdaf. One can also call the routine at a
    later time during an assimilation process to change parameters.
    The parameter will be set for the DA method that was specified
    in the call to :func:`pyPDAF.PDAF.init`.

    Parameters
    ----------
    idval : int
        Index of parameter
    value : float
        Parameter value
    flag : int
        Status flag: 0 for no error

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    c__pdaf_set_rparam(&idval, &value, &flag)

    return flag

def genobs_set_rparam(int idval, double value):
    """genobs_set_rparam(idval: int, value: float) -> int

    Set a real-valued GENOBS parameter.

    This routine is specific to PDAF's generated-observation mode. It updates
    a real-valued GENOBS control parameter after the PDAF method has been
    selected. The meaning of ``idval`` is defined by PDAF's GENOBS parameter
    table.

    Parameters
    ----------
    idval : int
        Identifier of the GENOBS real parameter to update.
    value : float
        New parameter value.

    Returns
    -------
    flag : int
        PDAF status flag. A value of ``0`` indicates success.
    """
    cdef int flag
    c__pdaf_genobs_set_rparam(&idval, &value, &flag)
    return flag


def set_seedset(int  seedset_in):
    """set_seedset(seedset_in: int) -> None

    Choose a seedset for the random number generator used in PDAF.

    Parameters
    ----------
    seedset_in : int
        Seedset index (1-20)

    """
    c__pdaf_set_seedset(&seedset_in)


def set_seed(int [::1] seedvec):
    """set_seed(seedvec: np.ndarray) -> None

    Set PDAF's four-integer random seed vector.

    This is an alias of :func:`set_seedvec`. Use it before routines that draw
    PDAF random numbers when an assimilation run has to be reproducible. The
    fourth entry of ``seedvec`` must be odd.

    Parameters
    ----------
    seedvec : ndarray[np.intc, ndim=1]
        Four-integer seed vector following the LAPACK convention used by PDAF.

    Returns
    -------
    None
    """
    c__pdaf_set_seed(&seedvec[0])

def set_seedvec(int [::1] seedvec):
    """set_seedvec(seedvec: np.ndarray) -> None

    Set PDAF's four-integer random seed vector.

    The seed vector controls random numbers generated by PDAF helper routines,
    including random ensemble or perturbation generation. Store the result of
    :func:`pyPDAF.PDAF.get_seedvec` and pass it back here to continue a random
    sequence from a known point. The fourth entry must be odd.

    Parameters
    ----------
    seedvec : ndarray[np.intc, ndim=1]
        Four-integer seed vector following the LAPACK convention used by PDAF.

    Returns
    -------
    None
    """
    c__pdaf_set_seedvec(&seedvec[0])


def set_smoother_ens(int  maxlag):
    """set_smoother_ens(maxlag: int) -> Tuple[np.ndarray, int]

    Get a pointer to smoother ensemble.

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
    cdef int  status
    c__pdaf_set_smootherens(<CFI_cdesc_t *> &sens_point_cfi, &maxlag, &status)

    cdef CFI_index_t sens_point_subscripts[3]
    sens_point_subscripts[0] = sens_point_cfi.dim[0].lower_bound
    sens_point_subscripts[1] = sens_point_cfi.dim[1].lower_bound
    sens_point_subscripts[2] = sens_point_cfi.dim[2].lower_bound
    cdef double * sens_point_ptr_np
    sens_point_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &sens_point_cfi, sens_point_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_point_np = np.asarray(<double [:sens_point_cfi.dim[0].extent:1,:sens_point_cfi.dim[1].extent,:sens_point_cfi.dim[2].extent]> sens_point_ptr_np, order="F")
    return sens_point_np, status




