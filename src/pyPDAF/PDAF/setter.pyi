# pylint: disable=unused-argument
"""stub file for setter.pyx"""
import typing
import numpy as np


def set_comm_pdaf(in_comm_pdaf: int) -> None:
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

def set_debug_flag(debugval: int) -> None:
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

def set_ens_pointer() -> typing.Tuple[np.ndarray, int]:
    """Return the ensemble in a numpy array.

    Here the internal array data has the same memoery address
    as PDAF ensemble array allowing for manual ensemble modification.

    Returns
    -------
    ens_ptr_np : np.ndarray
        Numpy array view of the ensemble
    status : int
        Status flag
    """

def set_iparam(idval: int, value: int, flag: int) -> int:
    """Set integer parameters for PDAF.

    The integer parameters specific to a DA method can be set in the array
    `param_int` that is an argument of :func:`pyPDAF.PDAF.init`
    (see the page on
    `initializing PDAF <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_).

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

def set_memberid(memberid: int) -> int:
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

def set_offline_mode(screen: int) -> None:
    """Activate offline mode of PDAF.

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    """

def set_rparam(idval: int, value: float, flag: int) -> int:
    """Set floating-point parameters for PDAF.

    The floating-point parameters specific to a DA method can be set in the array
    `filter_param_r` that is an argument of :func:`pyPDAF.PDAF.init`
    (see the page on
    `initializing PDAF <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_).

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

def set_seedset(seedset_in: int) -> None:
    """Choose a seedset for the random number generator used in PDAF.

    Parameters
    ----------
    seedset_in : int
        Seedset index (1-20)

    Returns
    -------
    """

def set_smoother_ens(maxlag: int) -> typing.Tuple[np.ndarray, int]:
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
