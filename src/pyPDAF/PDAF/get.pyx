import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t
from pyPDAF.cfi_binding cimport CFI_type_double
from pyPDAF.cfi_binding cimport CFI_cdesc_rank3

def get_assim_flag():
    r"""get_assim_flag() -> int

    Return the flag that
    indicates if the DA is performed in the last time step.
    It only works for online DA systems.


    Returns
    -------
    did_assim : int
        flag: (1) for assimilation; (0) else
    """
    cdef int  did_assim
    with nogil:
        c__pdaf_get_assim_flag(&did_assim)

    return did_assim


def get_localfilter():
    r"""get_localfilter() -> int

    Return whether a local filter is used.


    Returns
    -------
    lfilter : int
        whether the filter is domain-localized (1) or not (0)
        * 1 for local filters (including ENSRF/EAKF),
        * 0 for global filters (including LEnKF, which performs covariance localization)
    """
    cdef int  localfilter_out
    with nogil:
        c__pdaf_get_localfilter(&localfilter_out)

    return localfilter_out


def get_local_type():
    r"""get_local_type() -> int

    The routine returns the information on the localization type of the selected filter.

    With this one can distinguish filters using
    * domain localization (LESTKF, LETKF, LSEIK, LNETF),
    * covariance localization (LEnKF), or
    * covariance localization with observation handling like domain localization (ENSRF/EAKF).


    Returns
    -------
    localtype : int
        * (0) no localization; global filter
        * (1) domain localization (LESTKF, LETKF, LNETF, LSEIK)
        * (2) covariance localization (LEnKF)
        * (3) covariance loc. but observation handling like domain localization (ENSRF)
    """
    cdef int  localtype
    with nogil:
        c__pdaf_get_local_type(&localtype)

    return localtype


def get_memberid(int  memberid):
    """get_memberid(memberid: int) -> int

    Return the ensemble member id on the current process.

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
    """get_obsmemberid(memberid: int) -> int

    Return the ensemble member id
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


def get_smoother_ens():
    """get_smoother_ens() -> typing.Tuple[np.ndarray, int, int]

    Return the smoothed ensemble in earlier time steps.

    It is only used when the smoother options is used .

    Returns
    -------
    sens_point_np : ndarray[np.float64, ndim=3]
        Pointer to smoothed ensemble
    maxlag: int
        Maximum lag
    status: int
        Status flag
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


