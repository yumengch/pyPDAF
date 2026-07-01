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
    c__pdaf_get_obsmemberid(&memberid)

    return memberid

def get_seed():
    """get_seed() -> np.ndarray

    Return PDAF's current four-integer random seed vector.

    This is an alias of :func:`get_seedvec`. The returned seed can be stored
    and later passed to :func:`pyPDAF.PDAF.set_seed` or
    :func:`pyPDAF.PDAF.set_seedvec` to continue a reproducible PDAF random
    sequence.

    Returns
    -------
    seedvec : ndarray[np.intc, ndim=1]
        Four-integer PDAF random seed vector.
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] seedvec_np = np.zeros((4), dtype=np.intc, order="F")
    cdef int [::1] seedvec = seedvec_np
    c__pdaf_get_seed(&seedvec[0])
    return seedvec_np

def get_seedvec():
    """get_seedvec() -> np.ndarray

    Return PDAF's current four-integer random seed vector.

    The seed vector follows the LAPACK random-number convention used by PDAF.
    In particular, the fourth entry is expected to be odd when the vector is
    supplied back to PDAF.

    Returns
    -------
    seedvec : ndarray[np.intc, ndim=1]
        Four-integer PDAF random seed vector.

    See Also
    --------
    pyPDAF.PDAF.set_seedvec : Store a seed vector for PDAF random numbers.
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] seedvec_np = np.zeros((4), dtype=np.intc, order="F")
    cdef int [::1] seedvec = seedvec_np
    c__pdaf_get_seedvec(&seedvec[0])
    return seedvec_np

def get_rndcount():
    """get_rndcount() -> int

    Return PDAF's random-number generation counter.

    PDAF increments this internal counter when its random-number helper
    routines are used. It is mainly useful for diagnostics and for checking
    that a reproducible run consumed the expected number of random draws.

    Returns
    -------
    rndcount : int
        Number of PDAF random-number generation calls counted so far.
    """
    cdef int rndcount
    c__pdaf_get_rndcount(&rndcount)
    return rndcount


def reset_fcst_flag():
    """reset_fcst_flag() -> int

    Return PDAF's forecast-time reset flag.

    PDAF uses this flag in the flexible parallelization mode when an
    assimilation routine ``PDAF_assimilate_*`` is used instead of a
    ``PDAF_put_state_*`` routine. A nonzero value tells the user-side model
    integration to reset its model time before continuing the forecast.

    Returns
    -------
    reset_fcst_flag : int
        Internal PDAF flag indicating whether the user-side model time should
        be reset.
    """
    cdef int reset_fcst_flag_value
    c__pdaf_reset_fcst_flag(&reset_fcst_flag_value)
    return reset_fcst_flag_value


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
    cdef int  maxlag
    cdef int  status
    c__pdaf_get_smootherens(<CFI_cdesc_t *> &sens_point_cfi, &maxlag, &status)

    cdef CFI_index_t sens_point_subscripts[3]
    sens_point_subscripts[0] = sens_point_cfi.dim[0].lower_bound
    sens_point_subscripts[1] = sens_point_cfi.dim[1].lower_bound
    sens_point_subscripts[2] = sens_point_cfi.dim[2].lower_bound
    cdef double * sens_point_ptr_np
    sens_point_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &sens_point_cfi, sens_point_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_point_np = np.asarray(<double [:sens_point_cfi.dim[0].extent:1,:sens_point_cfi.dim[1].extent,:sens_point_cfi.dim[2].extent]> sens_point_ptr_np, order="F")
    return sens_point_np, maxlag, status


