import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3


def diag_dimobs():
    """diag_dimobs() -> np.ndarray

    Observation dimension for each observation type.

    Returns
    -------
    dim_obs: np.ndarray
        Observation dimension for each observation type. shape: (n_obs,)
    """
    cdef CFI_cdesc_rank1 dim_obs_ptr_cfi
    c__pdafomi_diag_dimobs(<CFI_cdesc_t *> &dim_obs_ptr_cfi)

    cdef CFI_index_t dim_obs_ptr_subscripts[1]
    dim_obs_ptr_subscripts[0] = dim_obs_ptr_cfi.dim[0].lower_bound
    cdef int * dim_obs_ptr_ptr_np
    dim_obs_ptr_ptr_np = <int *>CFI_address(<CFI_cdesc_t *> &dim_obs_ptr_cfi, dim_obs_ptr_subscripts)
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] dim_obs_ptr_np = np.asarray(<int [:dim_obs_ptr_cfi.dim[0].extent:1]> dim_obs_ptr_ptr_np, order="F")
    return dim_obs_ptr_np


def diag_get_hx(int  id_obs):
    """diag_get_hx(id_obs: int) -> Tuple[int, np.ndarray]

    Observed ensemble for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    hx_p_ptr : ndarray[np.float64, ndim=2]
        Pointer to observed ensemble mean
        Array shape: (dim_obs_p_diag, dim_ens)
    """
    cdef CFI_cdesc_rank2 hx_p_ptr_cfi
    cdef int  dim_obs_diag
    c__pdafomi_diag_get_hx(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &hx_p_ptr_cfi)

    cdef CFI_index_t hx_p_ptr_subscripts[2]
    hx_p_ptr_subscripts[0] = hx_p_ptr_cfi.dim[0].lower_bound
    hx_p_ptr_subscripts[1] = hx_p_ptr_cfi.dim[1].lower_bound
    cdef double * hx_p_ptr_ptr_np
    hx_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &hx_p_ptr_cfi, hx_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hx_p_ptr_np = np.asarray(<double [:hx_p_ptr_cfi.dim[0].extent:1,:hx_p_ptr_cfi.dim[1].extent]> hx_p_ptr_ptr_np, order="F")
    return dim_obs_diag, hx_p_ptr_np


def diag_get_hxmean(int  id_obs):
    """diag_get_hxmean(id_obs: int) -> Tuple[int, np.ndarray]

    Observed ensemble mean for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    hxmean_p_ptr : ndarray[np.float64, ndim=1]
        Pointer to observed ensemble mean
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 hxmean_p_ptr_cfi
    cdef int  dim_obs_diag
    c__pdafomi_diag_get_hxmean(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &hxmean_p_ptr_cfi)

    cdef CFI_index_t hxmean_p_ptr_subscripts[1]
    hxmean_p_ptr_subscripts[0] = hxmean_p_ptr_cfi.dim[0].lower_bound
    cdef double * hxmean_p_ptr_ptr_np
    hxmean_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &hxmean_p_ptr_cfi, hxmean_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxmean_p_ptr_np = np.asarray(<double [:hxmean_p_ptr_cfi.dim[0].extent:1]> hxmean_p_ptr_ptr_np, order="F")
    return dim_obs_diag, hxmean_p_ptr_np


def diag_get_ivar(int  id_obs):
    """diag_get_ivar(id_obs: int) -> Tuple[int, np.ndarray]

    Inverse of observation error variance for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    ivar_ptr : ndarray[np.float64, ndim=1]
        Pointer to inverse observation error variances
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 ivar_ptr_cfi
    cdef int  dim_obs_diag
    c__pdafomi_diag_get_ivar(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &ivar_ptr_cfi)

    cdef CFI_index_t ivar_ptr_subscripts[1]
    ivar_ptr_subscripts[0] = ivar_ptr_cfi.dim[0].lower_bound
    cdef double * ivar_ptr_ptr_np
    ivar_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &ivar_ptr_cfi, ivar_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] ivar_ptr_np = np.asarray(<double [:ivar_ptr_cfi.dim[0].extent:1]> ivar_ptr_ptr_np, order="F")
    return dim_obs_diag, ivar_ptr_np


def diag_get_obs(int  id_obs):
    """diag_get_obs(id_obs: int) -> Tuple[int, int, np.ndarray, np.ndarray]

    Observation vector and corresponding coordinates for specified observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    ncoord : int
        Number of observation dimensions
    obs_p_ptr : ndarray[np.float64, ndim=1]
        Pointer to observation vector
        Array shape: (:)
    ocoord_p_ptr : ndarray[np.float64, ndim=2]
        Pointer to Coordinate array
        Array shape: (:,:)
    """
    cdef CFI_cdesc_rank1 obs_p_ptr_cfi
    cdef CFI_cdesc_rank2 ocoord_p_ptr_cfi
    cdef int  dim_obs_diag
    cdef int  ncoord
    c__pdafomi_diag_get_obs(&id_obs, &dim_obs_diag, &ncoord,
                                <CFI_cdesc_t *> &obs_p_ptr_cfi, <CFI_cdesc_t *> &ocoord_p_ptr_cfi)

    cdef CFI_index_t obs_p_ptr_subscripts[1]
    obs_p_ptr_subscripts[0] = obs_p_ptr_cfi.dim[0].lower_bound
    cdef double * obs_p_ptr_ptr_np
    obs_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &obs_p_ptr_cfi, obs_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_p_ptr_np = np.asarray(<double [:obs_p_ptr_cfi.dim[0].extent:1]> obs_p_ptr_ptr_np, order="F")
    cdef CFI_index_t ocoord_p_ptr_subscripts[2]
    ocoord_p_ptr_subscripts[0] = ocoord_p_ptr_cfi.dim[0].lower_bound
    ocoord_p_ptr_subscripts[1] = ocoord_p_ptr_cfi.dim[1].lower_bound
    cdef double * ocoord_p_ptr_ptr_np
    ocoord_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &ocoord_p_ptr_cfi, ocoord_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ocoord_p_ptr_np = np.asarray(<double [:ocoord_p_ptr_cfi.dim[0].extent:1,:ocoord_p_ptr_cfi.dim[1].extent]> ocoord_p_ptr_ptr_np, order="F")
    return dim_obs_diag, ncoord, obs_p_ptr_np, ocoord_p_ptr_np


def diag_nobstypes(int  nobs):
    """diag_nobstypes(nobs: int) -> int

    The number of observation types that are active in an assimilation run.

    Parameters
    ----------
    nobs : int
        Number of observation types (input can be an arbitrary number)

    Returns
    -------
    nobs : int
        Number of observation types
    """
    c__pdafomi_diag_nobstypes(&nobs)

    return nobs


def diag_obs_rmsd(int  nobs, int  verbose):
    """diag_obs_rmsd(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    Root mean squared distance between observation and obseved model state
    for each observation type.

    Parameters
    ----------
    nobs : int
        Number of observation types
    verbose : int
        Verbosity flag

    Returns
    -------
    nobs : int
        Number of observation types
    rmsd_pointer : ndarray[np.float64, ndim=1]
        Vector of RMSD values
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 rmsd_pointer_cfi
    c__pdafomi_diag_obs_rmsd(&nobs, <CFI_cdesc_t *> &rmsd_pointer_cfi, &verbose)

    cdef CFI_index_t rmsd_pointer_subscripts[1]
    rmsd_pointer_subscripts[0] = rmsd_pointer_cfi.dim[0].lower_bound
    cdef double * rmsd_pointer_ptr_np
    rmsd_pointer_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &rmsd_pointer_cfi, rmsd_pointer_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] rmsd_pointer_np = np.asarray(<double [:rmsd_pointer_cfi.dim[0].extent:1]> rmsd_pointer_ptr_np, order="F")
    return nobs, rmsd_pointer_np


def diag_stats(int  nobs, int  verbose):
    """diag_stats(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    A selection of 6 statistics comparing the observations and the
    observed ensemble mean for each observation type.

    Parameters
    ----------
    nobs : int
        Number of observation types
    verbose : int
        Verbosity flag

    Returns
    -------
    nobs : int
        Number of observation types
    obsstats_ptr : ndarray[np.float64, ndim=2]
        Array of observation statistics
        Included statistics are:
        - (1,:) correlations between observation and observed ensemble mean
        - (2,:) centered RMS difference between observation and observed ensemble mean
        - (3,:) mean bias (observation minus observed ensemble mean)
        - (4,:) mean absolute difference between observation and observed ensemble mean
        - (5,:) variance of observations
        - (6,:) variance of observed ensemble mean
        Array shape: (:,:)
    """
    cdef CFI_cdesc_rank2 obsstats_ptr_cfi
    c__pdafomi_diag_stats(&nobs, <CFI_cdesc_t *> &obsstats_ptr_cfi, &verbose)

    cdef CFI_index_t obsstats_ptr_subscripts[2]
    obsstats_ptr_subscripts[0] = obsstats_ptr_cfi.dim[0].lower_bound
    obsstats_ptr_subscripts[1] = obsstats_ptr_cfi.dim[1].lower_bound
    cdef double * obsstats_ptr_ptr_np
    obsstats_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &obsstats_ptr_cfi, obsstats_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] obsstats_ptr_np = np.asarray(<double [:obsstats_ptr_cfi.dim[0].extent:1,:obsstats_ptr_cfi.dim[1].extent]> obsstats_ptr_ptr_np, order="F")
    return nobs, obsstats_ptr_np

def diag_rmsd(int nobs, int verbose):
    """diag_rmsd(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    Compute RMSD between observations and observed ensemble means.

    For each active OMI observation type, PDAF-OMI compares the observation
    vector ``y`` with the observed ensemble mean ``Hx`` and returns the root
    mean squared deviation
    ``sqrt(mean((y - Hx)**2))`` over all observations of that type. If only
    the observed ensemble ``H(X_i)`` was stored for diagnostics, PDAF-OMI first
    computes ``Hx`` from that ensemble.

    The routine only returns diagnostics when OMI diagnostic observation data
    are available. If no observation diagnostics are active, ``nobs`` is set
    to ``0``.

    Parameters
    ----------
    nobs : int
        Input/output number of OMI observation types. The input value can be
        arbitrary; PDAF-OMI overwrites it with the number of active
        observation types for which diagnostics are returned.
    verbose : int
        Verbosity flag. If greater than zero, PDAF-OMI prints the RMSD table.

    Returns
    -------
    nobs : int
        Number of observation-type entries returned by PDAF-OMI.
    rmsd : ndarray[np.float64, ndim=1]
        PDAF-owned RMSD values. ``rmsd[i]`` is the RMSD for observation type
        ``i + 1``.
    """
    cdef CFI_cdesc_rank1 rmsd_pointer_cfi
    c__pdafomi_diag_rmsd(&nobs, <CFI_cdesc_t *> &rmsd_pointer_cfi, &verbose)
    cdef CFI_index_t rmsd_pointer_subscripts[1]
    rmsd_pointer_subscripts[0] = rmsd_pointer_cfi.dim[0].lower_bound
    cdef double * rmsd_pointer_ptr_np
    rmsd_pointer_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &rmsd_pointer_cfi, rmsd_pointer_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] rmsd_pointer_np = np.asarray(<double [:rmsd_pointer_cfi.dim[0].extent:1]> rmsd_pointer_ptr_np, order="F")
    return nobs, rmsd_pointer_np

def diag_diffstats(int nobs, int verbose):
    """diag_diffstats(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    Compare observations with observed ensemble means by observation type.

    This is the OMI observation-diagnostic counterpart of
    :func:`pyPDAF.PDAF.diag_diffstats`. For each active observation type it
    compares the observation vector ``y`` with the observed ensemble mean
    ``Hx``. The returned matrix is organized as ``stats[statistic, obs_type]``:
    rows select the statistic and columns select the observation type.

    These diagnostics are useful for Taylor diagrams and for separating a
    mean offset from pattern or variability errors. Correlation and standard
    deviations describe the centered variability of ``y`` and ``Hx``; centered
    RMSD describes the RMS mismatch after removing the means; bias and mean
    absolute deviation describe the raw difference ``y - Hx``.

    Parameters
    ----------
    nobs : int
        Input/output number of OMI observation types. The input value can be
        arbitrary; PDAF-OMI overwrites it with the number of active
        observation types for which diagnostics are returned.
    verbose : int
        Verbosity flag. If greater than zero, PDAF-OMI prints one row per
        observation type.

    Returns
    -------
    nobs : int
        Number of observation-type entries returned by PDAF-OMI.
    stats : ndarray[np.float64, ndim=2]
        PDAF-owned statistics matrix with shape ``(6, nobs)``. The rows are:

        ``stats[0, :]``
            Pearson correlation between observation anomalies and observed
            ensemble-mean anomalies.
        ``stats[1, :]``
            Centered RMS deviation between ``y`` and ``Hx``.
        ``stats[2, :]``
            Bias ``mean(y) - mean(Hx)``.
        ``stats[3, :]``
            Mean absolute deviation ``mean(abs(y - Hx))``.
        ``stats[4, :]``
            Standard deviation of the observations ``y``.
        ``stats[5, :]``
            Standard deviation of the observed ensemble mean ``Hx``.
    """
    cdef CFI_cdesc_rank2 obsstats_ptr_cfi
    c__pdafomi_diag_diffstats(&nobs, <CFI_cdesc_t *> &obsstats_ptr_cfi, &verbose)
    cdef CFI_index_t obsstats_ptr_subscripts[2]
    obsstats_ptr_subscripts[0] = obsstats_ptr_cfi.dim[0].lower_bound
    obsstats_ptr_subscripts[1] = obsstats_ptr_cfi.dim[1].lower_bound
    cdef double * obsstats_ptr_ptr_np
    obsstats_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &obsstats_ptr_cfi, obsstats_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] obsstats_ptr_np = np.asarray(<double [:obsstats_ptr_cfi.dim[0].extent:1,:obsstats_ptr_cfi.dim[1].extent]> obsstats_ptr_ptr_np, order="F")
    return nobs, obsstats_ptr_np

def diag_crps(int nobs, int perturb, int verbose):
    """diag_crps(nobs: int, perturb: int, verbose: int) -> Tuple[int, np.ndarray]

    Compute CRPS diagnostics between observations and observed ensembles.

    For each active OMI observation type, this routine compares the observed
    ensemble ``H(X_i)`` with the observation vector. It returns CRPS and the
    Hersbach (2000) decomposition into reliability, potential CRPS, and
    uncertainty. CRPS is a probabilistic score: smaller values indicate that
    the ensemble distribution is closer to the verifying observation.

    If ``perturb=0``, PDAF-OMI uses the stored observations directly. For any
    nonzero value, PDAF-OMI creates perturbed observations by drawing Gaussian
    noise with variance from the stored inverse observation variance and adds
    it to the observations before computing CRPS.

    Parameters
    ----------
    nobs : int
        Input/output number of OMI observation types. The input value can be
        arbitrary; PDAF-OMI overwrites it with the number of active
        observation types for which diagnostics are returned.
    perturb : int
        Perturbation option. ``0`` uses the stored observations. Any nonzero
        value uses perturbed observations with Gaussian observation-error
        noise.
    verbose : int
        Verbosity flag. If greater than zero, PDAF-OMI prints one row per
        observation type.

    Returns
    -------
    nobs : int
        Number of observation-type entries returned by PDAF-OMI.
    crps : ndarray[np.float64, ndim=2]
        PDAF-owned CRPS diagnostic matrix with shape ``(4, nobs)``. The rows
        are CRPS, reliability, potential CRPS, and uncertainty.

    References
    ----------
    Hersbach, H. (2000). Decomposition of the continuous ranked probability
    score for ensemble prediction systems. Weather and Forecasting, 15,
    559-570.
    """
    cdef CFI_cdesc_rank2 crps_pointer_cfi
    c__pdafomi_diag_crps(&nobs, <CFI_cdesc_t *> &crps_pointer_cfi, &perturb, &verbose)
    cdef CFI_index_t crps_pointer_subscripts[2]
    crps_pointer_subscripts[0] = crps_pointer_cfi.dim[0].lower_bound
    crps_pointer_subscripts[1] = crps_pointer_cfi.dim[1].lower_bound
    cdef double * crps_pointer_ptr_np
    crps_pointer_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &crps_pointer_cfi, crps_pointer_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] crps_pointer_np = np.asarray(<double [:crps_pointer_cfi.dim[0].extent:1,:crps_pointer_cfi.dim[1].extent]> crps_pointer_ptr_np, order="F")
    return nobs, crps_pointer_np

