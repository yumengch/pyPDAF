# pylint: disable=unused-argument
"""Stub file for PDAFomi diag module
"""
from typing import Tuple
import numpy as np


def set_obs_diag(diag: int) -> None:
    """Activate or deactivate the observation diagnostics.

    By default, observation diagnostics are activated that stores
    additional information for diagnostics.
    However, as this functionality increases the required memory,
    it might be desirable to deactivate this functionality.

    This function is used deactivate the observation diagnostics.
    Once deactivated, one cannot use diagnostics in :mod:`pyPDAF.PDAFomi.diag`.
    It is also possible to re-activate the observation diagnostics at a later time.

    The function can be called by all processes, but it is sufficient to call it
    for those processes that handle observations, which usually are the filter processes.

    This function can be called after the initialization of PDAF in `pyPDAF.PDAF.init`.


    Parameters
    ----------
    diag : int
        Value for observation diagnostics mode
        - > 0: activates observation diagnostics
        - 0: deactivates observation diagnostics
    """

def diag_dimobs() -> np.ndarray:
    """Observation dimension for each observation type.

    This function is only useful after observation operators are performed.

    Returns
    -------
    dim_obs: np.ndarray
        Observation dimension for each observation type. shape: (n_obs,)
    """

def diag_get_hx(id_obs: int) -> Tuple[int, np.ndarray]:
    """Observed ensemble for given observation type.

    This function is only useful after observation operators are performed.

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

def diag_get_hxmean(id_obs: int) -> Tuple[int, np.ndarray]:
    """Observed ensemble mean for given observation type.

    This function is only useful after observation operators are performed.

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

def diag_get_ivar(id_obs: int) -> Tuple[int, np.ndarray]:
    """Inverse of observation error variance for given observation type.

    This function is only useful after observation operators are performed.

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

def diag_get_obs(id_obs: int) -> Tuple[int, int, np.ndarray, np.ndarray]:
    """Observation vector and corresponding coordinates for specified observation type.

    This function is only useful after observation operators are performed.

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

def diag_nobstypes(nobs: int) -> int:
    """The number of observation types that are active in an assimilation run.

    This function is only useful after observation operators are performed.

    Parameters
    ----------
    nobs : int
        Number of observation types (input can be an arbitrary number)

    Returns
    -------
    nobs : int
        Number of observation types
    """

def diag_obs_rmsd(nobs: int, verbose: int) -> Tuple[int, np.ndarray]:
    """Root mean squared distance between observation and obseved model state
    for each observation type.

    This function is only useful after observation operators are performed.

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

def diag_stats(nobs: int, verbose: int) -> Tuple[int, np.ndarray]:
    """A selection of 6 statistics comparing the observations and the
    observed ensemble mean for each observation type.

    This function is only useful after observation operators are performed.

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

def diag_rmsd(nobs: int, verbose: int) -> Tuple[int, np.ndarray]:
    """Compute RMSD between observations and observed ensemble means.

    For each active OMI observation type, PDAF-OMI compares the observation
    vector ``y`` with the observed ensemble mean ``Hx`` and returns the root
    mean squared deviation
    ``sqrt(mean((y - Hx)**2))`` over all observations of that type. If only
    the observed ensemble ``H(X_i)`` was stored for diagnostics, PDAF-OMI first
    computes ensemble mean.

    This function is only useful after observation operators are performed.

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

def diag_diffstats(nobs: int, verbose: int) -> Tuple[int, np.ndarray]:
    """Compare observations with observed ensemble means by observation type.

    This function is the same as :func:`diag_stats`.

    This function is only useful after observation operators are performed.

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

def diag_crps(nobs: int, perturb: int, verbose: int) -> Tuple[int, np.ndarray]:
    """Compute CRPS diagnostics between observations and observed ensembles.

    For each active OMI observation type, this routine compares the observed
    ensemble ``H(X_i)`` with the observation vector. It returns CRPS and the
    Hersbach (2000) decomposition into reliability, potential CRPS, and
    uncertainty. CRPS is a probabilistic score: smaller values indicate that
    the ensemble distribution is closer to the verifying observation.

    If ``perturb=0``, PDAF-OMI uses the stored observations directly. For any
    nonzero value, PDAF-OMI creates perturbed observations by drawing Gaussian
    noise with variance from the stored inverse observation variance and adds
    it to the observations before computing CRPS.

    This function is only useful after observation operators are performed.

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

def diag_omit_by_inno():
    """Set omitted observations with high observation error for diagnostics only.

    The high observation errors are only used for observation diagnostics in
    :func:`diag_get_ivar` and :func:`diag_crps`.

    This function is only useful after observation operators are performed.
    """
    c__pdafomi_diag_omit_by_inno()