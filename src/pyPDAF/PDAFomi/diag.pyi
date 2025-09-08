# pylint: disable=unused-argument
"""Stub file for PDAFomi diag module
"""
from typing import Tuple
import numpy as np


def diag_dimobs() -> np.ndarray:
    """Observation dimension for each observation type.

    Returns
    -------
    dim_obs: np.ndarray
        Observation dimension for each observation type. shape: (n_obs,)
    """

def diag_get_hx(id_obs: int) -> Tuple[int, np.ndarray]:
    """Observed ensemble for given observation type.

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
