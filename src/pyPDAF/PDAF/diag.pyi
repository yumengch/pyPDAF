# pylint: disable=unused-argument
import numpy as np
from typing import Tuple

def diag_ensmean(
    dim: int,
    dim_ens: int,
    state: np.ndarray[Tuple[int], float],
    ens: np.ndarray[Tuple[int, int], float]
) -> Tuple[np.ndarray, int]:
    """
    Compute the ensemble mean of the state ensemble.

    Parameters
    ----------
    dim : int
        State dimension.
    dim_ens : int
        Ensemble size.
    state : ndarray[np.float64, ndim=1], shape (dim,)
        State vector.
    ens : ndarray[np.float64, ndim=2], shape (dim, dim_ens)
        State ensemble.

    Returns
    -------
    state : ndarray[np.float64, ndim=1], shape (dim,)
        State vector (ensemble mean).
    status : int
        Status flag (0=success).
    """

def diag_stddev_nompi(
    dim: int,
    dim_ens: int,
    state: np.ndarray,
    ens: np.ndarray,
    do_mean: int
) -> Tuple[np.ndarray, float, int]: ...

def diag_stddev(
    dim_p: int,
    dim_ens: int,
    state_p: np.ndarray,
    ens_p: np.ndarray,
    do_mean: int,
    comm_filter: int
) -> Tuple[np.ndarray, float, int]: ...

def diag_variance_nompi(
    dim: int,
    dim_ens: int,
    state: np.ndarray,
    ens: np.ndarray,
    do_mean: int,
    do_stddev: int
) -> Tuple[np.ndarray, np.ndarray, float, int]: ...

def diag_variance(
    dim_p: int,
    dim_ens: int,
    state_p: np.ndarray,
    ens_p: np.ndarray,
    do_mean: int,
    do_stddev: int,
    comm_filter: int
) -> Tuple[np.ndarray, np.ndarray, float, int]: ...

def diag_rmsd_nompi(
    dim_p: int,
    statea_p: np.ndarray,
    stateb_p: np.ndarray
) -> Tuple[float, int]: ...

def diag_rmsd(
    dim_p: int,
    statea_p: np.ndarray,
    stateb_p: np.ndarray,
    comm_filter: int
) -> Tuple[float, int]: ...

def diag_crps(
    dim_p: int,
    dim_ens: int,
    element: int,
    oens: np.ndarray,
    obs: np.ndarray
) -> Tuple[float, float, float, float, int]: ...

def diag_crps_mpi(
    dim_p: int,
    dim_ens: int,
    element: int,
    oens: np.ndarray,
    obs: np.ndarray,
    comm_filter: int,
    mype_filter: int,
    npes_filter: int
) -> Tuple[float, float, float, float, int]: ...

def diag_crps_nompi(
    dim: int,
    dim_ens: int,
    element: int,
    oens: np.ndarray,
    obs: np.ndarray
) -> Tuple[float, float, float, float, int]: ...

def diag_effsample(
    dim_sample: int,
    weights: np.ndarray
) -> float: ...

def diag_ensstats(
    dim: int,
    dim_ens: int,
    element: int,
    state: np.ndarray,
    ens: np.ndarray
) -> Tuple[float, float, int]: ...

def diag_compute_moments(
    dim_p: int,
    dim_ens: int,
    ens: np.ndarray,
    kmax: int,
    bias: int
) -> np.ndarray: ...

def diag_histogram(
    ncall: int,
    dim: int,
    dim_ens: int,
    element: int,
    state: np.ndarray,
    ens: np.ndarray,
    hist: np.ndarray
) -> Tuple[np.ndarray, float, int]: ...

def diag_reliability_budget(
    n_times: int,
    dim_ens: int,
    dim_p: int,
    ens_p: np.ndarray,
    obsvar: np.ndarray,
    obs_p: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]: ...