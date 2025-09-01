# pylint: disable=unused-argument
import numpy as np
from typing import Tuple

def diag_ensmean(
    dim: int,
    dim_ens: int,
    state: np.ndarray,
    ens: np.ndarray
) -> Tuple[np.ndarray, int]:
    r"""Compute the ensemble mean of the state ensemble.

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
) -> Tuple[np.ndarray, float, int]:
    r"""Compute ensemble standard deviation and ensemble mean without MPI.

    Parameters
    ----------
    dim : int
        state dimension
    dim_ens : int
        Ensemble size
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)
    do_mean : int
        Whether to compute ensemble mean

    Returns
    -------
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    stddev : double
        Standard deviation of ensemble
    status : int
        Status flag (0=success)
    """

def diag_stddev(
    dim_p: int,
    dim_ens: int,
    state_p: np.ndarray,
    ens_p: np.ndarray,
    do_mean: int,
    comm_filter: int
) -> Tuple[np.ndarray, float, int]:
    r"""Compute ensemble standard deviation and ensemble mean.

    Parameters
    ----------
    dim_p : int
        state dimension
    dim_ens : int
        Ensemble size
    state_p : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim_p, dim_ens)
    do_mean : int
        Whether to compute ensemble mean
    comm_filter : int
        Filter communicator

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim_p)
    stddev_g : double
        Global mean standard deviation of ensemble
    status : int
        Status flag (0=success)
    """

def diag_variance_nompi(
    dim: int,
    dim_ens: int,
    state: np.ndarray,
    ens: np.ndarray,
    do_mean: int,
    do_stddev: int
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    r"""Compute ensemble variance/standard deviation and mean without MPI.

    Parameters
    ----------
    dim : int
        state dimension
    dim_ens : int
        Ensemble size
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)
    do_mean : int
        Whether to compute ensemble mean
    do_stddev : int
        Whether to compute the ensemble mean standard deviation

    Returns
    -------
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    variance : ndarray[np.float64, ndim=1]
        Variance state vector
        Array shape: (dim)
    stddev : double
        Standard deviation of ensemble
    status : int
        Status flag (0=success)
    """

def diag_variance(
    dim_p: int,
    dim_ens: int,
    state_p: np.ndarray,
    ens_p: np.ndarray,
    do_mean: int,
    do_stddev: int,
    comm_filter: int
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    r"""Compute ensemble variance/standard deviation and mean.

    Parameters
    ----------
    dim_p : int
        state dimension
    dim_ens : int
        Ensemble size
    state_p : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim_p, dim_ens)
    do_mean : int
        Whether to compute ensemble mean
    do_stddev : int
        Whether to compute the ensemble mean standard deviation
    comm_filter : int
        Filter communicator

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim_p)
    variance_p : ndarray[np.float64, ndim=1]
        Variance state vector
        Array shape: (dim_p)
    stddev_g : double
        Global standard deviation of ensemble
    status : int
        Status flag (0=success)
    """

def diag_rmsd_nompi(
    dim_p: int,
    statea_p: np.ndarray,
    stateb_p: np.ndarray
) -> Tuple[float, int]:
    r"""Compute the root mean squared distance between two vectors without MPI.

    Parameters
    ----------
    dim_p : int
        state dimension
    statea_p : ndarray[np.float64, ndim=1]
        State vector A
        Array shape: (dim_p)
    stateb_p : ndarray[np.float64, ndim=1]
        State vector B
        Array shape: (dim_p)

    Returns
    -------
    rmsd_p : double
        RSMD
    status : int
        Status flag (0=success)
    """

def diag_rmsd(
    dim_p: int,
    statea_p: np.ndarray,
    stateb_p: np.ndarray,
    comm_filter: int
) -> Tuple[float, int]:
    r"""Compute the root mean squared distance between two vectors.

    Parameters
    ----------
    dim_p : int
        state dimension
    statea_p : ndarray[np.float64, ndim=1]
        State vector A
        Array shape: (dim_p)
    stateb_p : ndarray[np.float64, ndim=1]
        State vector B
        Array shape: (dim_p)
    comm_filter : int
        Filter communicator

    Returns
    -------
    rmsd_g : double
        Global RSMD
    status : int
        Status flag (0=success)
    """

def diag_crps(
    dim_p: int,
    dim_ens: int,
    element: int,
    oens: np.ndarray,
    obs: np.ndarray
) -> Tuple[float, float, float, float, int]:
    r"""Obtain a continuous rank probability score for an ensemble.

    The implementation is based on [1]_.

    References
    ----------
    .. [1] Hersbach, H. (2000),
           Decomposition of the Continuous Ranked Probability
           Score for
           Ensemble Prediction Systems,
           Wea. Forecasting, 15, 559–570,
           doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2

    Parameters
    ----------
    dim_p: int
        Dimension of state vector
    dim_ens: int
        Ensemble size
    element : int
        ID of element to be used. If element=0, mean values over all elements are computed
    oens : ndarray[tuple[dim, dim_ens, ...], np.float64]
        State ensemble. shape: (dim_p, dim_ens)
    obs : ndarray[tuple[dim, ...], np.float64]
        State ensemble. shape: (dim_p)

    Returns
    -------
    CRPS : float
        CRPS
    reli : float
        Reliability
    resol : float
        resolution
    uncert : float
        uncertainty
    status : int
        Status flag (0=success)
    """

def diag_crps_nompi(
    dim: int,
    dim_ens: int,
    element: int,
    oens: np.ndarray,
    obs: np.ndarray
) -> Tuple[float, float, float, float, int]:
    r"""Obtain a continuous rank probability score for an ensemble without MPI.

    The implementation is based on [1]_.

    References
    ----------
    .. [1] Hersbach, H. (2000),
           Decomposition of the Continuous Ranked Probability
           Score for
           Ensemble Prediction Systems,
           Wea. Forecasting, 15, 559–570,
           doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2

    Parameters
    ----------
    dim_p: int
        Dimension of state vector
    dim_ens: int
        Ensemble size
    element : int
        ID of element to be used. If element=0, mean values over all elements are computed
    oens : ndarray[tuple[dim, dim_ens, ...], np.float64]
        State ensemble. shape: (dim_p, dim_ens)
    obs : ndarray[tuple[dim, ...], np.float64]
        State ensemble. shape: (dim_p)

    Returns
    -------
    CRPS : float
        CRPS
    reli : float
        Reliability
    resol : float
        resolution
    uncert : float
        uncertainty
    status : int
        Status flag (0=success)
    """

def diag_effsample(
    dim_sample: int,
    weights: np.ndarray
) -> float:
    r"""Calculating the effective sample size of a particle filter.

    Based on [1]_, it is defined as the
    inverse of the sum of the squared particle filter weights:
    :math:`N_{eff} = \frac{1}{\sum_{i=1}^{N} w_i^2}`
    where :math:`w_i` is the weight of particle with index i.
    and :math:`N` is the number of particles.

    If the :math:`N_{eff}=N`, all weights are identical,
    and the filter has no influence on the analysis.
    If :math:`N_{eff}=0`, the filter is collapsed.

    This is typically called during the analysis step
    of a particle filter,
    e.g. in the analysis step of NETF and LNETF.

    References
    ----------
    .. [1] Doucet, A., de Freitas, N., Gordon, N. (2001).
           An Introduction to Sequential Monte Carlo Methods.
           In: Doucet, A., de Freitas, N., Gordon, N. (eds)
           Sequential Monte Carlo Methods in Practice.
           Statistics for Engineering and Information Science.
           Springer, New York, NY.
           https://doi.org/10.1007/978-1-4757-3437-9_1

    Parameters
    ----------
    dim_sample : int
        Sample size
    weights : ndarray[np.float64, ndim=1]
        Weights of the samples
        Array shape: (dim_sample)

    Returns
    -------
    n_eff : double
        Effecfive sample size
    """

def diag_ensstats(
    dim: int,
    dim_ens: int,
    element: int,
    state: np.ndarray,
    ens: np.ndarray
) -> Tuple[float, float, int]:
    r"""Computing the skewness and kurtosis of
    the ensemble of a given element of the state vector.

    The definition used for kurtosis follows that used by [1]_.

    References
    ----------
    .. [1] Lawson, W. G., & Hansen, J. A. (2004).
           Implications of stochastic and deterministic
           filters as ensemble-based
           data assimilation methods in varying regimes
           of error growth.
           Monthly weather review, 132(8), 1966-1981.

    Parameters
    ----------
    dim : int
        PE-local state dimension
    dim_ens : int
        Ensemble size
    element : int
        ID of element to be used
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)

    Returns
    -------
    skewness : double
        Skewness of ensemble
    kurtosis : double
        Kurtosis of ensemble
    status : int
        Status flag (0=success)
    """

def diag_compute_moments(
    dim_p: int,
    dim_ens: int,
    ens: np.ndarray,
    kmax: int,
    bias: int
) -> np.ndarray:
    r"""Computes the mean, the unbiased variance, the unbiased skewness,
    and the unbiased excess kurtosis from an ensemble.

    Parameters
    ----------
    dim_p : int
        local size of the state
    dim_ens : int
        number of ensemble members/samples
    ens : ndarray[np.float64, ndim=2]
        ensemble matrix
        Array shape: (dim_p,dim_ens)
    kmax : int
        maximum order of central moment that is computed, maximum is 4
    bias : int
        if 0 bias correction is applied (default)
        otherwise, no bias correction

    Returns
    -------
    moments : ndarray[np.float64, ndim=2]
        The columns contain the moments of the ensemble
        - column 0: mean
        - column 1: variance
        - column 2: skewness
        - column 3: excess kurtosis
        Array shape: (dim_p, kmax)
    """

def diag_histogram(
    ncall: int,
    dim: int,
    dim_ens: int,
    element: int,
    state: np.ndarray,
    ens: np.ndarray,
    hist: np.ndarray
) -> Tuple[np.ndarray, float, int]:
    r"""Computing the rank histogram of an ensemble.

    A rank histogram is used to diagnose
    the reliability of the ensemble [1]_.
    A perfectly reliable ensemble should have
    a uniform rank histogram.

    The function can be called in the
    pre/poststep routine of PDAF
    both before and after the analysis step
    to collect the histogram information.

    References
    ----------
    .. [1] Hamill, T. M. (2001).
           Interpretation of rank histograms
           for verifying ensemble forecasts.
           Monthly Weather Review, 129(3), 550-560.

    Parameters
    ----------
    ncall : int
        Number of calls to routine
    dim : int
        State dimension
    dim_ens : int
        Ensemble size
    element : int
        Element of vector used for histogram
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)
    hist : ndarray[np.intc, ndim=1]
        Histogram about the state
        Array shape: (dim_ens+1)

    Returns
    -------
    hist : ndarray[np.intc, ndim=1]
        Histogram about the state
        Array shape: (dim_ens+1)
    delta : double
        deviation measure from flat histogram
    status : int
        Status flag (0=success)
    """

def diag_reliability_budget(
    n_times: int,
    dim_ens: int,
    dim_p: int,
    ens_p: np.ndarray,
    obsvar: np.ndarray,
    obs_p: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Compute ensemble reliability budget

    The diagnostics derives a balance relationship that decomposes the
    departure between the ensemble mean and observations:
    Depar^2 = Bias^2 + EnsVar + ObsUnc^2 + Residual
    under the assumption of a perfectly reliable ensemble [1]_. When the
    residual term is not small, one can identify the sources of
    problematic ensemble representation of the uncertainty by looking at
    each terms. One of the benefits of
    the reliability budget diagnostics is that it can identify spatially
    local issues in ensemble perturbations.

    In this subroutine, the budget array returns each term of the
    reliability budget is given at each given time step. The returned
    array can be used for significant t-test where each time step is used
    as a sample from the population. The actual budget is the mean over
    time steps except for Bias^2 term, which is given separately. This
    is because the unsquared bias is used in the t-test in the original
    paper.

    The subroutine does not comes with a t-test to ensure that
    the reliability budget is statistically significant as done
    in the original paper. However, this can be achieved with external
    libraries or software if the t-test is deemed important. The t-test
    should account for the temporal autocorrelation between time steps
    in the given trajectory.

    The diagnostics requires a trajectory of ensemble and observations,
    a trajectory of an ensemle of observation error variances. The
    ensemble of observation error variances can be the square of
    observation errors sampled from observation error distribution as
    done in stochastic ensemble Kalman filter. In deterministic ensemble
    systems, the ensemble of observation error variances can be the
    observation error variances where each ensemble member has the same
    value.

    References
    ----------
    .. [1] Rodwell, M. J., Lang, S. T. K., Ingleby, N. B., Bormann, N., Holm,
    E., Rabier, F., ... & Yamaguchi, M. (2016). Reliability in ensemble data assimilation.
    Quarterly Journal of the Royal Meteorological Society, 142(694), 443-454.

    Parameters
    ----------
    n_times : int
        Number of time steps
    dim_ens : int
        Number of ensemble members
    dim_p : int
        Dimension of the state vector
    ens_p : ndarray[np.float64, ndim=3]
        Ensemble matrix over times
        Array shape: (dim_p, dim_ens, n_times)
    obsvar : ndarray[np.float64, ndim=3]
        Squared observation error/variance at n_times
        Array shape: (dim_p, dim_ens, n_times)
    obs_p : ndarray[np.float64, ndim=2]
        Observation vector
        Array shape: (dim_p, n_times)

    Returns
    -------
    budget : ndarray[np.float64, ndim=3]
        Budget term for a single time step
        Array shape: (dim_p, n_times, 5)
    bias_2 : ndarray[np.float64, ndim=1]
        bias^2 uses
        Array shape: (dim_p)
    """