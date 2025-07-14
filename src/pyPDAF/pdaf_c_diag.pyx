import sys
import numpy as np
cimport numpy as cnp
from . cimport pdaf_c_cb_interface as pdaf_cb
from .cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from .cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from .cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

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

def diag_ensmean(int  dim, int  dim_ens, double [::1] state, 
    double [::1,:] ens):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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

    Returns
    -------
    state : ndarray[np.float64, ndim=1]
        State vector
        Array shape: (dim)
    status : int 
        Status flag (0=success)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")
    cdef int  status
    with nogil:
        c__pdaf_diag_ensmean(&dim, &dim_ens, &state[0], &ens[0,0], &status)

    return state_np, status


def diag_stddev_nompi(int  dim, int  dim_ens, double [::1] state, 
    double [::1,:] ens, int  do_mean):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")
    cdef double  stddev
    cdef int  status
    with nogil:
        c__pdaf_diag_stddev_nompi(&dim, &dim_ens, &state[0], &ens[0,0], 
                                  &stddev, &do_mean, &status)

    return state_np, stddev, status


def diag_stddev(int  dim_p, int  dim_ens, double [::1] state_p, 
    double [::1,:] ens_p, int  do_mean, int  comm_filter):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef double  stddev_g
    cdef int  status
    with nogil:
        c__pdaf_diag_stddev(&dim_p, &dim_ens, &state_p[0], &ens_p[0,0], 
                            &stddev_g, &do_mean, &comm_filter, &status)

    return state_p_np, stddev_g, status


def diag_variance_nompi(int  dim, int  dim_ens, double [::1] state, 
    double [::1,:] ens, int  do_mean, int  do_stddev):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] variance_np = np.zeros((dim), dtype=np.float64, order="F")
    cdef double [::1] variance = variance_np
    cdef double  stddev
    cdef int  status
    with nogil:
        c__pdaf_diag_variance_nompi(&dim, &dim_ens, &state[0], &ens[0,0], 
                                    &variance[0], &stddev, &do_mean, 
                                    &do_stddev, &status)

    return state_np, variance_np, stddev, status


def diag_variance(int  dim_p, int  dim_ens, double [::1] state_p, 
    double [::1,:] ens_p, int  do_mean, int  do_stddev, int  comm_filter):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] variance_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] variance_p = variance_p_np
    cdef double  stddev_g
    cdef int  status
    with nogil:
        c__pdaf_diag_variance(&dim_p, &dim_ens, &state_p[0], &ens_p[0,0], 
                              &variance_p[0], &stddev_g, &do_mean, 
                              &do_stddev, &comm_filter, &status)

    return state_p_np, variance_p_np, stddev_g, status


def diag_rmsd_nompi(int  dim_p, double [::1] statea_p, double [::1] stateb_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef double  rmsd_p
    cdef int  status
    with nogil:
        c__pdaf_diag_rmsd_nompi(&dim_p, &statea_p[0], &stateb_p[0], 
                                &rmsd_p, &status)

    return rmsd_p, status


def diag_rmsd(int  dim_p, double [::1] statea_p, double [::1] stateb_p, 
    int  comm_filter):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef double  rmsd_g
    cdef int  status
    with nogil:
        c__pdaf_diag_rmsd(&dim_p, &statea_p[0], &stateb_p[0], &rmsd_g, 
                          &comm_filter, &status)

    return rmsd_g, status


def diag_crps(int  dim_p, int  dim_ens, int  element, double [::1,:] oens, 
    double [::1] obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int 
        PE-local state dimension
    dim_ens : int 
        Ensemble size
    element : int 
        index of element in full state vector
    oens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim_p, dim_ens)
    obs : ndarray[np.float64, ndim=1]
        Observation / truth
        Array shape: (dim_p)

    Returns
    -------
    crps : double 
        CRPS
    reli : double 
        Reliability
    pot_crps : double 
        potential CRPS
    uncert : double 
        uncertainty
    status : int 
        Status flag (0=success)
    """
    cdef double  crps
    cdef double  reli
    cdef double  pot_crps
    cdef double  uncert
    cdef int  status
    with nogil:
        c__pdaf_diag_crps(&dim_p, &dim_ens, &element, &oens[0,0], &obs[0], 
                          &crps, &reli, &pot_crps, &uncert, &status)

    return crps, reli, pot_crps, uncert, status


def diag_crps_mpi(int  dim_p, int  dim_ens, int  element, 
    double [::1,:] oens, double [::1] obs, int  comm_filter, 
    int  mype_filter, int  npes_filter):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int 
        PE-local state dimension
    dim_ens : int 
        Ensemble size
    element : int 
        index of element in full state vector
    oens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim_p, dim_ens)
    obs : ndarray[np.float64, ndim=1]
        Observation / truth
        Array shape: (dim_p)
    comm_filter : int 
        MPI communicator for filter
    mype_filter : int 
        rank of MPI communicator
    npes_filter : int 
        size of MPI communicator

    Returns
    -------
    crps : double 
        CRPS
    reli : double 
        Reliability
    pot_crps : double 
        potential CRPS
    uncert : double 
        uncertainty
    status : int 
        Status flag (0=success)
    """
    cdef double  crps
    cdef double  reli
    cdef double  pot_crps
    cdef double  uncert
    cdef int  status
    with nogil:
        c__pdaf_diag_crps_mpi(&dim_p, &dim_ens, &element, &oens[0,0], 
                              &obs[0], &comm_filter, &mype_filter, 
                              &npes_filter, &crps, &reli, &pot_crps, 
                              &uncert, &status)

    return crps, reli, pot_crps, uncert, status


def diag_crps_nompi(int  dim, int  dim_ens, int  element, 
    double [::1,:] oens, double [::1] obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int 
        PE-local state dimension
    dim_ens : int 
        Ensemble size
    element : int 
        ID of element to be used
    oens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)
    obs : ndarray[np.float64, ndim=1]
        State ensemble
        Array shape: (dim)

    Returns
    -------
    crps : double 
        CRPS
    reli : double 
        Reliability
    resol : double 
        resolution
    uncert : double 
        uncertainty
    status : int 
        Status flag (0=success)
    """
    cdef double  crps
    cdef double  reli
    cdef double  resol
    cdef double  uncert
    cdef int  status
    with nogil:
        c__pdaf_diag_crps_nompi(&dim, &dim_ens, &element, &oens[0,0], 
                                &obs[0], &crps, &reli, &resol, &uncert, &status)

    return crps, reli, resol, uncert, status


def diag_effsample(int  dim_sample, double [::1] weights):
    """Calculating the effective sample size of a particle filter.

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
    cdef double  n_eff
    with nogil:
        c__pdaf_diag_effsample(&dim_sample, &weights[0], &n_eff)

    return n_eff


def diag_ensstats(int  dim, int  dim_ens, int  element, double [::1] state, 
    double [::1,:] ens):
    """Computing the skewness and kurtosis of
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
    cdef double  skewness
    cdef double  kurtosis
    cdef int  status
    with nogil:
        c__pdaf_diag_ensstats(&dim, &dim_ens, &element, &state[0], 
                              &ens[0,0], &skewness, &kurtosis, &status)

    return skewness, kurtosis, status


def diag_compute_moments(int  dim_p, int  dim_ens, double [::1,:] ens, 
    int  kmax, int  bias):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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

    Returns
    -------
    moments : ndarray[np.float64, ndim=2]
        The columns contain the moments of the ensemble
        Array shape: (dim_p, kmax)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] moments_np = np.zeros((dim_p, kmax), dtype=np.float64, order="F")
    cdef double [::1,:] moments = moments_np
    with nogil:
        c__pdaf_diag_compute_moments(&dim_p, &dim_ens, &ens[0,0], &kmax, 
                                     &moments[0,0], &bias)

    return moments_np


def diag_histogram(int  ncall, int  dim, int  dim_ens, int  element, 
    double [::1] state, double [::1,:] ens, int [::1] hist):
    """Computing the rank histogram of an ensemble.

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
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hist_np = np.asarray(hist, dtype=np.intc, order="F")
    cdef double  delta
    cdef int  status
    with nogil:
        c__pdaf_diag_histogram(&ncall, &dim, &dim_ens, &element, &state[0], 
                               &ens[0,0], &hist[0], &delta, &status)

    return hist_np, delta, status


def diag_reliability_budget(int  n_times, int  dim_ens, int  dim_p, 
    double [::1,:,:] ens_p, double [::1,:,:] obsvar, double [::1,:] obs_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] budget_np = np.zeros((dim_p, n_times, 5), dtype=np.float64, order="F")
    cdef double [::1,:,:] budget = budget_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] bias_2_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] bias_2 = bias_2_np
    with nogil:
        c__pdaf_diag_reliability_budget(&n_times, &dim_ens, &dim_p, 
                                        &ens_p[0,0,0], &obsvar[0,0,0], 
                                        &obs_p[0,0], &budget[0,0,0], &bias_2[0])

    return budget_np, bias_2_np


