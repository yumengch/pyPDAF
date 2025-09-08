# pylint: disable=unused-argument, too-many-lines
"""Stub file for user-supplied functions.
"""
from typing import Tuple
import numpy as np

def add_obs_err_pdaf(step: int, dim_obs_p: int, c_p: np.ndarray) -> np.ndarray:
    """Add the observation error covariance matrix to the matrix C.

    The input matrix is the projection of the ensemble
    covariance matrix onto the observation space that is computed
    during the analysis step of the stochastic EnKF. That is, HPH.T.
    The function returns HPH.T + R.

    The operation is for the global observation space.
    Thus, it is independent of whether the filter is executed with or
    without parallelization.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Size of observation vector
    C_p : ndarray[np.float64, ndim=2]
        Matrix to which the observation error covariance matrix is added
        shape: (dim_obs_p, dim_obs_p)

    Returns
    -------
    C_p : ndarray[np.float64, ndim=2]
        Matrix with added obs error covariance, i.e., HPH.T + R
        shape: (dim_obs_p, dim_obs_p)
    """

def init_ens_pdaf(filtertype: int, dim_p: int, dim_ens: int, state_p: np.ndarray,
    uinv: np.ndarray, ens_p: np.ndarray, flag: int
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Fill the ensemble array that is provided by PDAF with an initial ensemble of model states.

    This function is called by :func:`pyPDAF.PDAF.init`. The initialised
    ensemble array will be distributed to model by :func:`pyPDAF.PDAF.init_forecast`.

    Parameters
    ----------
    filtertype : int
            filter type given in PDAF_init
    dim_p : int
            PE-local state dimension given by PDAF_init
    dim_ens : int
            number of ensemble members
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag
    """

def next_observation_pdaf(stepnow: int, nsteps: int,
    doexit: int,
    time: float
) -> Tuple[int, int, float]:
    """Get the number of time steps to be computed in the forecast phase.

    At the beginning of a forecast phase, this is called once by
        * :func:`pyPDAF.PDAF.init_forecast`
        * :func:`pyPDAF.PDAF3.assimilate_X`
        * ...

    Parameters
    ----------
    stepnow : int
            the current time step given by PDAF

    Returns
    -------
    nsteps : int
            number of forecast time steps until next assimilation;
            this can also be interpreted as
            number of assimilation function calls
            to perform a new assimilation
    doexit : int
            whether to exit forecasting (1 for exit)
    time : double
            current model (physical) time
    """

def collect_state_pdaf(
    dim_p: int,
    state_p: np.ndarray
) -> np.ndarray:
    """Collect state vector from model/any arrays to pdaf arrays

    Parameters
    ----------
    dim_p : int
        pe-local state dimension

    state_p : ndarray[tuple[dim_p, ...], np.float64]
        local state vector

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        local state vector
    """

def distribute_state_pdaf(
    dim_p: int,
    state_p: np.ndarray
) -> np.ndarray:
    """Distribute a state vector from pdaf to the model/any arrays

    Parameters
    ----------
    dim_p : int
            PE-local state dimension
    state_p : ndarray[np.float64, ndim=1]
            PE-local state vector
            Array shape: (dim_p)

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            PE-local state vector
            Array shape: (dim_p)
    """

def prepoststep_pdaf(
    step: int,
    dim_p: int,
    dim_ens: int,
    dim_ens_l: int,
    dim_obs_p: int,
    state_p: np.ndarray,
    uinv: np.ndarray,
    ens_p: np.ndarray,
    flag: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Process ensemble before or after DA.

    Parameters
    ----------
    step : int
            current time step
            (negative for call before analysis/preprocessing)
    dim_p : int
            PE-local state vector dimension
    dim_ens : int
            number of ensemble members
    dim_ens_l : int
            number of ensemble members run serially
            on each model task
    dim_obs_p : int
            PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
            pe-local forecast/analysis state
            (the array 'state_p' is generally not
            initialised in the case of ESTKF/ETKF/EnKF/SEIK,
            so it can be used freely here.)
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            Inverse of the transformation matrix in ETKF and ESKTF;
            inverse of matrix formed by right singular vectors of error
            covariance matrix of ensemble perturbations in SEIK/SEEK.
            not used in EnKF.
            Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            pe-local forecast/analysis state
            (the array 'state_p' is generally not
            initialised in the case of ESTKF/ETKF/EnKF/SEIK,
            so it can be used freely here.)
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            Inverse of the transformation matrix in ETKF and ESKTF;
            inverse of matrix formed by right singular vectors of error
            covariance matrix of ensemble perturbations in SEIK/SEEK.
            not used in EnKF.
            Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    """

def init_dim_obs_pdaf(step: int, dim_obs_p: int) -> int:
    """Determine the size of the vector of observations

    The primary purpose of this function is to
    obtain the dimension of the observation vector.
    In OMI, in this function, one also sets the properties
    of `obs_f`, read the observation vector from
    files, setting the observation error variance
    when diagonal observation error covariance matrix
    is used. The `pyPDAF.PDAF.omi_gather_obs` function
    is also called here.

    Furthermore, in this user-supplied function, one also sets the interpolation
    coefficients used by observation operators.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.

    Returns
    -------
    dim_obs_p : int
        Dimension of the observation vector.
    """

def init_dim_obs_f_pdaf(step: int, dim_obs_f: int) -> int:
    """Determine the size of the full observations vector

    This function is used with domain localised filters to obtain the dimension of full
    observation vector that contains observations outside process-local domains.
    The smallest vector only needs observations within localisation radius of the
    process-local domain.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_f : int
        Dimension of the full observation vector.

    Returns
    -------
    dim_obs_f : int
        Dimension of the full observation vector.
    """

def init_obs_pdaf(step: int, dim_obs_p: int, observation_p: np.ndarray) -> np.ndarray:
    """Provide the observation vector for the current time step.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.
    observation_p : ndarray[np.float64, ndim=1]
        Observation vector.

    Returns
    -------
    observation_p : ndarray[np.float64, ndim=1]
        Filled observation vector.
    """

def init_obs_f_pdaf(step: int, dim_obs_f: int, observation_f: np.ndarray) -> np.ndarray:
    """Provide the observation vector for the current time step.

    This function is used with domain localised filters to obtain a full
    observation vector that contains observations outside process-local domains.
    The smallest vector only needs observations within localisation radius of the
    process-local domain.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_f : int
        Dimension of the observation vector.
    observation_f : ndarray[np.float64, ndim=1]
        Observation vector.

    Returns
    -------
    observation_f : ndarray[np.float64, ndim=1]
        Filled observation vector.
    """

def init_obs_covar_pdaf(step: int, dim_obs: int, dim_obs_p: int,
                        covar: np.ndarray, obs_p: np.ndarray,
                        isdiag: bool) -> Tuple[np.ndarray, bool]:
    """Provide observation error covariance matrix to PDAF.

    This function is used in stochastic EnKF for generating observation perturbations.

    Parameters
    ----------
    step: int
        current time step
    dim_obs : int
        dimension of global observation vector
    dim_obs_p: int
        dimension of process-local observation vector
    covar: np.ndarray[np.float64, dim=2]
        Observation error covariance matrix. shape: (dim_obs_p, dim_obs_p)
    obs_p: np.ndarray[np.float64, dim=1]
        Process-local observation vector. shape: dim_obs_p
    isdiag: bool
        Flag indicating if the covariance matrix is diagonal.

    Returns
    -------
    covar: np.ndarray[np.float64, dim=2]
        Observation error covariance matrix. shape: (dim_obs_p, dim_obs_p)
    is_diag: bool
        Flag indicating if the covariance matrix is diagonal.
    """

def init_obsvar_pdaf(step: int, dim_obs_p: int, obs_p: np.ndarray, meanvar: float) -> float:
    """Compute mean observation error variance.

    This is used by ETKF-variants for adaptive forgetting factor (type_forget=1).
    This can be global mean, or sub-domain mean.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_p: int
        Dimension of process-local observation vector
    obs_p: np.ndarray[np.float64, dim=1]
        Process-local observation vector. shape: (dim_obs_p,)
    meanvar: float
        Mean observation error variance.

    Returns
    -------
    meanvar: float
        Mean observation error variance.
    """

def init_obsvars_pdaf(step: int, dim_obs_f: int, var_f: np.ndarray) -> np.ndarray:
    """Provide a vector observation variance.

    This is used by EnSRF/EAKF.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_f: int
        Dimension of observation vector
    var_f: np.ndarray[np.float64, dim=1]
        Observation variance vector. shape: (dim_obs_f,)

    Returns
    -------
    var_f: np.ndarray[np.float64, dim=1]
        Observation variance vector. shape: (dim_obs_f,)
    """

def prodrinva_pdaf(step: int, dim_obs_p: int, rank: int,
                   obs_p: np.ndarray, a_p: np.ndarray, c_p: np.ndarray) -> np.ndarray:
    r"""Provide :math:`\mathbf{R}^{-1} \times \mathbf{A}`.

    Here, one should compute :math:`\mathbf{R}^{-1} \times \mathbf{A}` where
    :math:`\mathbf{R}` is observation error covariance matrix.
    The matrix :math:`\mathbf{A}` depends on the filter algorithm. In ESTKF,
    :math:`\mathbf{R}` can is ensemble perturbation in observation space.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of observation vector
    rank: int
        Rank of the matrix A (second dimension of A)
        This is ensemble size for ETKF and ensemble size - 1 for ESTKF.
    obs_p: np.ndarray[np.float, dim=1]
        Observation vector. shape: (dim_obs_p,)
    a_p: np.ndarray[np.float, dim=2]
        Input matrix A. shape: (dim_obs_p, rank)
    c_p: np.ndarray[np.float, dim=2]
        Output matrix :math:`\mathbf{C} = \mathbf{R}^{-1} \times \mathbf{A}`.
        shape: (dim_obs_p, rank)

    Returns
    -------
    c_p: np.ndarray[np.float64, dim=2]
        Output matrix :math:`\mathbf{C} = \mathbf{R}^{-1} \times \mathbf{A}`.
        shape: (dim_obs_p, rank)
    """

def obs_op_pdaf(step: int, dim_p: int, dim_obs_p: int,
                state_p: np.ndarray, m_state_p: np.ndarray) -> np.ndarray:
    r"""Apply observation operator

    This function computes :math:`\mathbf{H} \mathbf{x}`, where
    :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{x}` is state vector.

    Parameters
    ----------
    step : int
        Current step
    dim_p : int
        Dimension of state vector
    dim_obs_p: int
        Dimension of observation vector
    state_p: np.ndarray[np.float, dim=1]
        State vector. shape: (dim_p,)
    m_state_p: np.ndarray[np.float, dim=1]
        Observed state vector. shape: (dim_obs_p,)

    Returns
    -------
    m_state_p: np.ndarray[np.float64, dim=1]
        Observed state vector. shape: (dim_obs_p,)
    """

def obs_op_f_pdaf(step: int, dim_p: int, dim_obs_p: int,
                  state_p: np.ndarray, m_state_p: np.ndarray) -> np.ndarray:
    r"""Apply observation operator for full observed state vector

    This function computes :math:`\mathbf{H} \mathbf{x}`, where
    :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{x}` is state vector.
    See `here
    <https://pdaf.awi.de/trac/wiki/ImplementAnalysislestkf#Explanationoffullobservations>`_
    for the meaning of full observations.

    Parameters
    ----------
    step : int
        Current step
    dim_p : int
        Dimension of state vector
    dim_obs_p: int
        Dimension of observation vector
    state_p: np.ndarray[np.float, dim=1]
        State vector. shape: (dim_p,)
    m_state_p: np.ndarray[np.float, dim=1]
        Observed state vector. shape: (dim_obs_p,)

    Returns
    -------
    m_state_p: np.ndarray[np.float64, dim=1]
        Observed state vector. shape: (dim_obs_p,)
    """

def g2l_obs_pdaf(domain_p: int, step: int, dim_obs_f: int,
                 dim_obs_l: int, mstate_f: np.ndarray, mstate_l: np.ndarray) -> np.ndarray:
    """Convert global observed state vector to local vector.

    This is used by domain localisation methods. In these methods, each local
    domain has their own observation vector and observed state vector.

    Parameters
    ----------
    domain_p:int
        Current local domain index
    step: int
        Current time step
    dim_obs_f: int
        Global observation vector dimension.
    dim_obs_l: int
        Local observation vector dimension.
    mstate_f: np.ndarray[np.float64, dim=1]
        Global observed state vector. shape: (dim_obs_f,)
    mstate_l: np.ndarray[np.float64, dim=1]
        Local observed state vector. shape: (dim_obs_l,)

    Returns
    -------
    mstate_l: np.ndarray[np.float64, dim=1]
        Local observed state vector. shape: (dim_obs_l,)
    """

def g2l_state_pdaf(step: int, domain_p: int, dim_p: int,
                   state_p: np.ndarray, dim_l: int, state_l: np.ndarray) -> np.ndarray:
    """Get local state vector.

    Get the state vector for analysis local domain.

    Parameters
    ----------
    step: int
        Current time step
    domain_p: int
        Current local domain index
    dim_p: int
        Process-local state vector dimension.
    state_p: np.ndarray[np.float, dim=1]
        Process-local state vector.
    dim_l: int
        Local analysis domain state vector dimension.
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector. Shape: (dim_l,)

    Returns
    -------
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector. Shape: (dim_l,)
    """

def init_dim_l_pdaf(step: int, domain_p: int, dim_l: int) -> int:
    """Initialise local analysis domain state vector dimension.

    When PDAFlocal is used, one should call :func:`pyPDAF.PDAFlocal.set_indices` here.

    Parameters
    ----------
    step: int
        Current step
    domain_p: int
        Current local domain index
    dim_l: int
        Local analysis domain state vector dimension.

    Returns
    -------
    dim_l: int
        Local analysis domain state vector dimension.
    """

def init_dim_obs_l_pdaf(domain_p: int, step: int, dim_obs_f: int, dim_obs_l: int) -> int:
    """Initialise the dimension of local analysis domain observation vector.

    One can simplify this function by using :func:`pyPDAF.PDAFomi.init_dim_obs_l_xxx`.

    Parameters
    ----------
    domain_p: int
        Current local domain index
    step: int
        Current step
    dim_obs_f: int
        Global observation vector dimension.
    dim_obs_l: int
        Local observation vector dimension.

    Returns
    -------
    dim_obs_l: int
        Local observation vector dimension.
    """

def init_n_domains_p_pdaf(step: int, n_domains_p: int) -> int:
    """Get number of analysis domains.

    Parameters
    ----------
    step: int
        Current step
    n_domains_p: int
        Number of analysis domains.

    Returns
    -------
    n_domains_p: int
        Number of analysis domains.
    """

def init_obs_l_pdaf(domain_p: int, step: int, dim_obs_l: int,
                    observation_l: np.ndarray) -> np.ndarray:
    """Initialise observation vector for local analysis domain.

    Parameters
    ----------
    domain_p: int
        Current local domain index
    step: int
        Current time step
    dim_obs_l: int
        Local observation vector dimension.
    observation_l: np.ndarray[np.float, dim=1]
        Local observation vector.

    Returns
    -------
    observation_l: np.ndarray[np.float, dim=1]
        Local observation vector.
    """

def init_obsvar_l_pdaf(domain_p: int, step: int, dim_obs_l: int,
                       obs_l: np.ndarray, dim_obs_p: int, meanvar_l: float) -> float:
    """Get mean of analysis domain local observation variance.

    This is used by local adaptive forgetting factor (type_forget=2)

    Parameters
    ----------
    domain_p: int
        Current local domain index
    step: int
        Current time step
    dim_obs_l: int
        Local observation vector dimension.
    obs_l: np.ndarray[np.float, dim=1]
        Local observation vector.
    dim_obs_p: int
        Process-local observation vector dimension.
    meanvar_l: float
        Mean of analysis domain local observation variance.

    Returns
    -------
    meanvar_l: float
        Mean of analysis domain local observation variance.
    """

def init_obserr_f_pdaf(step: int, dim_obs_f: int,
                       obs_f: np.ndarray, obserr_f: np.ndarray) -> np.ndarray:
    """Initializes the full vector of observations error standard deviations.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_f : int
        Full observation vector dimension.
    obs_f : np.ndarray[np.float, dim=1]
        Full observation vector.
    obserr_f : np.ndarray[np.float, dim=1]
        Full observation error vector.

    Returns
    -------
    obserr_f : np.ndarray[np.float, dim=1]
        Full observation error vector.
    """

def l2g_state_pdaf(step: int, domain_p: int, dim_l: int,
                   state_l: np.ndarray, dim_p: int, state_p: np.ndarray) -> np.ndarray:
    """Assign local state vector to process-local global state vector.

    Parameters
    ----------
    step: int
        Current time step
    domain_p: int
        Current local domain index
    dim_l: int
        Local analysis domain state vector dimension.
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector.
    dim_p: int
        Process-local global state vector dimension.
    state_p: np.ndarray[np.float, dim=1]
        Process-local global state vector.

    Returns
    -------
    state_p: np.ndarray[np.float, dim=1]
        Process-local global state vector.
    """

def prodrinva_l_pdaf(domain_p: int, step: int, dim_obs_l: int,
                     rank: int, obs_l: np.ndarray, a_l: np.ndarray,
                     c_l: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    r"""Provide :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.

    Here, one should do :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.
    The matrix :math:`\mathbf{A}_l` depends on the filter algorithm.

    One can also perform observation localisation. This can be helped by the function
    :func:`pyPDAF.PDAF.local_weight` to get the observation weight.

    Parameters
    ----------
    domain_p: int
        Current local domain index.
    step: int
        Current time step
    dim_obs_l: int
        Dimension of observation vector in local analysis domain
    rank: int
        Rank of the local analysis domain
        The size of it dpends on the filter algorithms.
    obs_l: np.ndarray[np.float, dim=1]
        Observation vector in local analysis domain. shape: (dim_obs_l, )
    a_l: np.ndarray[np.float, dim=2]
        Matrix A in local analysis domain.
        shape: (dim_obs_l, rank)
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)

    Returns
    -------
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)
    """

def localize_covar_pdaf(dim_p: int, dim_obs: int, hp_p: np.ndarray,
                        hph: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Perform covariance localisation.

    This is only used for stochastic EnKF. The localisation is performed
    for HP and HPH.T.

    This can be helped by function :func:`pyPDAF.PDAFomi.localize_covar`.
    This is replaced by :func:`PDAFomi.set_localize_covar` in PDAF3.

    Parameters
    ----------
    dim_p: int
        Dimension of the state vector.
    dim_obs: int
        Dimension of the observation vector.
    hp_p: np.ndarray[np.float, dim=2]
        Matrix HP. Shape: (dim_obs, dim_p)
    hph: np.ndarray[np.float, dim=2]
        Matrix HPH.T. Shape: (dim_obs, dim_obs)

    Returns
    -------
    hp_p: np.ndarray[np.float, dim=2]
        Localised matrix HP. Shape: (dim_obs, dim_p)
    hph: np.ndarray[np.float, dim=2]
        Localised matrix HPH.T. Shape: (dim_obs, dim_obs)
    """

def localize_covar_serial_pdaf(iobs: int, dim_p: int, dim_obs:int,
                               hp_p: np.ndarray, hph: np.ndarray,
                               hxy_p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Apply covariance localisation in EnSRF/EAKF.

    The localisation is applied to each observation element. The weight can be
    obtained by :func:`pyPDAF.PDAF.local_weight`.

    Parameters
    ----------
    iobs: int
        Index of the observation element.
    dim_p: int
        Dimension of the state vector.
    dim_obs: int
        Dimension of the observation vector.
    hp_p: np.ndarray[np.float, dim=1]
        Matrix HP. Shape: (dim_p)
    hph: np.ndarray[np.float, dim=1]
        Matrix HPH.T. Shape: (dim_obs)

    Returns
    -------
    hp_p: np.ndarray[np.float, dim=1]
        Localised matrix HP. Shape: (dim_p)
    hph: np.ndarray[np.float, dim=1]
        Localised matrix HPH.T. Shape: (dim_obs)
    """

def likelihood_pdaf(step: int, dim_obs_p: int, obs_p: np.ndarray,
                    resid: np.ndarray, likely: float) -> float:
    r"""Compute the likelihood of the observation for a given ensemble member.

    The function is used with the nonlinear filter NETF and particle filter. The likelihood
    depends on the assumed observation error distribution.
    For a Gaussian observation error, the likelihood is
    :math:`\exp(-0.5(\mathbf{y}-\mathbf{H}\mathbf{x})^\mathrm{T}R^{-1}(\mathbf{y}-\mathbf{H}\mathbf{x}))`.
    The vector :math:`\mathbf{y}-\mathbf{H}\mathbf{x} = \mathrm{resid}` is
    provided as an input argument.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_p : int
        Dimension of the observation vector.
    obs_p: np.ndarray[np.float, dim=1]
        Observation vector. Shape: (dim_obs_p)
    resid: np.ndarray[np.float, dim=1]
        Residual vector between observations and state. Shape: (dim_obs_p)
    likely: float
        Likelihood of the observation

    Returns
    -------
    likely: float
        Likelihood of the observation
    """

def likelihood_l_pdaf(domain_p: int, step: int, dim_obs_l: int,
                      obs_l: np.ndarray, resid_l: np.ndarray,
                      double: float) -> float:
    r"""Compute the likelihood of the observation for a given ensemble member
    according to the observations used for the local analysis.

    The function is used in the localized nonlinear filter LNETF. The likelihood
    depends on the assumed observation error distribution.
    For a Gaussian observation error, the likelihood is
    :math:`\exp(-0.5(\mathbf{y}-\mathbf{H}\mathbf{x})^\mathrm{T}R^{-1}(\mathbf{y}-\mathbf{H}\mathbf{x}))`.
    The vector :math:`\mathbf{y}-\mathbf{H}\mathbf{x} = \mathrm{resid}` is
    provided as an input argument.

    This function is also the place to perform observation localisation.
    To initialize a vector of weights, the routine :func:`pyPDAF.PDAF.local_weight`
    can be called.

    Parameters
    ----------
    domain_p: int
        Current local analysis domain index
    step: int
        Current time step
    dim_obs_l: int
        Dimension of the local observation vector.
    obs_l: np.ndarray[np.float, dim=1]
        Observation vector. Shape: (dim_obs_l)
    resid_l: np.ndarray[np.float, dim=1]
        Residual vector between observations and state. Shape: (dim_obs_l)
    likely_l: float
        Likelihood of the local observation

    Returns
    -------
    likely_l: float
        Likelihood of the local observation
    """

def get_obs_f_pdaf(step: int, dim_obs_f: int, observation_f: np.ndarray) -> np.ndarray:
    """Receive synthetic observations from PDAF.

    This function is used in twin experiments for observation generations.
    One can, for example, save synthetic observations in this function.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_f: int
        Full observation vector dimension.
    observation_f: np.ndarray[np.float, dim=1]
        Full observation vector.

    Returns
    -------
    observation_f: np.ndarray[np.float, dim=1]
        Full observation vector.
    """

def cvt_adj_ens_pdaf(iter: int, dim_p: int, dim_ens: int, dim_cv_ens_p:int,
                     ens: np.ndarray, vcv_p: np.ndarray, cv_p: np.ndarray) -> np.ndarray:
    r"""The adjoint control variable transformation involving ensembles.

    Here, this function performs
    :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. Here, the input vector can be
    :math:`\mathbf{v}_h = \mathbf{H}^\mathrm{T}\mathbf{R}^{-1}\delta \mathbf{x}`
    with :math:`\delta \mathbf{x}` the innovation vector.

    The control vector transform is given by :math:`\mathbf{U}\mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DEnVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_ens: int
        Dimension of the ensemble.
    dim_cv_ens_p: int
        Dimension of the control vector.
    ens_p: np.ndarray[np.float, dim=2]
        Ensemble matrix. shape: (dim_p, dim_ens)
    vcv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}_h`. shape: (dim_p, )
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cv_ens_p, )

    Returns
    -------
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cv_ens_p, )
    """

def cvt_adj_pdaf(iter: int, dim_p: int, dim_cvec: int,
                 vcv_p: np.ndarray, cv_p: np.ndarray) -> np.ndarray:
    r"""The adjoint control variable transformation.

    Here, this function performs
    :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. Here, the input vector can be
    :math:`\mathbf{v}_h = \mathbf{H}^\mathrm{T}\mathbf{R}^{-1}\delta \mathbf{x}`
    with :math:`\delta \mathbf{x}` the innovation vector.

    The control vector transform is given by :math:`\mathbf{U}\mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_cvec: int
        Dimension of the control vector.
    vcv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}_h`. shape: (dim_p, )
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cvec, )

    Returns
    -------
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cvec, )
    """

def cvt_pdaf(iter: int, dim_p: int, dim_cvec: int,
             cv_p: np.ndarray, vv_p: np.ndarray) -> np.ndarray:
    r"""The control variable transformation.

    Here, this function performs :math:`\mathbf{U} \mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_cvec: int
        Dimension of the control vector.
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}`. shape: (dim_cvec, )
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )

    Returns
    -------
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )
    """

def cvt_ens_pdaf(iter: int, dim_p: int, dim_ens: int, dim_cv_ens_p: int,
                 v_p: np.ndarray, vv_p: np.ndarray) -> np.ndarray:
    r"""The control variable transformation involving ensembles.

    Here, this function performs :math:`\mathbf{U} \mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DEnVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_ens: int
        Dimension of the ensemble.
    dim_cv_ens_p: int
        Dimension of the control vector.
    ens_p: np.ndarray[np.float, dim=2]
        Ensemble matrix. shape: (dim_p, dim_ens)
    v_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}`. shape: (dim_cv_ens_p, )
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )

    Returns
    -------
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )
    """

def obs_op_adj_pdaf(step: int, dim_p: int, dim_obs_p: int,
                    m_state_p: np.ndarray, state_p: np.ndarray) -> np.ndarray:
    r"""Apply adjoint observation operator

    This function computes :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`, where
    :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{w}` is a vector in observation space.

    Parameters
    ----------
    step : int
        Current step
    dim_p : int
        Dimension of state vector
    dim_obs_p: int
        Dimension of observation vector
    m_state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{w}`. shape: (dim_obs_p,)
    state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`. shape: (dim_p,)

    Returns
    -------
    state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`. shape: (dim_p,)
    """

def obs_op_lin_pdaf(step: int, dim_p: int, dim_obs_p: int,
                    state_p: np.ndarray, m_state_p: np.ndarray) -> np.ndarray:
    r"""Apply linearised observation operator

    This function computes :math:`\mathbf{H} \mathbf{x}`, where
    :math:`\mathbf{H}` is the linearised observation operator and
    :math:`\mathbf{x}` is state vector.

    Parameters
    ----------
    step : int
        Current step
    dim_p : int
        Dimension of state vector
    dim_obs_p: int
        Dimension of observation vector
    state_p: np.ndarray[np.float, dim=1]
        State vector. shape: (dim_p,)
    m_state_p: np.ndarray[np.float, dim=1]
        Observed state vector. shape: (dim_obs_p,)

    Returns
    -------
    m_state_p: np.ndarray[np.float64, dim=1]
        Observed state vector. shape: (dim_obs_p,)
    """

def likelihood_hyb_l_pdaf(domain_p: int, step: int, dim_obs_l: int,
                          obs_l: np.ndarray, gamma: float, resid_l: np.ndarray,
                          likely_l: np.ndarray) -> float:
    r"""Compute the likelihood of the observation for a given ensemble member
    according to the observations used for the local analysis with hybrid weight.

    The function is used in the localized nonlinear filter LKNETF. The likelihood
    depends on the assumed observation error distribution.
    For a Gaussian observation error, the likelihood is
    :math:`\exp(-0.5(\mathbf{y}-\mathbf{H}\mathbf{x})^\mathrm{T}R^{-1}(\mathbf{y}-\mathbf{H}\mathbf{x}))`.
    The vector :math:`\mathbf{y}-\mathbf{H}\mathbf{x} = \mathrm{resid}` is
    provided as an input argument.

    The hybrid weight `gamma` is used weight between LNETF and LETKF. which is applied
    to :math:`R^{-1}(\mathbf{y}-\mathbf{H}\mathbf{x})`.

    This function is also the place to perform observation localisation.
    To initialize a vector of weights, the routine :func:`pyPDAF.PDAF.local_weight`
    can be called.

    Parameters
    ----------
    domain_p: int
        Current local analysis domain index
    step: int
        Current time step
    dim_obs_l: int
        Dimension of the local observation vector.
    obs_l: np.ndarray[np.float, dim=1]
        Observation vector. Shape: (dim_obs_l)
    gamma: float
        Hybrid weight provided by PDAF
    resid_l: np.ndarray[np.float, dim=1]
        Residual vector between observations and state. Shape: (dim_obs_l)
    likely_l: float
        Likelihood of the local observation

    Returns
    -------
    likely_l: float
        Likelihood of the local observation

    """

def prodrinva_hyb_l_pdaf(domain_p: int, step: int, dim_obs_l: int, rank: int,
                         obs_l: np.ndarray, gamma: float,
                         a_l: np.ndarray, c_l: np.ndarray) -> np.ndarray:
    r"""Provide :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` with weighting.

    Here, one should do :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.
    The matrix :math:`\mathbf{A}_l` depends on the filter algorithm.

    This function is used in LKNETF where `gamma` is multipled with `c_l` for
    weighting between LETKF and LNETF.

    One can also perform observation localisation. This can be helped by the function
    :func:`pyPDAF.PDAF.local_weight` to get the observation weight.

    Parameters
    ----------
    domain_p: int
        Current local domain index.
    step: int
        Current time step
    dim_obs_l: int
        Dimension of observation vector in local analysis domain
    rank: int
        Rank of the local analysis domain
        The size of it dpends on the filter algorithms.
    obs_l: np.ndarray[np.float, dim=1]
        Observation vector in local analysis domain.  shape: (dim_obs_l, )
    gamma: float
        Hybrid weight provided by PDAF
    a_l: np.ndarray[np.float, dim=2]
        Matrix A in local analysis domain. shape: (dim_obs_l, rank)
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)

    Returns
    -------
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)
    """
