.. function:: py__init_obs_covar_pdaf

    Provide observation error covariance matrix to PDAF.

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
