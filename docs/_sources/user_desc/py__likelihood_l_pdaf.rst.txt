.. function:: py__likelihood_l_pdaf

    Compute the likelihood of the observation for a given ensemble member
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
