.. function:: py__likelihood_pdaf

    Compute the likelihood of the observation for a given ensemble member.

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
