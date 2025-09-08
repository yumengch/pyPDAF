.. function:: py__init_obs_f_pdaf

    Provide the observation vector for the current time step.

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
