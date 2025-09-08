.. function:: py__init_obserr_f_pdaf

    Initializes the full vector of observations error standard deviations.

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
