.. function:: py__init_obs_l_pdaf

    Initialise observation vector for local analysis domain.

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
