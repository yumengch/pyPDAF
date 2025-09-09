py__init_obs_pdaf
=================

.. py:function:: py__init_obs_pdaf(step: int, dim_obs_p: int, observation_p: np.ndarray) -> np.ndarray

    Provide the observation vector for the current time step.

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