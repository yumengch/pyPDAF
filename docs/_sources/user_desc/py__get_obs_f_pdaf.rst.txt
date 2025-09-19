py__get_obs_f_pdaf
==================

.. py:function:: py__get_obs_f_pdaf(step: int, dim_obs_f: int, observation_f: np.ndarray) -> np.ndarray

    Receive synthetic observations from PDAF.

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
