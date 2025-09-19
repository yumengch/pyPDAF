py__init_obsvar_l_pdaf
=======================

.. py:function:: py__init_obsvar_l_pdaf(domain_p: int, step: int, dim_obs_l: int, obs_l: np.ndarray, dim_obs_p: int, meanvar_l: float) -> float

    Get mean of analysis domain local observation variance.

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
