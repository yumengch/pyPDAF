py__obs_op_f_pdaf
===================

.. py:function:: py__obs_op_f_pdaf(step: int, dim_p: int, dim_obs_p: int, state_p: np.ndarray, m_state_p: np.ndarray) -> np.ndarray

    Apply observation operator for full observed state vector

    This function computes :math:`\mathbf{H} \mathbf{x}`, where
    :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{x}` is state vector.
    See `here <https://pdaf.awi.de/trac/wiki/ImplementAnalysislestkf#Explanationoffullobservations>`_
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
