py__obs_op_adj_pdaf
====================

.. py:function:: py__obs_op_adj_pdaf(step: int, dim_p: int, dim_obs_p: int, m_state_p: np.ndarray, state_p: np.ndarray) -> np.ndarray

    Apply adjoint observation operator

    This function computes :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`, where
    :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{w}` is a vector in observation space.

    Parameters
    ----------
    step : int
        Current step
    dim_p : int
        Dimension of state vector
    dim_obs_p: int
        Dimension of observation vector
    m_state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{w}`. shape: (dim_obs_p,)
    state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`. shape: (dim_p,)

    Returns
    -------
    state_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{H}^\mathrm{T} \mathbf{w}`. shape: (dim_p,)
