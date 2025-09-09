py__prodrinva_pdaf
==================

.. py:function:: py__prodrinva_pdaf(step: int, dim_obs_p: int, rank: int, obs_p: np.ndarray, a_p: np.ndarray, c_p: np.ndarray) -> np.ndarray

    Provide :math:`\mathbf{R}^{-1} \times \mathbf{A}`.

    Here, one should compute :math:`\mathbf{R}^{-1} \times \mathbf{A}` where
    :math:`\mathbf{R}` is observation error covariance matrix.
    The matrix :math:`\mathbf{A}` depends on the filter algorithm. In ESTKF,
    :math:`\mathbf{R}` can is ensemble perturbation in observation space.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of observation vector
    rank: int
        Rank of the matrix A (second dimension of A)
        This is ensemble size for ETKF and ensemble size - 1 for ESTKF.
    obs_p: np.ndarray[np.float, dim=1]
        Observation vector. shape: (dim_obs_p,)
    a_p: np.ndarray[np.float, dim=2]
        Input matrix A. shape: (dim_obs_p, rank)
    c_p: np.ndarray[np.float, dim=2]
        Output matrix :math:`\mathbf{C} = \mathbf{R}^{-1} \times \mathbf{A}`.
        shape: (dim_obs_p, rank)

    Returns
    -------
    c_p: np.ndarray[np.float64, dim=2]
        Output matrix :math:`\mathbf{C} = \mathbf{R}^{-1} \times \mathbf{A}`.
        shape: (dim_obs_p, rank)
