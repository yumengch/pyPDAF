py__prodrinva_hyb_l_pdaf
========================

.. py:function:: py__prodrinva_hyb_l_pdaf(domain_p: int, step: int, dim_obs_l: int, rank: int, obs_l: np.ndarray, gamma: float, a_l: np.ndarray, c_l: np.ndarray) -> np.ndarray

    Provide :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` with weighting.

    Here, one should do :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.
    The matrix :math:`\mathbf{A}_l` depends on the filter algorithm.

    This function is used in LKNETF where `gamma` is multipled with `c_l` for
    weighting between LETKF and LNETF.

    One can also perform observation localisation. This can be helped by the function
    :func:`pyPDAF.PDAF.local_weight` to get the observation weight.

    Parameters
    ----------
    domain_p: int
        Current local domain index.
    step: int
        Current time step
    dim_obs_l: int
        Dimension of observation vector in local analysis domain
    rank: int
        Rank of the local analysis domain
        The size of it dpends on the filter algorithms.
    obs_l: np.ndarray[np.float, dim=1]
        Observation vector in local analysis domain.  shape: (dim_obs_l, )
    gamma: float
        Hybrid weight provided by PDAF
    a_l: np.ndarray[np.float, dim=2]
        Matrix A in local analysis domain. shape: (dim_obs_l, rank)
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)

    Returns
    -------
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)
