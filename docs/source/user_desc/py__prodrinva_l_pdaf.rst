.. function:: py__prodrinva_l_pdaf

    Provide :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.

    Here, one should do :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l`.
    The matrix :math:`\mathbf{A}_l` depends on the filter algorithm.

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
        Observation vector in local analysis domain. shape: (dim_obs_l, )
    a_l: np.ndarray[np.float, dim=2]
        Matrix A in local analysis domain.
        shape: (dim_obs_l, rank)
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)

    Returns
    -------
    c_l: np.ndarray[np.float, dim=2]
        :math:`\mathbf{R}^{-1}_l \times \mathbf{A}_l` in local analysis domain.
        shape: (dim_obs_l, rank)
