py__localize_covar_pdaf
========================

.. py:function:: py__localize_covar_pdaf(dim_p: int, dim_obs: int, hp_p: np.ndarray, hph: np.ndarray) -> Tuple[np.ndarray, np.ndarray]

    Perform covariance localisation.

    This is only used for stochastic EnKF. The localisation is performed
    for HP and HPH.T.

    This can be helped by function :func:`pyPDAF.PDAFomi.localize_covar`.
    This is replaced by :func:`PDAFomi.set_localize_covar` in PDAF3.

    Parameters
    ----------
    dim_p: int
        Dimension of the state vector.
    dim_obs: int
        Dimension of the observation vector.
    hp_p: np.ndarray[np.float, dim=2]
        Matrix HP. Shape: (dim_obs, dim_p)
    hph: np.ndarray[np.float, dim=2]
        Matrix HPH.T. Shape: (dim_obs, dim_obs)

    Returns
    -------
    hp_p: np.ndarray[np.float, dim=2]
        Localised matrix HP. Shape: (dim_obs, dim_p)
    hph: np.ndarray[np.float, dim=2]
        Localised matrix HPH.T. Shape: (dim_obs, dim_obs)
