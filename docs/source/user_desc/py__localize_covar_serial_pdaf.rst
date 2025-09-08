.. function:: py__localize_covar_serial_pdaf

    Apply covariance localisation in EnSRF/EAKF.

    The localisation is applied to each observation element. The weight can be
    obtained by :func:`pyPDAF.PDAF.local_weight`.

    Parameters
    ----------
    iobs: int
        Index of the observation element.
    dim_p: int
        Dimension of the state vector.
    dim_obs: int
        Dimension of the observation vector.
    hp_p: np.ndarray[np.float, dim=1]
        Matrix HP. Shape: (dim_p)
    hph: np.ndarray[np.float, dim=1]
        Matrix HPH.T. Shape: (dim_obs)

    Returns
    -------
    hp_p: np.ndarray[np.float, dim=1]
        Localised matrix HP. Shape: (dim_p)
    hph: np.ndarray[np.float, dim=1]
        Localised matrix HPH.T. Shape: (dim_obs)
