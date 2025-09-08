.. function:: py__init_obsvar_pdaf

    Compute mean observation error variance.

    This is used by ETKF-variants for adaptive forgetting factor (type_forget=1).
    This can be global mean, or sub-domain mean.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_p: int
        Dimension of process-local observation vector
    obs_p: np.ndarray[np.float64, dim=1]
        Process-local observation vector. shape: (dim_obs_p,)
    meanvar: float
        Mean observation error variance.

    Returns
    -------
    meanvar: float
        Mean observation error variance.
