py__init_obsvars_pdaf
=======================

.. py:function:: py__init_obsvars_pdaf(step: int, dim_obs_f: int, var_f: np.ndarray) -> np.ndarray

    Provide a vector observation variance.

    This is used by EnSRF/EAKF.

    Parameters
    ----------
    step: int
        Current time step
    dim_obs_f: int
        Dimension of observation vector
    var_f: np.ndarray[np.float64, dim=1]
        Observation variance vector. shape: (dim_obs_f,)

    Returns
    -------
    var_f: np.ndarray[np.float64, dim=1]
        Observation variance vector. shape: (dim_obs_f,)
