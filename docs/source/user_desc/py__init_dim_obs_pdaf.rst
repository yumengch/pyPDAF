py__init_dim_obs_pdaf
=====================

.. py:function:: py__init_dim_obs_pdaf(step: int, dim_obs_p: int) -> int

    Determine the size of the vector of observations

    The primary purpose of this function is to
    obtain the dimension of the observation vector.
    In OMI, in this function, one also sets the properties
    of `obs_f`, read the observation vector from
    files, setting the observation error variance
    when diagonal observation error covariance matrix
    is used. The `pyPDAF.PDAF.omi_gather_obs` function
    is also called here.

    Furthermore, in this user-supplied function, one also sets the interpolation
    coefficients used by observation operators.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.

    Returns
    -------
    dim_obs_p : int
        Dimension of the observation vector.
