py__init_dim_obs_f_pdaf
=======================

.. py:function:: py__init_dim_obs_f_pdaf(step: int, dim_obs_f: int) -> int

    Determine the size of the full observations vector

    This function is used with domain localised filters to obtain the dimension of full
    observation vector that contains observations outside process-local domains.
    The smallest vector only needs observations within localisation radius of the
    process-local domain.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_f : int
        Dimension of the full observation vector.

    Returns
    -------
    dim_obs_f : int
        Dimension of the full observation vector.
