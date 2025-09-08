.. function:: py__init_dim_obs_l_pdaf

    Initialise the dimension of local analysis domain observation vector.

    One can simplify this function by using :func:`pyPDAF.PDAFomi.init_dim_obs_l_xxx`.

    Parameters
    ----------
    domain_p: int
        Current local domain index
    step: int
        Current step
    dim_obs_f: int
        Global observation vector dimension.
    dim_obs_l: int
        Local observation vector dimension.

    Returns
    -------
    dim_obs_l: int
        Local observation vector dimension.
