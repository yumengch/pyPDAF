.. function:: py__init_dim_l_pdaf

    Initialise local analysis domain state vector dimension.

    When PDAFlocal is used, one should call :func:`pyPDAF.PDAFlocal.set_indices` here.

    Parameters
    ----------
    step: int
        Current step
    domain_p: int
        Current local domain index
    dim_l: int
        Local analysis domain state vector dimension.

    Returns
    -------
    dim_l: int
        Local analysis domain state vector dimension.
