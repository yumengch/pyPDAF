py__next_observation_pdaf
=========================

.. py:function:: py__next_observation_pdaf(stepnow: int, nsteps: int, doexit: int, time: float) -> Tuple[int, int, float]

    Get the number of time steps to be computed in the forecast phase.

    At the beginning of a forecast phase, this is called once by
        * :func:`pyPDAF.PDAF.init_forecast`
        * :func:`pyPDAF.PDAF3.assimilate_X`
        * ...

    Parameters
    ----------
    stepnow : int
            the current time step given by PDAF

    Returns
    -------
    nsteps : int
            number of forecast time steps until next assimilation;
            this can also be interpreted as
            number of assimilation function calls
            to perform a new assimilation
    doexit : int
            whether to exit forecasting (1 for exit)
    time : double
            current model (physical) time
