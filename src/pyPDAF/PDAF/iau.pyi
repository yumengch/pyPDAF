import numpy as np
from typing import Callable, Tuple

def iau_init(type_iau_in: int, nsteps_iau_in: int) -> int:
    """Initialise parameters for incremental analysis updates, IAU.

    It further allocates the array in which the ensemble increments are stored.
    This array exists on all processes that are part of model tasks.

    This is usually called after :func:`pyPDAF.PDAF.init`.
    It is important that this is called by all model processes because all these
    processes need the information on the IAU configuration and
    need to allocate the increment array.

    Parameters
    ----------
    type_iau_in : int
        - (0) no IAU
        - (1) constant increment weight 1/nsteps_iau
        - (2) Linear IAU weight with maximum in middle of IAU period. This weight
              linearly increases and decreases. This type is usually used if
              the IAU is applied in re-running over the previous observation period.
        - (3) Zero weights for null mode (can be used to apply IAU on user side).
              This stores the increment information, but does not apply the increment.
              One can use :func:`pyPDAF.PDAF.iau_set_pointer` to access the increment
              array.
    nsteps_iau_in : int
        number of time steps in IAU

    Returns
    -------
    flag : int
        Status flag
    """

def iau_reset(type_iau_in: int, nsteps_iau_in: int) -> int:
    """Modify the IAU type and the number of IAU time steps during a run

    While :func:`pyPDAF.PDAF.iau_init` sets the IAU type and number of IAU steps
    initially, one can change these settings during a model run.

    This function has to be called by all processes that are model processes.
    A common place is to call the function after an analysis step or
    in `distribute_state_pdaf`.

    Parameters
    ----------
    type_iau_in : int
        - (0) no IAU
        - (1) constant increment weight 1/nsteps_iau
        - (2) Linear IAU weight with maximum in middle of IAU period. This weight
              linearly increases and decreases. This type is usually used if
              the IAU is applied in re-running over the previous observation period.
        - (3) Zero weights for null mode (can be used to apply IAU on user side).
              This stores the increment information, but does not apply the increment.
              One can use :func:`pyPDAF.PDAF.iau_set_pointer` to access the increment
              array.
    nsteps_iau_in : int
        number of time steps in IAU

    Returns
    -------
    flag : int
        Status flag
    """

def iau_set_weights(iweights: int, weights: np.ndarray) -> None:
    """Provide a user-specified vector of increment weights.

    While :func:`pyPDAF.PDAF.iau_init` allows to choose among
    pre-defined weight functions, one might like to use a different function
    and the corresponding weights can be set here.

    All model processes must call the routine.
    A common place is to call the function after an analysis step or in `distribute_state_pdaf`.

    Parameters
    ----------
    iweights : int
        Length of weights input vector
        If iweights is different from the number of IAU steps set in
        :func:`pyPDAF.PDAF.iau_init` or :func:`pyPDAF.PDAF.iau_reset`,
        only the minimum of iweights and the set IAU steps is filled with
        the provided weights vector.
    weights : ndarray[np.float64, ndim=1]
        Input weight vector
        Array shape: (iweights)
    """

def iau_set_pointer() -> Tuple[np.ndarray, int]:
    """Set a pointer to the ensemble increments array.

    This gives direct access to the increment array,
    e.g. to analyze it or to write it into a file for restarting.

    If it is called by each single process, but it only provides a pointer to
    the process-local part of the increment array.

    For domain-decomposed models, this array only includes the state vector
    part for the process domain. In addition, it usually only contains a
    sub-ensemble unless one uses the flexible parallelization mode with a
    single model task. For the fully parallel mode, the process(es) of a
    single model task only hold a single ensemble state.

    Returns
    -------
    iau_ptr_np : np.ndarray
        The increment array (process-local part)
    flag : int
        Status flag
    """

def iau_init_inc(dim_p: int, dim_ens_l: int, ens_inc: np.ndarray) -> int:
    """Fill the process-local increment array.

    A common use is when the IAU should be applied from the initial time of a run,
    for example if the increment was computed in a previous assimilation run
    and then stored for restarting. Since PDAF can only compute an increment
    in an analysis step, the user needs to provide the increment.

    The function is called after :func:`pyPDAF.PDAF.iau_init`.
    It has to be called by all processes that are model processes and one needs
    to provide the task-local ensemble
    (i.e. with local ensemble size `dim_ens_l=1` for the fully parallel mode,
    and usually `dim_ens_l>1` for the flexible parallelization mode).
    The function cannot be called in `init_ens_pdaf` since this function is only
    executed by filter processes instead of all model processes.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens_l : int
        Number of ensemble members that are run by a loop for each model task
    ens_inc : ndarray[np.float64, ndim=2]
        PE-local increment ensemble
        Array shape: (dim_p, dim_ens_l)

    Returns
    -------
    flag : int
        Status flag
    """

def iau_add_inc(py__collect_state_pdaf: Callable, py__distribute_state_pdaf: Callable) -> None:
    """Apply IAU to model forecasts in flexible parallel mode.

    PDAF automatically apply IAU in fully parallel model. However, one has to
    apply IAU with this function in flexible parallel mode.

    This has to be used in each model time stepping.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Collect a state vector for PDAF

    py__distribute_state_pdaf : Callable
        Distribute a state vector for PDAF
    """
