# pylint: disable=unused-argument
from typing import Tuple, Callable, Any
import numpy as np
from numpy.typing import NDArray


def correlation_function(ctype: int, length: float, distance: float) -> float:
    """The value of the chosen correlation function according to the specified
    length scale.

    Parameters
    ----------
    ctype : int
        Type of correlation function
        * ctype=1: Gaussian with f(0)=1.0
        * ctype=2: 5th-order polynomial (Gaspace/Cohn, 1999)
    length : float
        Length scale of function
        * ctype=1: standard deviation
        * ctype=2: support length f=0 for distance>length
    distance : float
        Distance at which the function is evaluated

    Returns
    -------
    value : float
        Value of the function
    """


def deallocate() -> None:
    """Finalise the PDAF systems including freeing some of the memory used by PDAF.

    This function cannot be used to free all allocated PDAF memory.
    Therefore, one should not use :func:`pyPDAF.PDAF.init` afterwards.
    """


def eofcovar(
    dim: int,
    nstates: int,
    nfields: int,
    dim_fields: NDArray[np.int32],
    offsets: NDArray[np.int32],
    remove_mstate: int,
    do_mv: int,
    states: NDArray[np.float64],
    meanstate: NDArray[np.float64],
    verbose: int,
) -> Tuple[
    NDArray[np.float64],  # states
    NDArray[np.float64],  # stddev
    NDArray[np.float64],  # svals
    NDArray[np.float64],  # svec
    NDArray[np.float64],  # meanstate
    int,                  # status
]:
    """
    EOF analysis of an ensemble of state vectors by singular value decomposition.

    Typically, this function is used with :func:`pyPDAF.PDAF.SampleEns`
    to generate an ensemble of a chosen size (up to the number of EOFs plus one).

    Here, the function performs a singular value decomposition
    of the ensemble anomaly of the input matrix,
    which is usually an ensemble formed by state vectors
    at multiple time steps.
    The singular values and corresponding singular vectors
    can be used to construct a covariance matrix.
    This can be used as the initial error covariance for the initial ensemble.

    A multivariate scaling can be performed to ensure that all fields
    in the state vectors have unit variance.

    It can be useful to store more EOFs than one finally
    might want to use to have the flexibility
    to carry the ensemble size.


    See Also
    --------
    `PDAF webpage <https://pdaf.awi.de/trac/wiki/EnsembleGeneration>`_

    Parameters
    ----------
    dim : int
        Dimension of state vector
    nstates : int
        Number of state vectors
    nfields : int
        Number of fields in state vector
    dim_fields : ndarray[np.intc, ndim=1]
        Size of each field
        Array shape: (nfields)
    offsets : ndarray[np.intc, ndim=1]
        Start position of each field
        Array shape: (nfields)
    remove_mstate : int
        1: subtract mean state from states
    do_mv : int
        1: Do multivariate scaling; 0: no scaling
    states : ndarray[np.float64, ndim=2]
        State perturbations
        Array shape: (dim, nstates)
    meanstate : ndarray[np.float64, ndim=1]
        Mean state (only changed if remove_mstate=1)
        Array shape: (dim)
    verbose : int
        Verbosity flag

    Returns
    -------
    states : ndarray[np.float64, ndim=2]
        State perturbations
        Array shape: (dim, nstates)
    stddev : ndarray[np.float64, ndim=1]
        Standard deviation of field variability
        Array shape: (nfields)
    svals : ndarray[np.float64, ndim=1]
        Singular values divided by sqrt(nstates-1)
        Array shape: (nstates)
    svec : ndarray[np.float64, ndim=2]
        Singular vectors
        Array shape: (dim, nstates)
    meanstate : ndarray[np.float64, ndim=1]
        Mean state (only changed if remove_mstate=1)
        Array shape: (dim)
    status : int
        Status flag
    """


def force_analysis() -> None:
    """Perform assimilation after this function call.

    This function overwrite member index of the ensemble state
    by local_dim_ens (number of ensembles for current process,
    in full parallel setup, this is 1.) and the counter
    cnt_steps by nsteps-1.
    This forces that the analysis step is executed at
    the next call to PDAF assimilation functions.
    """


def gather_dim_obs_f(dim_obs_p: int) -> int:
    """Gather the dimension of observation vector
    across multiple local domains/filter processors.

    This function is typically used in deprecated PDAF functions
    without OMI.

    This function can be used in the user-supplied function of
    :func:`py__init_dim_obs_f_pdaf`,
    but it is recommended to use :func:`pyPDAF.PDAF.omi_gather_obs`
    with OMI.

    This function does two things:
        1. Receiving observation dimension on each local process.
        2. Gather the total dimension of observation
           across local process and the displacement of PE-local
           observations relative to the total observation vector

    The dimension of observations are used to allocate observation
    arrays. Therefore, it must be used before
    :func:`pyPDAF.PDAF.gather_obs_f` or :func:`pyPDAF.PDAF.gather_obs_f2`.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension

    Returns
    -------
    dim_obs_f : int
        Full observation dimension
    """


def gather_obs_f(
    obs_p: NDArray[np.float64],
    dimobs_f: int,
) -> Tuple[NDArray[np.float64], int]:
    """In the local filters (LESKTF, LETKF, LSEIK, LNETF) this function
    returns the total observation vector from process-local observations.
    The function depends on :func:`pyPDAF.PDAF.gather_dim_obs_f` which defines
    the process-local observation dimensions. Further,
    the related routine :func:`pyPDAF.PDAF.gather_obs_f2` is used to
    gather the associated 2D observation coordinates

    Parameters
    ----------
    obs_p : ndarray[np.float64, ndim=1]
        PE-local vector
        Array shape: (dimobs_p)
    dimobs_f : int
        Full observation dimension, dimension of `obs_f`

    Returns
    -------
    obs_f : ndarray[np.float64, ndim=1]
        Full gathered vector
        Array shape: (dimobs_f)
    status : int
        Status flag:
    """


def gather_obs_f2(
    coords_p: NDArray[np.float64],
    nrows: int,
    dimobs_f: int,
) -> Tuple[NDArray[np.float64], int]:
    """In the local filters (LESKTF, LETKF, LSEIK, LNETF)
    this function returns the full observation coordinates from process-local
    observation coordinates. The function depends on
    func:`pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local
    observation dimensions. Further, the related routine
    func:`pyPDAF.PDAF.gather_obs_f` is used to gather the associated
    observation vectors.

    The routine is typically used in the routines `py__init_dim_obs_f_pdaf`
    if the analysis step of the local filters is parallelized.

    Parameters
    ----------
    coords_p : ndarray[np.float64, ndim=2]
        PE-local array
        Array shape: (nrows, dimobs_p)
    nrows : int
        Number of rows in array
    dimobs_f : int
        Full observation dimension, dimension of `coords_f`

    Returns
    -------
    coords_f : ndarray[np.float64, ndim=2]
        Full gathered array
        Array shape: (nrows, dimobs_f)
    status : int
        Status flag:
    """


def gather_obs_f_flex(
    dim_obs_p: int,
    dim_obs_f: int,
    obs_p: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], int]:
    """Gather full observation from processor
    local observation without PDAF-internal info.

    In the local filters (LESKTF, LETKF, LSEIK, LNETF)
    this function returns the full observation
    from process-local observation.

    This function has a similar functionality as
    :func:`pyPDAF.PDAF.gather_obs_f`, but it does not depend
    on PDAF-internal observation dimension information,
    which is specified once :func:`pyPDAF.PDAF.gather_dim_obs_f`
    is executed.

    This function can also be used as a generic function
    to gather any 2D arrays
    where the second dimension should be concatenated.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension
    dim_obs_f : int
        Full observation dimension
    obs_p : ndarray[np.float64, ndim=1]
        PE-local vector
        Array shape: (dim_obs_p)

    Returns
    -------
    obs_f : ndarray[np.float64, ndim=1]
        Full gathered vector
        Array shape: (dim_obs_f)
    status : int
        Status flag: (0) no error
    """


def gather_obs_f2_flex(
    dim_obs_p: int,
    dim_obs_f: int,
    coords_p: NDArray[np.float64],
    nrows: int,
) -> Tuple[NDArray[np.float64], int]:
    """Gather full observation coordinates from processor
    local observation coordinates without PDAF-internal info.

    In the local filters (LESKTF, LETKF, LSEIK, LNETF)
    this function returns the full observation coordinates
    from process-local observation coordinates.

    This function has a similar functionality as
    :func:`pyPDAF.PDAF.gather_obs_f2`, but it does not depend
    on PDAF-internal observation dimension information,
    which is specified once :func:`pyPDAF.PDAF.gather_dim_obs_f`
    is executed.

    This function can also be used as a generic function
    to gather any 2D arrays
    where the second dimension should be concatenated.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension
    dim_obs_f : int
        Full observation dimension
    coords_p : ndarray[np.float64, ndim=2]
        PE-local array
        Array shape: (nrows, dim_obs_p)
    nrows : int
        Number of rows in array

    Returns
    -------
    coords_f : ndarray[np.float64, ndim=2]
        Full gathered array
        Array shape: (nrows, dim_obs_f)
    status : int
        Status flag: (0) no error
    """


def init(
    filtertype: int,
    subtype: int,
    stepnull: int,
    param_int: NDArray[np.int32],
    dim_pint: int,
    param_real: NDArray[np.float64],
    dim_preal: int,
    comm_model: int,
    comm_filter: int,
    comm_couple: int,
    task_id: int,
    n_modeltasks: int,
    in_filterpe: bool,
    py__init_ens_pdaf: Callable[..., Any],
    in_screen: int,
) -> Tuple[NDArray[np.int32], NDArray[np.float64], int]:
    """Initialise the PDAF system.

    It is called once at the beginning of the assimilation.

    The function specifies the type of DA methods,
    parameters of the filters, the MPI communicators,
    and other parallel options.
    The filter options including `filtertype`, `subtype`,
    `param_int`, and `param_real`
    are introduced in
    `PDAF filter options wiki page <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_.
    Note that the size of `param_int` and `param_real` depends on
    the filter type and subtype. However, for most filters,
    they require at least the state vector size and ensemble size
    for `param_int`, and the forgetting factor for `param_real`.

    The MPI communicators asked by this function depends on
    the parallelisation strategy.
    For the default parallelisation strategy, the user
    can use the parallelisation module
    provided under in `example directory <https://github.com/yumengch/pyPDAF/blob/main/example>`_
    without modifications.
    The parallelisation can differ based on online and offline cases.
    Users can also refer to `parallelisation documentation <https://yumengch.github.io/pyPDAF/parallel.html>`_ for
    explanations or modifications.

    This function also asks for a user-supplied function
    :func:`py__init_ens_pdaf`.
    This function is designed to provides an initial ensemble
    to the internal PDAF ensemble array.
    The internal PDAF ensemble then can be distributed to
    initialise the model forecast using
    :func:`pyPDAF.PDAF.get_state`.
    This user-supplied function can be empty if the model
    has already read the ensemble from restart files.

    Parameters
    ----------
    filtertype : int
        Type of filter
    subtype : int
        Sub-type of filter
    stepnull : int
        Initial time step of assimilation
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameter
    comm_model : int
        Model communicator
    comm_filter : int
        Filter communicator
    comm_couple : int
        Coupling communicator
    task_id : int
        Id of my ensemble task
    n_modeltasks : int
        Number of parallel model tasks
    in_filterpe : bint
        Is my PE a filter-PE?
    py__init_ens_pdaf : Callable
        Initialise ensemble array in PDAF
    in_screen : int
        Control screen output:

    Returns
    -------
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    outflag : int
        Status flag, 0: no error, error codes:
    """


def init_forecast(
    py__next_observation_pdaf: Callable[..., Any],
    py__distribute_state_pdaf: Callable[..., Any],
    py__prepoststep_pdaf: Callable[..., Any],
    outflag: int,
) -> int:
    """The routine PDAF_init_forecast has to be called once at the end of the
    initialization of PDAF/start of DA cycles.

    This function has the purpose to initialize the model fields to be
    propagated from the array holding the ensemble states. In addition,
    the function initializes the information on how many time steps have to be
    performed in the upcoming forecast phase before the next assimilation step,
    and an exit flag indicating whether further model integrations have to be computed.
    These variables are used internally in PDAF and can be retrieved by
    the user by calling PDAF_get_fcst_info.

    Parameters
    ----------
    py__next_observation_pdaf : Callable
        Get the number of time steps to be computed in the forecast phase.
        See details for :func:`pyPDAF.pdaf_c_cb_interface.c__next_observation_pdaf`.
    py__distribute_state_pdaf : Callable
        Distribute a state vector from pdaf to the model/any arrays
        See details for :func:`pyPDAF.pdaf_c_cb_interface.c__distribute_state_pdaf`.
    py__prepoststep_pdaf : Callable
        Process ensemble before or after DA.
        See details for :func:`pyPDAF.pdaf_c_cb_interface.c__prepoststep_pdaf`.
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """


def local_weight(
    wtype: int,
    rtype: int,
    cradius: float,
    sradius: float,
    distance: float,
    nrows: int,
    ncols: int,
    a: NDArray[np.float64],
    var_obs: float,
    verbose: int,
) -> float:
    """Get localisation weight for given distance,
    cut-off radius, support radius, weighting type,
    and weighting function.

    This function is used in the analysis step of a filter
    to computes a localisation weight.

    Typically, in domain-localised filters, the function
    is called in user-supplied :func:`py__prodRinvA_l_pdaf`.
    In LEnKF, this function is called
    in user-supplied :func:`py__localize_covar_pdaf`.

    This function is usually only used without PDAF-OMI.

    Parameters
    ----------
    wtype : int
        type of weight function:
        * `wtype=0`: unit weight
            (`weight=1` up to distance=cradius)
        * `wtype=1`: exponential decrease
            (`weight=1/e` at distance=sradius;
            `weight=0` for distance>cradius)
        * `wtype=2`: 5th order polynomial
            (Gaspari and Cohn 1999; `weight=0` for distance>cradius)
    rtype : int
        type of regulated weighting:
        * `rtype/=1`: no regulation
        * `rtype=1`: regulated by variance of the matrix A and
            the observation variance
    cradius : float
        cut-off radius where weight = 0 beyond the cradius
    sradius : float
        support radius of localisation function. This depends on `wtype`:
        * `wtype=0`: sradius is not used
        * `wtype=1`: weight = :math:`e^{-\frac{distance}{sradius}}`
        * `wtype=2`: weight = 0 if distance > sradius
            else weight = f(distance ,sradius)

        See also: `PDAF-OMI wiki <https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE>`_)
    distance : float
        distance to observation
    nrows : int
        Number of rows in matrix A
    ncols : int
        Number of columns in matrix A
    a : ndarray[tuple[nrows, ncols, ...], np.float64]
        ensemble perturbation/anomaly matrix;
        this matrix is used when weighting is regulated
        by mean variance, i.e., rtype = 1.
        Array shape: (nrows, ncols)
    var_obs : double
        Observation variance
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : double
        localisation weights
    """


def local_weights(
    wtype: int,
    cradius: float,
    sradius: float,
    dim: int,
    distance: NDArray[np.float64],
    verbose: int,
) -> NDArray[np.float64]:
    """Get a vector of localisation weights for given distances,
    cut-off radius, support radius, weighting type,
    and weighting function.

    This function is used in the analysis step of a filter
    to computes a localisation weight.

    Typically, in domain-localised filters, the function
    is called in user-supplied :func:`py__prodRinvA_l_pdaf`.
    In LEnKF, this function is called
    in user-supplied :func:`py__localize_covar_pdaf`.

    This function is usually only used without PDAF-OMI.

    This function is a vectorised version of
    :func:`pyPDAF.PDAF.local_weight` without any regulations.

    Parameters
    ----------
    wtype : int
        type of weight function:
        * `wtype=0`: unit weight
            (`weight=1` up to distance=cradius)
        * `wtype=1`: exponential decrease
            (`weight=1/e` at distance=sradius;
            `weight=0` for distance>cradius)
        * `wtype=2`: 5th order polynomial
            (Gaspari and Cohn 1999; `weight=0` for distance>cradius)
    rtype : int
        type of regulated weighting:
        * `rtype/=1`: no regulation
        * `rtype=1`: regulated by variance of the matrix A and
            the observation variance
    cradius : float
        cut-off radius where weight = 0 beyond the cradius
    sradius : float
        support radius of localisation function. This depends on `wtype`:
        * `wtype=0`: sradius is not used
        * `wtype=1`: weight = :math:`e^{-\frac{distance}{sradius}}`
        * `wtype=2`: weight = 0 if distance > sradius
            else weight = f(distance ,sradius)
    dim : int
        Size of distance and weight arrays
    distance : ndarray[np.float64, ndim=1]
        distance to observation
        Array shape: (dim)
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : ndarray[np.float64, ndim=1]
        Array for localisation weights
        Array shape: (dim)
    """


def print_filter_types(verbose: int) -> None:
    """Print the list of available named DA method types and their IDs.

    Wrapper of the Fortran routine ``PDAF_print_DA_types(verbose)``. If
    ``verbose > 0``, it writes to stdout the names and corresponding integer
    identifiers of supported data assimilation method types (same list as
    for filters). This helps when specifying method types by integer or by
    named parameter.

    The printed names include:
    - PDAF_DA_LESTKF
    - PDAF_DA_ESTKF
    - PDAF_DA_LETKF
    - PDAF_DA_ETKF
    - PDAF_DA_LENKF
    - PDAF_DA_ENKF
    - PDAF_DA_LSEIK
    - PDAF_DA_SEIK
    - PDAF_DA_ENSRF
    - PDAF_DA_LNETF
    - PDAF_DA_NETF
    - PDAF_DA_PF
    - PDAF_DA_LKNETF
    - PDAF_DA_GENOBS
    - PDAF_DA_3DVAR

    Parameters
    ----------
    verbose : int
        Verbosity flag. If 0, no output is printed; if > 0, the list is
        printed to stdout.

    Returns
    -------
    None
    """


def print_da_types(verbose: int) -> None:
    """Print PDAF timing and memory information.

    Wrapper of the Fortran routine ``PDAF_print_info(printtype)``. Prints
    aggregated timing information and/or memory usage depending on
    ``printtype``. Call this near the end of your DA program.

    Parameters
    ----------
    printtype : int
        Type of screen output:
        - 1: general timings
        - 3: timers focused on call-back routines (recommended)
        - 4: detailed timers (analyze filters)
        - 5: very detailed timers (deep filter analysis)
        - 10: allocated memory of the calling MPI task
        - 11: globally used memory (call from all processes)
    """


def print_info(printtype: int) -> None:
    """Print PDAF timing and memory information.

    Wrapper of the Fortran routine ``PDAF_print_info(printtype)``. Prints
    aggregated timing information and/or memory usage depending on
    ``printtype``. Call this near the end of your DA program.

    Parameters
    ----------
    printtype : int
        Type of information to be printed
        * printtype=1: Basic timers
        * printtype=3: Timers showing the time spent in the different call-back routines
            (this variant was added with PDAF 1.15)
        * printtype=4: More detailed timers about parts of the filter algorithm
            (before PDAF 1.15, this was timer level 3)
        * printtype=5: Very detailed timers about various operations in the filter algorithm
            (before PDAF 1.15, this was timer level 4)
        * printtype=10: Memory usage (The value 10 is valid since PDAF V2.1. For older versions use 2)

            - Memory required for the ensemble array,
                state vector, and transform matrix
            - Memory required by the analysis step
            - Memory required to perform the ensemble transformation
    """


def reset_forget(forget_in: float) -> None:
    """Reset the forgetting factor manually
    during the assimilation process.

    For the local ensemble Kalman filters
    the forgetting factor can be set either globally
    if this function is called outside of the loop over
    local domains,
    or
    the forgetting factor can be set differently
    for each local analysis domain within the loop over
    local domains.

    For the LNETF and the global filters
    only a global setting of the forgeting factor is possible.
    In addition, the implementation of adaptive choices
    for the forgetting factor (beyond what is implemented in PDAF) are possible.

    Parameters
    ----------
    forget_in : double
        New value of forgetting factor

    """


def sampleens(
    dim: int,
    dim_ens: int,
    modes: NDArray[np.float64],
    svals: NDArray[np.float64],
    state: NDArray[np.float64],
    verbose: int,
    flag: int,
) -> Tuple[
    NDArray[np.float64],  # modes
    NDArray[np.float64],  # state
    NDArray[np.float64],  # ens
    int,                  # flag
]:
    r"""Generate an ensemble from singular values and
    their vectors (EOF modes) of an ensemble anomaly matrix.

    The singular values and vectors are derived from
    the ensemble anomalies. This ensemble anomaly can be
    obtained from a time anomaly of a model trajectory using
    :func:`pyPDAF.PDAF.eofcovar`.

    Parameters
    ----------
    dim: int
        Size of the state vector
    dim_ens : int
        Ensemble size
    modes : ndarray[tuple[dim, dim_ens-1, ...], np.float64]
        array of EOF modes/matrix of singular vectors.
    svals : ndarray[tuple[dim_ens-1, ...], np.float64]
        singular values.
    state : ndarray[tuple[dim, ...], np.float64]
        PE-local model mean state.
    verbose : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    modes : ndarray[tuple[dim, dim_ens-1, ...], np.float64]
        array of EOF modes/matrix of singular vectors

        The 1st-th dimension dim is size of state vector
    state : ndarray[tuple[dim, ...], np.float64]
        PE-local model mean state

        The array dimension `dim` is size of state vector
    ens : ndarray[tuple[dim, dim_ens, ...], np.float64]
        State ensemble

        The 1st-th dimension dim is size of state vector
        The 2nd-th dimension dim_ens is size of ensemble
    flag : int
        Status flag
    """
