import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb

def init(int  filtertype, int  subtype, int  stepnull, int [::1] param_int,
    int  dim_pint, double [::1] param_real, int  dim_preal,
    py__init_ens_pdaf, int  in_screen):
    """init(filtertype: int, subtype: int, stepnull: int, param_int: np.ndarray, dim_pint: int, param_real: np.ndarray, dim_preal: int, py__init_ens_pdaf: Callable, in_screen: int) -> Tuple[np.ndarray, np.ndarray, int]

    Initialise the PDAF system.

    It is called once at the beginning of the assimilation.
    The function has to be used in tandem with :func:`pyPDAF.PDAF3.set_parallel`.

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
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    pdaf_cb.init_ens_pdaf = py__init_ens_pdaf
    cdef int  outflag
    c__pdaf3_init(&filtertype, &subtype, &stepnull, &param_int[0],
                     &dim_pint, &param_real[0], &dim_preal,
                     pdaf_cb.c__init_ens_pdaf, &in_screen, &outflag)

    return param_int_np, param_real_np, outflag

def init_forecast(py__next_observation_pdaf, py__distribute_state_pdaf,
    py__prepoststep_pdaf, int  outflag):
    """init_forecast(py__next_observation_pdaf: Callable, py__distribute_state_pdaf: Callable, py__prepoststep_pdaf: Callable, outflag: int) -> int

    The routine PDAF_init_forecast has to be called once at the end of the
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
    pdaf_cb.next_observation_pdaf = py__next_observation_pdaf
    pdaf_cb.distribute_state_pdaf = py__distribute_state_pdaf
    pdaf_cb.prepoststep_pdaf = py__prepoststep_pdaf
    c__pdaf3_init_forecast(pdaf_cb.c__next_observation_pdaf,
                              pdaf_cb.c__distribute_state_pdaf,
                              pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag

def init_parallel(int screen, int type_parallel, int online_coupling,
    int dim_ens, int n_modeltasks, int COMM_model):
    """init_parallel(screen: int, type_parallel: int, online_coupling: int, dim_ens: int, n_modeltasks: int, COMM_model: int) -> Tuple[int, int, int, int, int, int, int, int]

    Initialize PDAF3's flexible parallelization layout.

    This helper delegates communicator setup to PDAF3. It is intended for
    applications that use the PDAF3 interface and want PDAF to derive the model
    and assimilation communicators from the selected parallelization strategy.
    Call it during PDAF3 setup before routines that rely on PDAF's internal
    parallel configuration.

    The configuration of the parallelization is what we found to be
    'typically working'. If one has a particular case, e.g., the ensemble
    parallelization is readily prepared by the model, one can do the
    communicator configuration in the user code and then call
    :func:`pyPDAF.set_parallel` to provide PDAF with the configuration.


    Parameters
    ----------
    screen : int
        Verbosity flag for PDAF screen output.
    type_parallel : int
        PDAF3 parallelization strategy selector. (currently not used).

        Type of parallelization:
            - 0: common setup using task 1 for assimilation
            - 1: setup using separate task 0 for assimilation
    online_coupling : int
        Coupling mode flag. Use the value expected by PDAF for online or
        offline coupling.
        1: online DA coupling, 0: offline DA coupling
    dim_ens : int
        Ensemble size.
    n_modeltasks : int
        Number of parallel model tasks supplied to PDAF3.
        If online coupling is used, this value must match the number of model tasks
        that are running the model code. If offline coupling is used, this value
        is always 1 regardless of the input.
    COMM_model : int
        MPI communicator handle. This should be the communicator that contains
        all processors used to run the entire ensemble. This can be MPI_COMM_WORLD.
        In some cases, e.g., NEMO, the processors are split into model and XIOS
        communicators. In that case, this should be model communicator.

    Returns
    -------
    n_modeltasks : int
        Number of parallel model tasks supplied to PDAF3.
    COMM_model : int
        This will be the model communicator after splitting by the number
        of model tasks. This communicator contains processors that run a single
        model task instead of the entire ensemble. This is different from the input
        argument ``COMM_model``.
    mype_model : int
        Rank of the current process in ``COMM_model``.
    npes_model : int
        Number of processes in ``COMM_model``.
    COMM_assim : int
        MPI communicator contains all processors that run actual DA filters.
    mype_assim : int
        Rank of the current process in ``COMM_assim``.
    npes_assim : int
        Number of processes in ``COMM_assim``.
    task_id : int
        Ensemble task identifier of the current process.
    """
    cdef int mype_model
    cdef int npes_model
    cdef int COMM_assim
    cdef int mype_assim
    cdef int npes_assim
    cdef int task_id
    c__pdaf3_init_parallel(&screen, &type_parallel, &online_coupling,
        &dim_ens, &n_modeltasks, &COMM_model, &mype_model, &npes_model,
        &COMM_assim, &mype_assim, &npes_assim, &task_id)

    return n_modeltasks, COMM_model, mype_model, npes_model, COMM_assim, \
            mype_assim, npes_assim, task_id

def set_parallel(int in_comm_pdaf,
    int in_comm_model, int in_comm_filter, int in_comm_couple, int in_task_id,
    int in_n_modeltasks, bint in_filterpe, int flag):
    """set_parallel(in_comm_pdaf: int, in_comm_model: int, in_comm_filter: int, in_comm_couple: int, in_task_id: int, in_n_modeltasks: int, in_filterpe: bool, flag: int) -> int

    Set MPI communicators and parallelisation in PDAF.

    The MPI communicators asked by this function depends on
    the parallelisation strategy.
    For the default parallelisation strategy, the user
    can use the parallelisation module
    provided under in `example directory <https://github.com/yumengch/pyPDAF/blob/main/example>`_
    without modifications.
    The parallelisation can differ based on online and offline cases.
    Users can also refer to `parallelisation documentation <https://yumengch.github.io/pyPDAF/parallel.html>`_ for
    explanations or modifications.

    Parameters
    ----------
    in_comm_model : int
        Model communicator
    in_comm_filter : int
        Filter communicator
    in_comm_couple : int
        Coupling communicator
    in_task_id : int
        Id of my ensemble task
    in_n_modeltasks : int
        Number of parallel model tasks
    in_filterpe : bint
        Is my PE a filter-PE?
    flag: int
        Status flag

    Returns
    -------
    flag : int
        Status flag, 0: no error, error codes:
    """
    c__pdaf3_set_parallel(&in_comm_pdaf, &in_comm_model, &in_comm_filter,
                     &in_comm_couple, &in_task_id, &in_n_modeltasks,
                     &in_filterpe, &flag)

    return flag
