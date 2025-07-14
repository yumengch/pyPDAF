import sys
import numpy as np
cimport numpy as cnp
from . cimport pdaf_c_cb_interface as pdaf_cb
from .cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from .cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from .cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

try:
    import mpi4py
    mpi4py.rc.initialize = False
except ImportError:
    pass

# Global error handler
def global_except_hook(exctype, value, traceback):
    from traceback import print_exception
    try:
        import mpi4py.MPI

        if mpi4py.MPI.Is_initialized():
            try:
                sys.stderr.write('Uncaught exception was ''detected on rank {}.\n'.format(
                    mpi4py.MPI.COMM_WORLD.Get_rank()))
                print_exception(exctype, value, traceback)
                sys.stderr.write("\n")
                sys.stderr.flush()
            finally:
                try:
                    mpi4py.MPI.COMM_WORLD.Abort(1)
                except Exception as e:
                    sys.stderr.write('MPI Abort failed, this process will hang.\n')
                    sys.stderr.flush()
                    raise e
        else:
            sys.__excepthook__(exctype, value, traceback)
    except ImportError:
        sys.__excepthook__(exctype, value, traceback)

sys.excepthook = global_except_hook

def correlation_function(int  ctype, double  length, double  distance):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    ctype : int
        Type of correlation function
    length : double
        Length scale of function
    distance : double
        Distance at which the function is evaluated

    Returns
    -------
    value : double
        Value of the function
    """
    cdef double  value
    with nogil:
        c__pdaf_correlation_function(&ctype, &length, &distance, &value)

    return value


def deallocate():
    """Finalise the PDAF systems
    including freeing some of
    the memory used by PDAF.

    This function cannot be used to
    free all allocated PDAF memory.
    Therefore, one should not use
    :func:`pyPDAF.PDAF.init` afterwards.
    """
    with nogil:
        c__pdaf_deallocate()



def eofcovar(int  dim, int  nstates, int  nfields, int [::1] dim_fields,
    int [::1] offsets, int  remove_mstate, int  do_mv,
    double [::1,:] states, double [::1] meanstate, int  verbose):
    """EOF analysis of an ensemble of state vectors by singular value decomposition.

    Typically, this function is used with
    :func:`pyPDAF.PDAF.SampleEns`
    to generate an ensemble of a chosen size
    (up to the number of EOFs plus one).

    Here, the function performs a singular value decomposition
    of the ensemble anomaly of the input matrix,
    which is usually an ensemble formed by state vectors
    at multiple time steps.
    The singular values and corresponding singular vectors
    can be used to
    construct a covariance matrix.
    This can be used as the initial error covariance
    for the initial ensemble.

    A multivariate scaling can be performed to ensure that all fields in the state vectors have unit variance.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] states_np = np.asarray(states, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] stddev_np = np.zeros((nfields), dtype=np.float64, order="F")
    cdef double [::1] stddev = stddev_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] svals_np = np.zeros((nstates), dtype=np.float64, order="F")
    cdef double [::1] svals = svals_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] svec_np = np.zeros((dim, nstates), dtype=np.float64, order="F")
    cdef double [::1,:] svec = svec_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] meanstate_np = np.asarray(meanstate, dtype=np.float64, order="F")
    cdef int  status
    with nogil:
        c__pdaf_eofcovar(&dim, &nstates, &nfields, &dim_fields[0],
                         &offsets[0], &remove_mstate, &do_mv, &states[0,0],
                         &stddev[0], &svals[0], &svec[0,0], &meanstate[0],
                         &verbose, &status)

    return states_np, stddev_np, svals_np, svec_np, meanstate_np, status


def force_analysis():
    """Perform assimilation after this function call.

    This function overwrite member index of the ensemble state
    by local_dim_ens (number of ensembles for current process,
    in full parallel setup, this is 1.) and the counter
    cnt_steps by nsteps-1.
    This forces that the analysis step is executed at
    the next call to PDAF assimilation functions.
    """
    with nogil:
        c__pdaf_force_analysis()



def gather_dim_obs_f(int  dim_obs_p):
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
    :func:`pyPDAF.PDAF.gather_obs_f` or
    :func:`pyPDAF.PDAF.gather_obs_f2`.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension

    Returns
    -------
    dim_obs_f : int
        Full observation dimension
    """
    cdef int  dim_obs_f
    with nogil:
        c__pdaf_gather_dim_obs_f(&dim_obs_p, &dim_obs_f)

    return dim_obs_f


def gather_obs_f(double [::1] obs_p, int  dimobs_f):
    """In the local filters (LESKTF, LETKF, LSEIK, LNETF) this function returns the total observation vector from process-local observations. The function depends on `pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. Further, the related routine `pyPDAF.PDAF.gather_obs_f2` is used to
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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_np = np.zeros((dimobs_f), dtype=np.float64, order="F")
    cdef double [::1] obs_f = obs_f_np
    cdef int  status
    with nogil:
        c__pdaf_gather_obs_f(&obs_p[0], &obs_f[0], &status)

    return obs_f_np, status


def gather_obs_f2(double [::1,:] coords_p, int  nrows, int  dimobs_f):
    """In the local filters (LESKTF, LETKF, LSEIK, LNETF)
    this function returns the full observation coordinates from process-local observation coordinates. The function depends on `pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. Further, the related routine `pyPDAF.PDAF.gather_obs_f` is used to gather the associated observation vectors.

    The routine is typically used in the routines `py__init_dim_obs_f_pdaf` if the analysis step of the local filters is parallelized.

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
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] coords_f_np = np.zeros((nrows, dimobs_f), dtype=np.float64, order="F")
    cdef double [::1,:] coords_f = coords_f_np
    cdef int  status
    with nogil:
        c__pdaf_gather_obs_f2(&coords_p[0,0], &coords_f[0,0], &nrows, &status)

    return coords_f_np, status


def gather_obs_f_flex(int  dim_obs_p, int  dim_obs_f, double [::1] obs_p):
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
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_np = np.zeros((dim_obs_f), dtype=np.float64, order="F")
    cdef double [::1] obs_f = obs_f_np
    cdef int  status
    with nogil:
        c__pdaf_gather_obs_f_flex(&dim_obs_p, &dim_obs_f, &obs_p[0],
                                  &obs_f[0], &status)

    return obs_f_np, status


def gather_obs_f2_flex(int  dim_obs_p, int  dim_obs_f,
    double [::1,:] coords_p, int  nrows):
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
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] coords_f_np = np.zeros((nrows, dim_obs_f), dtype=np.float64, order="F")
    cdef double [::1,:] coords_f = coords_f_np
    cdef int  status
    with nogil:
        c__pdaf_gather_obs_f2_flex(&dim_obs_p, &dim_obs_f, &coords_p[0,0],
                                   &coords_f[0,0], &nrows, &status)

    return coords_f_np, status


def init(int  filtertype, int  subtype, int  stepnull, int [::1] param_int,
    int  dim_pint, double [::1] param_real, int  dim_preal,
    int  comm_model, int  comm_filter, int  comm_couple, int  task_id,
    int  n_modeltasks, bint  in_filterpe, py__init_ens_pdaf, int  in_screen):
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
        User-supplied routine for ensemble initialization

        Callback Parameters
        -------------------
        filtertype : int
                filter type given in PDAF_init
        dim_p : int
                PE-local state dimension given by PDAF_init
        dim_ens : int
                number of ensemble members
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                This array must be filled with the initial
                state of the model for SEEK, but it is not
                used for ensemble-based filters.
                One can still make use of this array within
                this function.
                Array shape: (dim_p)
        uinv : ndarray[np.float64, ndim=2]
                This array is the inverse of matrix
                formed by right singular vectors of error
                covariance matrix of ensemble perturbations.
                This array has to be filled in SEEK, but it is
                not used for ensemble-based filters.
                Nevertheless, one can still make use of this
                array within this function e.g.,
                for generating an initial ensemble perturbation
                from a given covariance matrix.
                Dimension of this array is determined by the
                filter type.
                - (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
                - (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
                - (1, 1) for (L)EnKF, particle filters and gen_obs
                Array shape: (dim_ens - 1, dim_ens-1)
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        flag : int
                pdaf status flag

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                This array must be filled with the initial
                state of the model for SEEK, but it is not
                used for ensemble-based filters.
                One can still make use of this array within
                this function.
                Array shape: (dim_p)
        uinv : ndarray[np.float64, ndim=2]
                This array is the inverse of matrix
                formed by right singular vectors of error
                covariance matrix of ensemble perturbations.
                This array has to be filled in SEEK, but it is
                not used for ensemble-based filters.
                Nevertheless, one can still make use of this
                array within this function e.g.,
                for generating an initial ensemble perturbation
                from a given covariance matrix.
                Dimension of this array is determined by the
                filter type.
                - (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
                - (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
                - (1, 1) for (L)EnKF, particle filters and gen_obs
                Array shape: (dim_ens - 1, dim_ens-1)
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        flag : int
                pdaf status flag

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
    pdaf_cb.init_ens_pdaf = <void*>py__init_ens_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_init(&filtertype, &subtype, &stepnull, &param_int[0],
                     &dim_pint, &param_real[0], &dim_preal, &comm_model,
                     &comm_filter, &comm_couple, &task_id, &n_modeltasks,
                     &in_filterpe, pdaf_cb.c__init_ens_pdaf, &in_screen,
                     &outflag)

    return param_int_np, param_real_np, outflag


def init_forecast(py__next_observation_pdaf, py__distribute_state_pdaf,
    py__prepoststep_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__next_observation_pdaf : Callable
        Provide information on next forecast

        Callback Parameters
        -------------------
        stepnow : int
                the current time step given by PDAF

        Callback Returns
        ----------------
        nsteps : int
                number of forecast time steps until next assimilation;
                this can also be interpreted as
                number of assimilation function calls
                to perform a new assimilation
        doexit : int
                whether to exit forecasting (1 for exit)
        time : double
                current model (physical) time

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

        Callback Parameters
        -------------------
        dim_p : int
                PE-local state dimension
        state_p : ndarray[np.float64, ndim=1]
                PE-local state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local state vector
                Array shape: (dim_p)

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

        Callback Parameters
        -------------------
        step : int
                current time step
                (negative for call before analysis/preprocessing)
        dim_p : int
                PE-local state vector dimension
        dim_ens : int
                number of ensemble members
        dim_ens_l : int
                number of ensemble members run serially
                on each model task
        dim_obs_p : int
                PE-local dimension of observation vector
        state_p : ndarray[np.float64, ndim=1]
                pe-local forecast/analysis state
                (the array 'state_p' is generally not
                initialised in the case of ESTKF/ETKF/EnKF/SEIK,
                so it can be used freely here.)
                Array shape: (dim_p)
        uinv : ndarray[np.float64, ndim=2]
                Inverse of the transformation matrix in ETKF and ESKTF;
                inverse of matrix formed by right singular vectors of error
                covariance matrix of ensemble perturbations in SEIK/SEEK.
                not used in EnKF.
                Array shape: (dim_ens-1, dim_ens-1)
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        flag : int
                pdaf status flag

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local forecast/analysis state
                (the array 'state_p' is generally not
                initialised in the case of ESTKF/ETKF/EnKF/SEIK,
                so it can be used freely here.)
                Array shape: (dim_p)
        uinv : ndarray[np.float64, ndim=2]
                Inverse of the transformation matrix in ETKF and ESKTF;
                inverse of matrix formed by right singular vectors of error
                covariance matrix of ensemble perturbations in SEIK/SEEK.
                not used in EnKF.
                Array shape: (dim_ens-1, dim_ens-1)
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf_init_forecast(pdaf_cb.c__next_observation_pdaf,
                              pdaf_cb.c__distribute_state_pdaf,
                              pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def local_weight(int  wtype, int  rtype, double  cradius, double  sradius,
    double  distance, int  nrows, int  ncols, double [::1,:] a,
    double  var_obs, int  verbose):
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
        Type of weight function
    rtype : int
        Type of regulated weighting
    cradius : double
        Cut-off radius
    sradius : double
        Support radius
    distance : double
        Distance to observation
    nrows : int
        Number of rows in matrix A
    ncols : int
        Number of columns in matrix A
    a : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (nrows, ncols)
    var_obs : double
        Observation variance
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : double
        Weights
    """
    cdef double  weight
    with nogil:
        c__pdaf_local_weight(&wtype, &rtype, &cradius, &sradius, &distance,
                             &nrows, &ncols, &a[0,0], &var_obs, &weight,
                             &verbose)

    return weight


def local_weights(int  wtype, double  cradius, double  sradius, int  dim,
    double [::1] distance, int  verbose):
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
        Type of weight function
    cradius : double
        Parameter for cut-off
    sradius : double
        Support radius
    dim : int
        Size of distance and weight arrays
    distance : ndarray[np.float64, ndim=1]
        Array holding distances
        Array shape: (dim)
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : ndarray[np.float64, ndim=1]
        Array for weights
        Array shape: (dim)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] weight_np = np.zeros((dim), dtype=np.float64, order="F")
    cdef double [::1] weight = weight_np
    with nogil:
        c__pdaf_local_weights(&wtype, &cradius, &sradius, &dim,
                              &distance[0], &weight[0], &verbose)

    return weight_np


def print_filter_types(int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    verbose : int


    Returns
    -------
    """
    with nogil:
        c__pdaf_print_filter_types(&verbose)



def print_da_types(int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    verbose : int


    Returns
    -------
    """
    with nogil:
        c__pdaf_print_da_types(&verbose)



def print_info(int  printtype):
    """Print the wallclock time and memory measured by PDAF.

    This should be called at the end of the DA program.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_print_info(&printtype)



def reset_forget(double  forget_in):
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

    Returns
    -------
    """
    with nogil:
        c__pdaf_reset_forget(&forget_in)



def sampleens(int  dim, int  dim_ens, double [::1,:] modes,
    double [::1] svals, double [::1] state, int  verbose, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        Size of state vector
    dim_ens : int
        Size of ensemble
    modes : ndarray[np.float64, ndim=2]
        Array of EOF modes
        Array shape: (dim, dim_ens-1)
    svals : ndarray[np.float64, ndim=1]
        Vector of singular values
        Array shape: (dim_ens-1)
    state : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim)
    verbose : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    modes : ndarray[np.float64, ndim=2]
        Array of EOF modes
        Array shape: (dim, dim_ens-1)
    state : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        State ensemble
        Array shape: (dim, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] modes_np = np.asarray(modes, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.zeros((dim, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ens = ens_np
    with nogil:
        c__pdaf_sampleens(&dim, &dim_ens, &modes[0,0], &svals[0],
                          &state[0], &ens[0,0], &verbose, &flag)

    return modes_np, state_np, ens_np, flag


