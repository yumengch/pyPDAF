import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

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

def mpi_init():
    """Initialise MPI
    """
    with nogil:
        c__pdaf_mpi_init()

def timeit(int timerid, str operation):
    """PDAF timer

    Parameters
    ----------
    timerid : int
        Timer ID
    operation : str
        Operation name
    """
    operation_byte = operation.encode('UTF-8')
    cdef char* operation_ptr = operation_byte
    with nogil:
        c__pdaf_timeit(&timerid, operation_ptr)

def set_forget(int  step, int  localfilter, int  dim_obs_p, int  dim_ens,
    double [::1,:] mens_p, double [::1] mstate_p, double [::1] obs_p,
    py__init_obsvar_pdaf, double  forget_in, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    localfilter : int
        Whether filter is domain-local
    dim_obs_p : int
        Dimension of observation vector
    dim_ens : int
        Ensemble size
    mens_p : ndarray[np.float64, ndim=2]
        Observed PE-local ensemble
        Array shape: (dim_obs_p, dim_ens)
    mstate_p : ndarray[np.float64, ndim=1]
        Observed PE-local mean state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        Observation vector
        Array shape: (dim_obs_p)
    py__init_obsvar_pdaf : Callable
        Initialize mean obs. error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    forget_in : double
        Prescribed forgetting factor
    screen : int
        Verbosity flag

    Returns
    -------
    forget_out : double
        Adaptively estimated forgetting factor
    """
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef double  forget_out
    with nogil:
        c__pdaf_set_forget(&step, &localfilter, &dim_obs_p, &dim_ens,
                           &mens_p[0,0], &mstate_p[0], &obs_p[0],
                           pdaf_cb.c__init_obsvar_pdaf, &forget_in,
                           &forget_out, &screen)

    return forget_out

def set_iparam_filters(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_set_iparam_filters(&id, &value, &flag)

    return flag


def set_rparam_filters(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_set_rparam_filters(&id, &value, &flag)

    return flag

def set_forget_local(int  domain, int  step, int  dim_obs_l, int  dim_ens,
    double [::1,:] hx_l, double [::1] hxbar_l, double [::1] obs_l,
    py__init_obsvar_l_pdaf, double  forget):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Dimension of local observation vector
    dim_ens : int
        Ensemble size
    hx_l : ndarray[np.float64, ndim=2]
        Local observed ensemble
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed state estimate
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    py__init_obsvar_l_pdaf : Callable
        Initialize local mean obs. error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    forget : double
        Prescribed forgetting factor

    Returns
    -------
    aforget : double
        Adaptive forgetting factor
    """
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef double  aforget
    with nogil:
        c__pdaf_set_forget_local(&domain, &step, &dim_obs_l, &dim_ens,
                                 &hx_l[0,0], &hxbar_l[0], &obs_l[0],
                                 pdaf_cb.c__init_obsvar_l_pdaf, &forget,
                                 &aforget)

    return aforget


def fcst_operations(int  step, py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Time step in current forecast phase
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

        Callback Parameters
        -------------------
        dim_p : int
                pe-local state dimension
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize PE-local observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    with nogil:
        c__pdaf_fcst_operations(&step, pdaf_cb.c__collect_state_pdaf,
                                pdaf_cb.c__distribute_state_pdaf,
                                pdaf_cb.c__init_dim_obs_pdaf,
                                pdaf_cb.c__obs_op_pdaf,
                                pdaf_cb.c__init_obs_pdaf, &outflag)

    return outflag


def letkf_ana_t(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, double [::1] state_l, double [::1,:] ens_l,
    double [::1,:] hz_l, double [::1] hxbar_l, double [::1] obs_l,
    double [::1,:] rndmat, double  forget, py__prodrinva_l_pdaf,
    int  type_trans, int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    type_trans : int
        Type of ensemble transformation
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        on exit: local weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv_l = ainv_l_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_l_np = np.asarray(hz_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] rndmat_np = np.asarray(rndmat, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_letkf_ana_t(&domain_p, &step, &dim_l, &dim_obs_l, &dim_ens,
                            &state_l[0], &ainv_l[0,0], &ens_l[0,0],
                            &hz_l[0,0], &hxbar_l[0], &obs_l[0],
                            &rndmat[0,0], &forget,
                            pdaf_cb.c__prodrinva_l_pdaf, &type_trans,
                            &screen, &debug, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hz_l_np, rndmat_np, forget, flag


def seik_update(int  step, int  dim_p, int  dim_ens, int  rank,
    double [::1] state_p, double [::1,:] uinv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prodrinva_pdaf, py__init_obsvar_pdaf, py__prepoststep_pdaf,
    int  screen, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for SEIK analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafseik_update(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                           &state_p[0], &uinv[0,0], &ens_p[0,0],
                           pdaf_cb.c__init_dim_obs_pdaf,
                           pdaf_cb.c__obs_op_pdaf,
                           pdaf_cb.c__init_obs_pdaf,
                           pdaf_cb.c__prodrinva_pdaf,
                           pdaf_cb.c__init_obsvar_pdaf,
                           pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                           &flag)

    return dim_obs_p, state_p_np, uinv_np, ens_p_np, flag


def _3dvar_update(int  step, int  dim_p, int  dim_ens, int  dim_cvec,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prodrinva_pdaf, py__prepoststep_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_cvec : int
        Size of control vector (parameterized part)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Not used in 3D-Var
        Array shape: (1, 1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for 3DVAR analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

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

    py__cvt_pdaf : Callable
        Apply control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Not used in 3D-Var
        Array shape: (1, 1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdaf3dvar_update(&step, &dim_p, &dim_obs_p, &dim_ens, &dim_cvec,
                            &state_p[0], &ainv[0,0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__prodrinva_pdaf,
                            pdaf_cb.c__prepoststep_pdaf,
                            pdaf_cb.c__cvt_pdaf, pdaf_cb.c__cvt_adj_pdaf,
                            pdaf_cb.c__obs_op_lin_pdaf,
                            pdaf_cb.c__obs_op_adj_pdaf, &screen,
                            &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag


def en3dvar_update_estkf(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec_ens, double [::1] state_p, double [::1,:] ainv,
    double [::1,:] ens_p, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__prepoststep_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_obsvar_pdaf, int  screen,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_cvec_ens : int
        Size of control vector (ensemble part)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for 3DVAR analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

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

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble var)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafen3dvar_update_estkf(&step, &dim_p, &dim_obs_p, &dim_ens,
                                    &dim_cvec_ens, &state_p[0], &ainv[0,0],
                                    &ens_p[0,0],
                                    pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__prodrinva_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__cvt_ens_pdaf,
                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf,
                                    pdaf_cb.c__init_obsvar_pdaf, &screen,
                                    &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag


def en3dvar_update_lestkf(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec_ens, double [::1] state_p, double [::1,:] ainv,
    double [::1,:] ens_p, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__prepoststep_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf,
    py__init_obs_f_pdaf, py__init_obs_l_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_cvec_ens : int
        Size of control vector (ensemble part)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for 3DVAR analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

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

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble var)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    py__init_dim_obs_f_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_f_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_f_pdaf : Callable
        Initialize PE-local observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_f : int
                Size of the full observation vector

        Callback Returns
        ----------------
        observation_f : ndarray[np.float64, ndim=1]
                Full vector of observations
                Array shape: (dim_obs_f)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from full state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    pdaf_cb.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    pdaf_cb.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafen3dvar_update_lestkf(&step, &dim_p, &dim_obs_p, &dim_ens,
                                     &dim_cvec_ens, &state_p[0],
                                     &ainv[0,0], &ens_p[0,0],
                                     pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__init_obs_pdaf,
                                     pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     pdaf_cb.c__init_dim_obs_f_pdaf,
                                     pdaf_cb.c__obs_op_f_pdaf,
                                     pdaf_cb.c__init_obs_f_pdaf,
                                     pdaf_cb.c__init_obs_l_pdaf,
                                     pdaf_cb.c__prodrinva_l_pdaf,
                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                     pdaf_cb.c__init_dim_l_pdaf,
                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                     pdaf_cb.c__g2l_state_pdaf,
                                     pdaf_cb.c__l2g_state_pdaf,
                                     pdaf_cb.c__g2l_obs_pdaf,
                                     pdaf_cb.c__init_obsvar_pdaf,
                                     pdaf_cb.c__init_obsvar_l_pdaf,
                                     &screen, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag


def etkf_update(int  step, int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__prodrinva_pdaf,
    py__init_obsvar_pdaf, py__prepoststep_pdaf, int  screen, int  subtype,
    int  dim_lag, double [::1,:,:] sens_p, int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for ETKF analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafetkf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                           &state_p[0], &ainv[0,0], &ens_p[0,0],
                           pdaf_cb.c__init_dim_obs_pdaf,
                           pdaf_cb.c__obs_op_pdaf,
                           pdaf_cb.c__init_obs_pdaf,
                           pdaf_cb.c__prodrinva_pdaf,
                           pdaf_cb.c__init_obsvar_pdaf,
                           pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                           &dim_lag, &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def netf_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1,:] ens_p, double [::1,:] rndmat, double [::1,:] t,
    int  type_forget, double  forget, int  type_winf, double  limit_winf,
    int  type_noise, double  noise_amp, double [::1,:] hz_p,
    double [::1] obs_p, py__likelihood_pdaf, int  screen, int  debug,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Orthogonal random matrix
        Array shape: (dim_ens, dim_ens)
    t : ndarray[np.float64, ndim=2]
        Ensemble transform matrix
        Array shape: (dim_ens, dim_ens)
    type_forget : int
        Type of forgetting factor
    forget : double
        Forgetting factor
    type_winf : int
        Type of weights inflation
    limit_winf : double
        Limit for weights inflation
    type_noise : int
        Type of pertubing noise
    noise_amp : double
        Amplitude of noise
    hz_p : ndarray[np.float64, ndim=2]
        Temporary matrices for analysis
        Array shape: (dim_obs_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__likelihood_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        resid : ndarray[np.float64, ndim=1]
                Input vector holding the residual
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        likely : double
                Output value of the likelihood

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local forecast state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    t : ndarray[np.float64, ndim=2]
        Ensemble transform matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] t_np = np.asarray(t, dtype=np.float64, order="F")
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    with nogil:
        c__pdaf_netf_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &state_p[0],
                         &ens_p[0,0], &rndmat[0,0], &t[0,0], &type_forget,
                         &forget, &type_winf, &limit_winf, &type_noise,
                         &noise_amp, &hz_p[0,0], &obs_p[0],
                         pdaf_cb.c__likelihood_pdaf, &screen, &debug, &flag)

    return state_p_np, ens_p_np, t_np, flag


def netf_smoothert(int  step, int  dim_p, int dim_obs_p, int  dim_ens,
    double [::1,:] ens_p, double [::1,:] rndmat, double [::1,:] ta,
    double [::1,:] hx_p, double[::1] obs_p,
    py__likelihood_pdaf, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p: int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Orthogonal random matrix
        Array shape: (dim_ens, dim_ens)
    ta : ndarray[np.float64, ndim=2]
        Ensemble transform matrix
        Array shape: (dim_ens, dim_ens)
    hx_p : ndarray[np.float64, ndim=1]
        Temporary matrices for analysis
        Array shape: (dim_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__likelihood_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        resid : ndarray[np.float64, ndim=1]
                Input vector holding the residual
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        likely : double
                Output value of the likelihood

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    t : ndarray[np.float64, ndim=2]
        Ensemble transform matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ta_np = np.asarray(ta, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hx_p_np = np.asarray(hx_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_p_np = np.asarray(obs_p, dtype=np.float64, order="F")

    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf

    with nogil:
        c__pdaf_netf_smoothert(&step, &dim_p, &dim_obs_p, &dim_ens,
                               &ens_p[0,0], &rndmat[0,0], &ta[0,0],
                               &hx_p[0,0], &obs_p[0],
                               pdaf_cb.c__likelihood_pdaf, &screen, &flag)

    return ens_p_np, ta_np, flag


def smoother_netf(int  dim_p, int  dim_ens, int  dim_lag,
    double [::1,:] ainv, double [::1,:,:] sens_p, int  cnt_maxlag, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ainv : ndarray[np.float64, ndim=2]
        Weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    screen : int
        Verbosity flag

    Returns
    -------
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_smoother_netf(&dim_p, &dim_ens, &dim_lag, &ainv[0,0],
                              &sens_p[0,0,0], &cnt_maxlag, &screen)

    return sens_p_np, cnt_maxlag


def lnetf_ana(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, double [::1,:] ens_l, double [::1,:] hx_l,
    double [::1] obs_l, double [::1,:] rndmat, py__likelihood_l_pdaf,
    int  type_forget, double  forget, int  type_winf, double  limit_winf,
    int  cnt_small_svals, double [::1] eff_dimens, double [::1,:] t,
    int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hx_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    py__likelihood_l_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    type_forget : int
        Typ eof forgetting factor
    forget : double
        Forgetting factor
    type_winf : int
        Type of weights inflation
    limit_winf : double
        Limit for weights inflation
    cnt_small_svals : int
        Number of small eigen values
    eff_dimens : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    t : ndarray[np.float64, ndim=2]
        local ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    cnt_small_svals : int
        Number of small eigen values
    eff_dimens : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    t : ndarray[np.float64, ndim=2]
        local ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] eff_dimens_np = np.asarray(eff_dimens, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] t_np = np.asarray(t, dtype=np.float64, order="F")
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    with nogil:
        c__pdaf_lnetf_ana(&domain_p, &step, &dim_l, &dim_obs_l, &dim_ens,
                          &ens_l[0,0], &hx_l[0,0], &obs_l[0], &rndmat[0,0],
                          pdaf_cb.c__likelihood_l_pdaf, &type_forget,
                          &forget, &type_winf, &limit_winf,
                          &cnt_small_svals, &eff_dimens[0], &t[0,0],
                          &screen, &debug, &flag)

    return ens_l_np, cnt_small_svals, eff_dimens_np, t_np, flag


def lnetf_smoothert(int  domain_p, int  step, int  dim_obs_f,
    int  dim_obs_l, int  dim_ens, double [::1,:] hx_f,
    double [::1,:] rndmat, py__g2l_obs_pdaf, py__init_obs_l_pdaf,
    py__likelihood_l_pdaf, int  screen, double [::1,:] t, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_f : int
        PE-local dimension of full observation vector
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    hx_f : ndarray[np.float64, ndim=2]
        PE-local full observed state ens.
        Array shape: (dim_obs_f, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__likelihood_l_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    screen : int
        Verbosity flag
    t : ndarray[np.float64, ndim=2]
        local ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag

    Returns
    -------
    t : ndarray[np.float64, ndim=2]
        local ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] t_np = np.asarray(t, dtype=np.float64, order="F")
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    with nogil:
        c__pdaf_lnetf_smoothert(&domain_p, &step, &dim_obs_f, &dim_obs_l,
                                &dim_ens, &hx_f[0,0], &rndmat[0,0],
                                pdaf_cb.c__g2l_obs_pdaf,
                                pdaf_cb.c__init_obs_l_pdaf,
                                pdaf_cb.c__likelihood_l_pdaf, &screen,
                                &t[0,0], &flag)

    return t_np, flag


def smoother_lnetf(int  domain_p, int  step, int  dim_p, int  dim_l,
    int  dim_ens, int  dim_lag, double [::1,:] ainv, double [::1,:] ens_l,
    double [::1,:,:] sens_p, int  cnt_maxlag, py__g2l_state_pdaf,
    py__l2g_state_pdaf, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_l : int
        State dimension on local analysis domain
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ainv : ndarray[np.float64, ndim=2]
        Weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        local past ensemble (temporary)
        Array shape: (dim_l, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    screen : int
        Verbosity flag

    Returns
    -------
    ens_l : ndarray[np.float64, ndim=2]
        local past ensemble (temporary)
        Array shape: (dim_l, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    with nogil:
        c__pdaf_smoother_lnetf(&domain_p, &step, &dim_p, &dim_l, &dim_ens,
                               &dim_lag, &ainv[0,0], &ens_l[0,0],
                               &sens_p[0,0,0], &cnt_maxlag,
                               pdaf_cb.c__g2l_state_pdaf,
                               pdaf_cb.c__l2g_state_pdaf, &screen)

    return ens_l_np, sens_p_np, cnt_maxlag


def memcount_ini(int  ncounters):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    ncounters : int
        Number of memory counters

    Returns
    -------
    """
    with nogil:
        c__pdaf_memcount_ini(&ncounters)



def memcount_define(str  stortype, int  wordlength):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    stortype : char
        Type of variable
    wordlength : int
        Word length for chosen type

    Returns
    -------
    """
    stortype_byte = stortype.encode('UTF-8')
    cdef char*  stortype_ptr = stortype_byte
    with nogil:
        c__pdaf_memcount_define(stortype_ptr, &wordlength)



def memcount(int  id, str stortype, int  dim):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Id of the counter
    stortype : char
        Type of variable
    dim : int
        Dimension of allocated variable

    Returns
    -------
    """
    stortype_byte = stortype.encode('UTF-8')
    cdef char*  stortype_ptr = stortype_byte
    with nogil:
        c__pdaf_memcount(&id, stortype_ptr, &dim)



def init_filters(int  type_filter, int  subtype, int [::1] param_int,
    int  dim_pint, double [::1] param_real, int  dim_preal, int  screen,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    type_filter : int
        Type of filter
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    screen : int
        Control screen output
    flag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    filterstr : char
        Name of filter algorithm
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef char*  filterstr_ptr
    cdef str  filterstr
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_init_filters(&type_filter, &subtype, &param_int[0],
                             &dim_pint, &param_real[0], &dim_preal,
                             filterstr_ptr, &ensemblefilter, &fixedbasis,
                             &screen, &flag)
    filterstr = filterstr_ptr.decode('UTF-8')
    return subtype, param_int_np, param_real_np, filterstr, ensemblefilter, fixedbasis, flag


def alloc_filters(str  filterstr, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    filterstr : char
        Name of filter algorithm
    subtype : int
        Sub-type of filter
    flag : int
        Status flag

    Returns
    -------
    flag : int
        Status flag
    """
    filterstr_byte = filterstr.encode('UTF-8')
    cdef char* filterstr_ptr = filterstr_byte
    with nogil:
        c__pdaf_alloc_filters(filterstr_ptr, &subtype, &flag)

    return flag


def configinfo_filters(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_configinfo_filters(&subtype, &verbose)

    return subtype


def options_filters(int  type_filter):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    type_filter : int
        Type of filter

    Returns
    -------
    """
    with nogil:
        c__pdaf_options_filters(&type_filter)



def print_info_filters(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_print_info_filters(&printtype)

def allreduce(int  val_p, int  mpitype, int  mpiop):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    val_p : int
        PE-local value
    mpitype : int
        MPI data type
    mpiop : int
        MPI operator

    Returns
    -------
    val_g : int
        reduced global value
    status : int
        Status flag: (0) no error
    """
    cdef int  val_g
    cdef int  status
    with nogil:
        c__pdaf_allreduce(&val_p, &val_g, &mpitype, &mpiop, &status)

    return val_g, status


def lseik_update(int  step, int  dim_p, int  dim_ens, int  rank,
    double [::1] state_p, double [::1,:] uinv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__prepoststep_pdaf, int  screen,
    int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Compute product of R^(-1) with HV

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaflseik_update(&step, &dim_p, &dim_obs_f, &dim_ens, &rank,
                            &state_p[0], &uinv[0,0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__init_obs_l_pdaf,
                            pdaf_cb.c__prodrinva_l_pdaf,
                            pdaf_cb.c__init_n_domains_p_pdaf,
                            pdaf_cb.c__init_dim_l_pdaf,
                            pdaf_cb.c__init_dim_obs_l_pdaf,
                            pdaf_cb.c__g2l_state_pdaf,
                            pdaf_cb.c__l2g_state_pdaf,
                            pdaf_cb.c__g2l_obs_pdaf,
                            pdaf_cb.c__init_obsvar_pdaf,
                            pdaf_cb.c__init_obsvar_l_pdaf,
                            pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                            &flag)

    return dim_obs_f, state_p_np, uinv_np, ens_p_np, flag


def ensrf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_ensrf_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def ensrf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_ensrf_alloc(&outflag)

    return outflag


def ensrf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_ensrf_config(&subtype, &verbose)

    return subtype


def ensrf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_ensrf_set_iparam(&id, &value, &flag)

    return flag


def ensrf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_ensrf_set_rparam(&id, &value, &flag)

    return flag


def ensrf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_ensrf_options()



def ensrf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_ensrf_memtime(&printtype)



def estkf_ana_fixed(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    int  rank, double [::1] state_p, double [::1,:] ainv,
    double [::1,:] ens_p, double [::1,:] hl_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  type_sqrt, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast mean state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix A - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 with some matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    type_sqrt : int
        Type of square-root of A
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast mean state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix A - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_p_np = np.asarray(hl_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_estkf_ana_fixed(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                                &state_p[0], &ainv[0,0], &ens_p[0,0],
                                &hl_p[0,0], &hxbar_p[0], &obs_p[0],
                                &forget, pdaf_cb.c__prodrinva_pdaf,
                                &screen, &type_sqrt, &debug, &flag)

    return state_p_np, ainv_np, ens_p_np, hl_p_np, flag


def etkf_ana_fixed(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1,:] ens_p, double [::1,:] hz_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        on exit: weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv = ainv_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_p_np = np.asarray(hz_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_etkf_ana_fixed(&step, &dim_p, &dim_obs_p, &dim_ens,
                               &state_p[0], &ainv[0,0], &ens_p[0,0],
                               &hz_p[0,0], &hxbar_p[0], &obs_p[0], &forget,
                               pdaf_cb.c__prodrinva_pdaf, &screen, &debug,
                               &flag)

    return state_p_np, ainv_np, ens_p_np, hz_p_np, flag


def estkf_update(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prodrinva_pdaf, py__init_obsvar_pdaf, py__prepoststep_pdaf,
    int  screen, int  subtype, int  envar_mode, int  dim_lag,
    double [::1,:,:] sens_p, int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of transform matrix A
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for ESTKF analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    envar_mode : int
        Flag whether routine is called from 3DVar for special functionality
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of transform matrix A
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafestkf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                            &state_p[0], &ainv[0,0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__prodrinva_pdaf,
                            pdaf_cb.c__init_obsvar_pdaf,
                            pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                            &envar_mode, &dim_lag, &sens_p[0,0,0],
                            &cnt_maxlag, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def lknetf_update_step(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_hyb_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__likelihood_l_pdaf, py__likelihood_hyb_l_pdaf, py__prepoststep_pdaf,
    int  screen, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_hyb_l_pdaf : Callable
        Compute product of R^(-1) with HV with hybrid weight

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        dim_ens : int
                Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        gamma : double
                Hybrid weight provided by PDAF
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, dim_ens)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, dim_ens)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, dim_ens)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    py__likelihood_l_pdaf : Callable
        Compute likelihood

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    py__likelihood_hyb_l_pdaf : Callable
        Compute likelihood with hybrid weight

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                Input vector holding the local residual
                Array shape: (dim_obs_l)
        gamma : double
                Hybrid weight provided by PDAF

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                Input vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaflknetf_update_step(&step, &dim_p, &dim_obs_f, &dim_ens,
                                  &state_p[0], &ainv[0,0], &ens_p[0,0],
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__init_obs_pdaf,
                                  pdaf_cb.c__init_obs_l_pdaf,
                                  pdaf_cb.c__prodrinva_hyb_l_pdaf,
                                  pdaf_cb.c__init_n_domains_p_pdaf,
                                  pdaf_cb.c__init_dim_l_pdaf,
                                  pdaf_cb.c__init_dim_obs_l_pdaf,
                                  pdaf_cb.c__g2l_state_pdaf,
                                  pdaf_cb.c__l2g_state_pdaf,
                                  pdaf_cb.c__g2l_obs_pdaf,
                                  pdaf_cb.c__init_obsvar_pdaf,
                                  pdaf_cb.c__init_obsvar_l_pdaf,
                                  pdaf_cb.c__likelihood_l_pdaf,
                                  pdaf_cb.c__likelihood_hyb_l_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf, &screen,
                                  &subtype, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, flag


def letkf_update(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__prepoststep_pdaf, int  screen,
    int  subtype, int  dim_lag, double [::1,:,:] sens_p, int  cnt_maxlag,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from full state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdafletkf_update(&step, &dim_p, &dim_obs_f, &dim_ens,
                            &state_p[0], &ainv[0,0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__init_obs_l_pdaf,
                            pdaf_cb.c__prodrinva_l_pdaf,
                            pdaf_cb.c__init_n_domains_p_pdaf,
                            pdaf_cb.c__init_dim_l_pdaf,
                            pdaf_cb.c__init_dim_obs_l_pdaf,
                            pdaf_cb.c__g2l_state_pdaf,
                            pdaf_cb.c__l2g_state_pdaf,
                            pdaf_cb.c__g2l_obs_pdaf,
                            pdaf_cb.c__init_obsvar_pdaf,
                            pdaf_cb.c__init_obsvar_l_pdaf,
                            pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                            &dim_lag, &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def lseik_ana_trans(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, int  rank, double [::1] state_l, double [::1,:] uinv_l,
    double [::1,:] ens_l, double [::1,:] hl_l, double [::1] hxbar_l,
    double [::1] obs_l, double [::1,:] omegat_in, double  forget,
    py__prodrinva_l_pdaf, int  nm1vsn, int  type_sqrt, int  screen,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_l : ndarray[np.float64, ndim=1]
        on exit: state on local analysis domain
        Array shape: (dim_l)
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    omegat_in : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    nm1vsn : int
        Whether covariance is normalized with 1/N or 1/(N-1)
    type_sqrt : int
        Type of square-root of A
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        on exit: state on local analysis domain
        Array shape: (dim_l)
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    omegat_in : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_l_np = np.asarray(uinv_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_l_np = np.asarray(hl_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] omegat_in_np = np.asarray(omegat_in, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_lseik_ana_trans(&domain_p, &step, &dim_l, &dim_obs_l,
                                &dim_ens, &rank, &state_l[0], &uinv_l[0,0],
                                &ens_l[0,0], &hl_l[0,0], &hxbar_l[0],
                                &obs_l[0], &omegat_in[0,0], &forget,
                                pdaf_cb.c__prodrinva_l_pdaf, &nm1vsn,
                                &type_sqrt, &screen, &debug, &flag)

    return state_l_np, uinv_l_np, ens_l_np, hl_l_np, omegat_in_np, forget, flag


def en3dvar_optim_lbfgs(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec_p, int  dim_obs_p, double [::1,:] ens_p,
    double [::1] obs_p, double [::1] dy_p, double [::1] v_p,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_en3dvar_optim_lbfgs(&step, &dim_p, &dim_ens, &dim_cvec_p,
                                    &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                    &dy_p[0], &v_p[0],
                                    pdaf_cb.c__prodrinva_pdaf,
                                    pdaf_cb.c__cvt_ens_pdaf,
                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf,
                                    &opt_parallel, &screen)

    return v_p_np


def en3dvar_optim_cgplus(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec_p, int  dim_obs_p, double [::1,:] ens_p,
    double [::1] obs_p, double [::1] dy_p, double [::1] v_p,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_en3dvar_optim_cgplus(&step, &dim_p, &dim_ens, &dim_cvec_p,
                                     &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                     &dy_p[0], &v_p[0],
                                     pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     &opt_parallel, &screen)

    return v_p_np


def en3dvar_optim_cg(int  step, int  dim_p, int  dim_ens, int  dim_cvec_p,
    int  dim_obs_p, double [::1,:] ens_p, double [::1] obs_p,
    double [::1] dy_p, double [::1] v_p, py__prodrinva_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_en3dvar_optim_cg(&step, &dim_p, &dim_ens, &dim_cvec_p,
                                 &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                 &dy_p[0], &v_p[0],
                                 pdaf_cb.c__prodrinva_pdaf,
                                 pdaf_cb.c__cvt_ens_pdaf,
                                 pdaf_cb.c__cvt_adj_ens_pdaf,
                                 pdaf_cb.c__obs_op_lin_pdaf,
                                 pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel,
                                 &screen)

    return v_p_np


def en3dvar_costf_cvt(int  step, int  iter, int  dim_p, int  dim_ens,
    int  dim_cvec_p, int  dim_obs_p, double [::1,:] ens_p,
    double [::1] obs_p, double [::1] dy_p, double [::1] v_p,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        Optimization iteration
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cvec_p : int
        PE-local size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector

    Returns
    -------
    j_tot : double
        on exit: Value of cost function
    gradj : ndarray[np.float64, ndim=1]
        on exit: PE-local gradient of J
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] gradj = gradj_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_en3dvar_costf_cvt(&step, &iter, &dim_p, &dim_ens,
                                  &dim_cvec_p, &dim_obs_p, &ens_p[0,0],
                                  &obs_p[0], &dy_p[0], &v_p[0], &j_tot,
                                  &gradj[0], pdaf_cb.c__prodrinva_pdaf,
                                  pdaf_cb.c__cvt_ens_pdaf,
                                  pdaf_cb.c__cvt_adj_ens_pdaf,
                                  pdaf_cb.c__obs_op_lin_pdaf,
                                  pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel)

    return j_tot, gradj_np


def en3dvar_costf_cg_cvt(int  step, int  iter, int  dim_p, int  dim_ens,
    int  dim_cvec_p, int  dim_obs_p, double [::1,:] ens_p,
    double [::1] obs_p, double [::1] dy_p, double [::1] v_p,
    double [::1] d_p, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    int  opt_parallel):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        Optimization iteration
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cvec_p : int
        PE-local size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    d_p : ndarray[np.float64, ndim=1]
        CG descent direction
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector

    Returns
    -------
    d_p : ndarray[np.float64, ndim=1]
        CG descent direction
        Array shape: (dim_cvec_p)
    j_tot : double
        on exit: Value of cost function
    gradj : ndarray[np.float64, ndim=1]
        on exit: gradient of J
        Array shape: (dim_cvec_p)
    hessjd : ndarray[np.float64, ndim=1]
        on exit: Hessian of J times d_p
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] d_p_np = np.asarray(d_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] gradj = gradj_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hessjd_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] hessjd = hessjd_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_en3dvar_costf_cg_cvt(&step, &iter, &dim_p, &dim_ens,
                                     &dim_cvec_p, &dim_obs_p, &ens_p[0,0],
                                     &obs_p[0], &dy_p[0], &v_p[0], &d_p[0],
                                     &j_tot, &gradj[0], &hessjd[0],
                                     pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel)

    return d_p_np, j_tot, gradj_np, hessjd_np


def gather_ens(int  dim_p, int  dim_ens_p, double [::1,:] ens,
    double [::1] state, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens_p : int
        Size of ensemble
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    screen : int
        Verbosity flag

    Returns
    -------
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 ens_cfi
    cdef CFI_cdesc_t *ens_ptr = <CFI_cdesc_t *> &ens_cfi
    cdef size_t ens_nbytes = ens.nbytes
    cdef CFI_index_t ens_extent[2]
    ens_extent[0] = ens.shape[0]
    ens_extent[1] = ens.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.asarray(ens, dtype=np.float64, order="F")

    cdef CFI_cdesc_rank1 state_cfi
    cdef CFI_cdesc_t *state_ptr = <CFI_cdesc_t *> &state_cfi
    cdef size_t state_nbytes = state.nbytes
    cdef CFI_index_t state_extent[1]
    state_extent[0] = state.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")

    with nogil:
        CFI_establish(ens_ptr, &ens[0,0], CFI_attribute_other,
                      CFI_type_double , ens_nbytes, 2, ens_extent)

        CFI_establish(state_ptr, &state[0], CFI_attribute_other,
                      CFI_type_double , state_nbytes, 1, state_extent)

        c__pdaf_gather_ens(&dim_p, &dim_ens_p, ens_ptr, state_ptr, &screen)

    return ens_np, state_np


def scatter_ens(int  dim_p, int  dim_ens_p, double [::1,:] ens,
    double [::1] state, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens_p : int
        Size of ensemble
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    state : ndarray[np.float64, ndim=1]
        PE-local state vector (for SEEK)
        Array shape: (:)
    screen : int
        Verbosity flag

    Returns
    -------
    ens : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (:, :)
    state : ndarray[np.float64, ndim=1]
        PE-local state vector (for SEEK)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank2 ens_cfi
    cdef CFI_cdesc_t *ens_ptr = <CFI_cdesc_t *> &ens_cfi
    cdef size_t ens_nbytes = ens.nbytes
    cdef CFI_index_t ens_extent[2]
    ens_extent[0] = ens.shape[0]
    ens_extent[1] = ens.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.asarray(ens, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_cfi
    cdef CFI_cdesc_t *state_ptr = <CFI_cdesc_t *> &state_cfi
    cdef size_t state_nbytes = state.nbytes
    cdef CFI_index_t state_extent[1]
    state_extent[0] = state.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_np = np.asarray(state, dtype=np.float64, order="F")
    with nogil:
        CFI_establish(ens_ptr, &ens[0,0], CFI_attribute_other,
                      CFI_type_double , ens_nbytes, 2, ens_extent)

        CFI_establish(state_ptr, &state[0], CFI_attribute_other,
                      CFI_type_double , state_nbytes, 1, state_extent)

        c__pdaf_scatter_ens(&dim_p, &dim_ens_p, ens_ptr, state_ptr, &screen)

    return ens_np, state_np


def hyb3dvar_optim_lbfgs(int  step, int  dim_p, int  dim_ens,
    int  dim_cv_par_p, int  dim_cv_ens_p, int  dim_obs_p,
    double [::1,:] ens_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_par_p, double [::1] v_ens_p, py__prodrinva_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel,
    double  beta_3dvar, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cv_par_p : int
        Size of control vector (parameterized)
    dim_cv_ens_p : int
        Size of control vector (ensemble)
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    beta_3dvar : double
        Hybrid weight
    screen : int
        Verbosity flag

    Returns
    -------
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_par_p_np = np.asarray(v_par_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_ens_p_np = np.asarray(v_ens_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_hyb3dvar_optim_lbfgs(&step, &dim_p, &dim_ens,
                                     &dim_cv_par_p, &dim_cv_ens_p,
                                     &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                     &dy_p[0], &v_par_p[0], &v_ens_p[0],
                                     pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     &opt_parallel, &beta_3dvar, &screen)

    return v_par_p_np, v_ens_p_np


def hyb3dvar_optim_cgplus(int  step, int  dim_p, int  dim_ens,
    int  dim_cv_par_p, int  dim_cv_ens_p, int  dim_obs_p,
    double [::1,:] ens_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_par_p, double [::1] v_ens_p, py__prodrinva_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel,
    double  beta_3dvar, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cv_par_p : int
        Size of control vector (parameterized)
    dim_cv_ens_p : int
        Size of control vector (ensemble)
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    beta_3dvar : double
        Hybrid weight
    screen : int
        Verbosity flag

    Returns
    -------
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_par_p_np = np.asarray(v_par_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_ens_p_np = np.asarray(v_ens_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_hyb3dvar_optim_cgplus(&step, &dim_p, &dim_ens,
                                      &dim_cv_par_p, &dim_cv_ens_p,
                                      &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                      &dy_p[0], &v_par_p[0], &v_ens_p[0],
                                      pdaf_cb.c__prodrinva_pdaf,
                                      pdaf_cb.c__cvt_pdaf,
                                      pdaf_cb.c__cvt_adj_pdaf,
                                      pdaf_cb.c__cvt_ens_pdaf,
                                      pdaf_cb.c__cvt_adj_ens_pdaf,
                                      pdaf_cb.c__obs_op_lin_pdaf,
                                      pdaf_cb.c__obs_op_adj_pdaf,
                                      &opt_parallel, &beta_3dvar, &screen)

    return v_par_p_np, v_ens_p_np


def hyb3dvar_optim_cg(int  step, int  dim_p, int  dim_ens,
    int  dim_cv_par_p, int  dim_cv_ens_p, int  dim_obs_p,
    double [::1,:] ens_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_par_p, double [::1] v_ens_p, py__prodrinva_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel,
    double  beta_3dvar, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cv_par_p : int
        Size of control vector (parameterized)
    dim_cv_ens_p : int
        Size of control vector (ensemble)
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    beta_3dvar : double
        Hybrid weight
    screen : int
        Verbosity flag

    Returns
    -------
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_par_p_np = np.asarray(v_par_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_ens_p_np = np.asarray(v_ens_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_hyb3dvar_optim_cg(&step, &dim_p, &dim_ens, &dim_cv_par_p,
                                  &dim_cv_ens_p, &dim_obs_p, &ens_p[0,0],
                                  &obs_p[0], &dy_p[0], &v_par_p[0],
                                  &v_ens_p[0], pdaf_cb.c__prodrinva_pdaf,
                                  pdaf_cb.c__cvt_pdaf,
                                  pdaf_cb.c__cvt_adj_pdaf,
                                  pdaf_cb.c__cvt_ens_pdaf,
                                  pdaf_cb.c__cvt_adj_ens_pdaf,
                                  pdaf_cb.c__obs_op_lin_pdaf,
                                  pdaf_cb.c__obs_op_adj_pdaf,
                                  &opt_parallel, &beta_3dvar, &screen)

    return v_par_p_np, v_ens_p_np


def hyb3dvar_costf_cvt(int  step, int  iter, int  dim_p, int  dim_ens,
    int  dim_cv_p, int  dim_cv_par_p, int  dim_cv_ens_p, int  dim_obs_p,
    double [::1,:] ens_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_par_p, double [::1] v_ens_p, double [::1] v_p,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    int  opt_parallel, double  beta):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        Optimization iteration
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cv_p : int
        Size of control vector (full)
    dim_cv_par_p : int
        Size of control vector (parameterized part)
    dim_cv_ens_p : int
        Size of control vector (ensemble part)
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        background innovation
        Array shape: (dim_obs_p)
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector (full)
        Array shape: (dim_cv_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    beta : double
        Hybrid weight

    Returns
    -------
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    j_tot : double
        on exit: Value of cost function
    gradj : ndarray[np.float64, ndim=1]
        on exit: PE-local gradient of J (full)
        Array shape: (dim_cv_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_par_p_np = np.asarray(v_par_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_ens_p_np = np.asarray(v_ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_np = np.zeros((dim_cv_p), dtype=np.float64, order="F")
    cdef double [::1] gradj = gradj_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_hyb3dvar_costf_cvt(&step, &iter, &dim_p, &dim_ens,
                                   &dim_cv_p, &dim_cv_par_p, &dim_cv_ens_p,
                                   &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                   &dy_p[0], &v_par_p[0], &v_ens_p[0],
                                   &v_p[0], &j_tot, &gradj[0],
                                   pdaf_cb.c__prodrinva_pdaf,
                                   pdaf_cb.c__cvt_pdaf,
                                   pdaf_cb.c__cvt_adj_pdaf,
                                   pdaf_cb.c__cvt_ens_pdaf,
                                   pdaf_cb.c__cvt_adj_ens_pdaf,
                                   pdaf_cb.c__obs_op_lin_pdaf,
                                   pdaf_cb.c__obs_op_adj_pdaf,
                                   &opt_parallel, &beta)

    return v_par_p_np, v_ens_p_np, j_tot, gradj_np


def hyb3dvar_costf_cg_cvt(int  step, int  iter, int  dim_p, int  dim_ens,
    int  dim_cv_par_p, int  dim_cv_ens_p, int  dim_obs_p,
    double [::1,:] ens_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_par_p, double [::1] v_ens_p, double [::1] d_par_p,
    double [::1] d_ens_p, py__prodrinva_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, double  beta):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        Optimization iteration
    dim_p : int
        PE-local state dimension
    dim_ens : int
        ensemble size
    dim_cv_par_p : int
        Size of control vector (parameterized part)
    dim_cv_ens_p : int
        Size of control vector (ensemble part)
    dim_obs_p : int
        PE-local dimension of observation vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_par_p : ndarray[np.float64, ndim=1]
        Control vector (parameterized part)
        Array shape: (dim_cv_par_p)
    v_ens_p : ndarray[np.float64, ndim=1]
        Control vector (ensemble part)
        Array shape: (dim_cv_ens_p)
    d_par_p : ndarray[np.float64, ndim=1]
        CG descent direction (parameterized part)
        Array shape: (dim_cv_par_p)
    d_ens_p : ndarray[np.float64, ndim=1]
        CG descent direction (ensemble part)
        Array shape: (dim_cv_ens_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    beta : double
        Hybrid weight

    Returns
    -------
    d_par_p : ndarray[np.float64, ndim=1]
        CG descent direction (parameterized part)
        Array shape: (dim_cv_par_p)
    d_ens_p : ndarray[np.float64, ndim=1]
        CG descent direction (ensemble part)
        Array shape: (dim_cv_ens_p)
    j_tot : double
        on exit: Value of cost function
    gradj_par : ndarray[np.float64, ndim=1]
        on exit: gradient of J (parameterized part)
        Array shape: (dim_cv_par_p)
    gradj_ens : ndarray[np.float64, ndim=1]
        on exit: gradient of J (ensemble part)
        Array shape: (dim_cv_ens_p)
    hessjd_par : ndarray[np.float64, ndim=1]
        on exit: Hessian of J times d_p (parameterized part)
        Array shape: (dim_cv_par_p)
    hessjd_ens : ndarray[np.float64, ndim=1]
        on exit: Hessian of J times d_p (ensemble part)
        Array shape: (dim_cv_ens_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] d_par_p_np = np.asarray(d_par_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] d_ens_p_np = np.asarray(d_ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_par_np = np.zeros((dim_cv_par_p), dtype=np.float64, order="F")
    cdef double [::1] gradj_par = gradj_par_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_ens_np = np.zeros((dim_cv_ens_p), dtype=np.float64, order="F")
    cdef double [::1] gradj_ens = gradj_ens_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hessjd_par_np = np.zeros((dim_cv_par_p), dtype=np.float64, order="F")
    cdef double [::1] hessjd_par = hessjd_par_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hessjd_ens_np = np.zeros((dim_cv_ens_p), dtype=np.float64, order="F")
    cdef double [::1] hessjd_ens = hessjd_ens_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_hyb3dvar_costf_cg_cvt(&step, &iter, &dim_p, &dim_ens,
                                      &dim_cv_par_p, &dim_cv_ens_p,
                                      &dim_obs_p, &ens_p[0,0], &obs_p[0],
                                      &dy_p[0], &v_par_p[0], &v_ens_p[0],
                                      &d_par_p[0], &d_ens_p[0], &j_tot,
                                      &gradj_par[0], &gradj_ens[0],
                                      &hessjd_par[0], &hessjd_ens[0],
                                      pdaf_cb.c__prodrinva_pdaf,
                                      pdaf_cb.c__cvt_pdaf,
                                      pdaf_cb.c__cvt_adj_pdaf,
                                      pdaf_cb.c__cvt_ens_pdaf,
                                      pdaf_cb.c__cvt_adj_ens_pdaf,
                                      pdaf_cb.c__obs_op_lin_pdaf,
                                      pdaf_cb.c__obs_op_adj_pdaf,
                                      &opt_parallel, &beta)

    return d_par_p_np, d_ens_p_np, j_tot, gradj_par_np, gradj_ens_np, hessjd_par_np, hessjd_ens_np


def print_version():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_print_version()



def en3dvar_analysis_cvt(int  step, int  dim_p, int  dim_obs_p,
    int  dim_ens, int  dim_cvec_ens, double [::1,:] ens_p,
    double [::1] state_inc_p, double [::1] hxbar_p, double [::1] obs_p,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  screen, int  type_opt,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    dim_cvec_ens : int
        Size of control vector
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    state_inc_p : ndarray[np.float64, ndim=1]
        PE-local state analysis increment
        Array shape: (dim_p)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    screen : int
        Verbosity flag
    type_opt : int
        Type of minimizer for 3DVar
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    state_inc_p : ndarray[np.float64, ndim=1]
        PE-local state analysis increment
        Array shape: (dim_p)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_inc_p_np = np.asarray(state_inc_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdafen3dvar_analysis_cvt(&step, &dim_p, &dim_obs_p, &dim_ens,
                                    &dim_cvec_ens, &state_p[0],
                                    &ens_p[0,0], &state_inc_p[0],
                                    &hxbar_p[0], &obs_p[0],
                                    pdaf_cb.c__prodrinva_pdaf,
                                    pdaf_cb.c__cvt_ens_pdaf,
                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf, &screen,
                                    &type_opt, &debug, &flag)

    return state_p_np, ens_p_np, state_inc_p_np, flag


def sisort(int  n, double [::1] veca):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    n : int

    veca : ndarray[np.float64, ndim=1]

        Array shape: (n)

    Returns
    -------
    veca : ndarray[np.float64, ndim=1]

        Array shape: (n)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] veca_np = np.asarray(veca, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_sisort(&n, &veca[0])

    return veca_np

def enkf_ana_rlm(int  step, int  dim_p, int  dim_obs_p, int dim_obs, int  dim_ens,
    int  rank_ana, double [::1] state_p, double [::1,:] ens_p,
    double [::1,:] hzb, double [::1,:] hx_p, double [::1] hxbar_p,
    double [::1] obs_p, py__add_obs_err_pdaf, py__init_obs_covar_pdaf,
    int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_obs : int
        Global dimension of observation vector
    dim_ens : int
        Size of state ensemble
    rank_ana : int
        Rank to be considered for inversion of HPH
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hzb : ndarray[np.float64, ndim=2]
        Ensemble tranformation matrix
        Array shape: (dim_ens, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Dimension of observation vector
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hzb : ndarray[np.float64, ndim=2]
        Ensemble tranformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hzb_np = np.asarray(hzb, dtype=np.float64, order="F")
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    with nogil:
        c__pdaf_enkf_ana_rlm(&step, &dim_p, &dim_obs_p, &dim_obs, &dim_ens,
                             &rank_ana, &state_p[0], &ens_p[0,0],
                             &hzb[0,0], &hx_p[0,0], &hxbar_p[0], &obs_p[0],
                             pdaf_cb.c__add_obs_err_pdaf,
                             pdaf_cb.c__init_obs_covar_pdaf, &screen,
                             &debug, &flag)

    return state_p_np, ens_p_np, hzb_np, flag


def smoother_enkf(int  dim_p, int  dim_ens, int  dim_lag,
    double [::1,:] ainv, double [::1,:,:] sens_p, int  cnt_maxlag,
    double  forget, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ainv : ndarray[np.float64, ndim=2]
        Weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    forget : double
        Forgetting factor
    screen : int
        Verbosity flag

    Returns
    -------
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_smoother_enkf(&dim_p, &dim_ens, &dim_lag, &ainv[0,0],
                              &sens_p[0,0,0], &cnt_maxlag, &forget, &screen)

    return sens_p_np, cnt_maxlag


def ensrf_update(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__init_obsvar_pdaf,
    py__localize_covar_serial_pdaf, py__prepoststep_pdaf, int  screen,
    int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of state ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obsvar_pdaf : Callable
        Initialize vector of observation error variances

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__localize_covar_serial_pdaf : Callable
        Apply localization for single-observation vectors

        Callback Parameters
        -------------------
        iobs : int
                Index of current observation
        dim_p : int
                Process-local state dimension
        dim_obs : int
                Number of observations
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

        Callback Returns
        ----------------
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

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

    screen : int
        Verbosity flag
    subtype : int
        Specification of filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafensrf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                            &state_p[0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__init_obsvar_pdaf,
                            pdaf_cb.c__localize_covar_serial_pdaf,
                            pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                            &flag)

    return dim_obs_p, state_p_np, ens_p_np, flag


def pf_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ens_p, int  type_resample,
    int  type_winf, double  limit_winf, int  type_noise, double  noise_amp,
    double [::1,:] hz_p, double [::1] obs_p, py__likelihood_pdaf,
    int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local forecast mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    type_resample : int
        Type of resampling scheme
    type_winf : int
        Type of weights inflation
    limit_winf : double
        Limit for weights inflation
    type_noise : int
        Type of pertubing noise
    noise_amp : double
        Amplitude of noise
    hz_p : ndarray[np.float64, ndim=2]
        Temporary matrices for analysis
        Array shape: (dim_obs_p, dim_ens)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__likelihood_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        resid : ndarray[np.float64, ndim=1]
                Input vector holding the residual
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        likely : double
                Output value of the likelihood

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local forecast mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    with nogil:
        c__pdaf_pf_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &state_p[0],
                       &ens_p[0,0], &type_resample, &type_winf,
                       &limit_winf, &type_noise, &noise_amp, &hz_p[0,0],
                       &obs_p[0], pdaf_cb.c__likelihood_pdaf, &screen,
                       &debug, &flag)

    return state_p_np, ens_p_np, flag


def pf_resampling(int  method, int  nin, int  nout, double [::1] weights,
    int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    method : int
        Choose resampling method
    nin : int
        number of particles
    nout : int
        number of particles to be resampled
    weights : ndarray[np.float64, ndim=1]
        Weights
        Array shape: (nin)
    screen : int
        Verbosity flag

    Returns
    -------
    ids : ndarray[np.intc, ndim=1]
        Indices of resampled ensmeble states
        Array shape: (nout)
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] ids_np = np.zeros((nout), dtype=np.intc, order="F")
    cdef int [::1] ids = ids_np
    with nogil:
        c__pdaf_pf_resampling(&method, &nin, &nout, &weights[0], &ids[0],
                              &screen)

    return ids_np


def mvnormalize(int  mode, int  dim_state, int  dim_field, int  offset,
    int  ncol, double [::1,:] states, double  stddev):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    mode : int
        Mode: (1) normalize, (2) re-scale
    dim_state : int
        Dimension of state vector
    dim_field : int
        Dimension of a field in state vector
    offset : int
        Offset of field in state vector
    ncol : int
        Number of columns in array states
    states : ndarray[np.float64, ndim=2]
        State vector array
        Array shape: (dim_state, ncol)
    stddev : double
        Standard deviation of field

    Returns
    -------
    states : ndarray[np.float64, ndim=2]
        State vector array
        Array shape: (dim_state, ncol)
    stddev : double
        Standard deviation of field
    status : int
        Status flag (0=success)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] states_np = np.asarray(states, dtype=np.float64, order="F")
    cdef int  status
    with nogil:
        c__pdaf_mvnormalize(&mode, &dim_state, &dim_field, &offset, &ncol,
                            &states[0,0], &stddev, &status)

    return states_np, stddev, status


def _3dvar_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_3dvar_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def _3dvar_alloc(int  subtype, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_3dvar_alloc(&subtype, &outflag)

    return outflag


def _3dvar_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_3dvar_config(&subtype, &verbose)

    return subtype


def _3dvar_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_3dvar_set_iparam(&id, &value, &flag)

    return flag


def _3dvar_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_3dvar_set_rparam(&id, &value, &flag)

    return flag


def _3dvar_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_3dvar_options()



def _3dvar_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_3dvar_memtime(&printtype)



def reset_dim_ens(int  dim_ens_in, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_ens_in : int
        Sub-type of filter
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_reset_dim_ens(&dim_ens_in, &outflag)

    return outflag


def reset_dim_p(int  dim_p_in, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p_in : int
        Sub-type of filter
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_reset_dim_p(&dim_p_in, &outflag)

    return outflag


def _3dvar_optim_lbfgs(int  step, int  dim_p, int  dim_cvec_p,
    int  dim_obs_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_p, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_3dvar_optim_lbfgs(&step, &dim_p, &dim_cvec_p, &dim_obs_p,
                                  &obs_p[0], &dy_p[0], &v_p[0],
                                  pdaf_cb.c__prodrinva_pdaf,
                                  pdaf_cb.c__cvt_pdaf,
                                  pdaf_cb.c__cvt_adj_pdaf,
                                  pdaf_cb.c__obs_op_lin_pdaf,
                                  pdaf_cb.c__obs_op_adj_pdaf,
                                  &opt_parallel, &screen)

    return v_p_np


def _3dvar_optim_cgplus(int  step, int  dim_p, int  dim_cvec_p,
    int  dim_obs_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_p, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_3dvar_optim_cgplus(&step, &dim_p, &dim_cvec_p, &dim_obs_p,
                                   &obs_p[0], &dy_p[0], &v_p[0],
                                   pdaf_cb.c__prodrinva_pdaf,
                                   pdaf_cb.c__cvt_pdaf,
                                   pdaf_cb.c__cvt_adj_pdaf,
                                   pdaf_cb.c__obs_op_lin_pdaf,
                                   pdaf_cb.c__obs_op_adj_pdaf,
                                   &opt_parallel, &screen)

    return v_p_np


def _3dvar_optim_cg(int  step, int  dim_p, int  dim_cvec_p, int  dim_obs_p,
    double [::1] obs_p, double [::1] dy_p, double [::1] v_p,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local state dimension
    dim_cvec_p : int
        Size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector
    screen : int
        Verbosity flag

    Returns
    -------
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] v_p_np = np.asarray(v_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf_3dvar_optim_cg(&step, &dim_p, &dim_cvec_p, &dim_obs_p,
                               &obs_p[0], &dy_p[0], &v_p[0],
                               pdaf_cb.c__prodrinva_pdaf,
                               pdaf_cb.c__cvt_pdaf,
                               pdaf_cb.c__cvt_adj_pdaf,
                               pdaf_cb.c__obs_op_lin_pdaf,
                               pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel,
                               &screen)

    return v_p_np


def _3dvar_costf_cvt(int  step, int  iter, int  dim_p, int  dim_cvec_p,
    int  dim_obs_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_p, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  opt_parallel):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        Optimization iteration
    dim_p : int
        PE-local state dimension
    dim_cvec_p : int
        PE-local size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        control vector
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector

    Returns
    -------
    j_tot : double
        on exit: Value of cost function
    gradj : ndarray[np.float64, ndim=1]
        on exit: PE-local gradient of J
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] gradj = gradj_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_3dvar_costf_cvt(&step, &iter, &dim_p, &dim_cvec_p,
                                &dim_obs_p, &obs_p[0], &dy_p[0], &v_p[0],
                                &j_tot, &gradj[0],
                                pdaf_cb.c__prodrinva_pdaf,
                                pdaf_cb.c__cvt_pdaf,
                                pdaf_cb.c__cvt_adj_pdaf,
                                pdaf_cb.c__obs_op_lin_pdaf,
                                pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel)

    return j_tot, gradj_np


def _3dvar_costf_cg_cvt(int  step, int  iter, int  dim_p, int  dim_cvec_p,
    int  dim_obs_p, double [::1] obs_p, double [::1] dy_p,
    double [::1] v_p, double [::1] d_p, py__prodrinva_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    int  opt_parallel):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    iter : int
        CG iteration
    dim_p : int
        PE-local state dimension
    dim_cvec_p : int
        PE-local size of control vector
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        Vector of observations
        Array shape: (dim_obs_p)
    dy_p : ndarray[np.float64, ndim=1]
        Background innovation
        Array shape: (dim_obs_p)
    v_p : ndarray[np.float64, ndim=1]
        Control vector
        Array shape: (dim_cvec_p)
    d_p : ndarray[np.float64, ndim=1]
        CG descent direction
        Array shape: (dim_cvec_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    opt_parallel : int
        Whether to use a decomposed control vector

    Returns
    -------
    d_p : ndarray[np.float64, ndim=1]
        CG descent direction
        Array shape: (dim_cvec_p)
    j_tot : double
        on exit: Value of cost function
    gradj : ndarray[np.float64, ndim=1]
        on exit: gradient of J
        Array shape: (dim_cvec_p)
    hessjd : ndarray[np.float64, ndim=1]
        on exit: Hessian of J times d_p
        Array shape: (dim_cvec_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] d_p_np = np.asarray(d_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gradj_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] gradj = gradj_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hessjd_np = np.zeros((dim_cvec_p), dtype=np.float64, order="F")
    cdef double [::1] hessjd = hessjd_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    cdef double  j_tot
    with nogil:
        c__pdaf_3dvar_costf_cg_cvt(&step, &iter, &dim_p, &dim_cvec_p,
                                   &dim_obs_p, &obs_p[0], &dy_p[0],
                                   &v_p[0], &d_p[0], &j_tot, &gradj[0],
                                   &hessjd[0], pdaf_cb.c__prodrinva_pdaf,
                                   pdaf_cb.c__cvt_pdaf,
                                   pdaf_cb.c__cvt_adj_pdaf,
                                   pdaf_cb.c__obs_op_lin_pdaf,
                                   pdaf_cb.c__obs_op_adj_pdaf, &opt_parallel)

    return d_p_np, j_tot, gradj_np, hessjd_np


def lknetf_analysis_t(int  domain_p, int  step, int  dim_l,
    int  dim_obs_l, int  dim_ens, double [::1] state_l,
    double [::1,:] ens_l, double [::1,:] hx_l, double [::1] hxbar_l,
    double [::1] obs_l, double [::1,:] rndmat, double  forget,
    py__prodrinva_l_pdaf, py__init_obsvar_l_pdaf, py__likelihood_l_pdaf,
    int  screen, int  type_forget, double [::1] eff_dimens, int  type_hyb,
    double  hyb_g, double  hyb_k, double [::1] gamma,
    double [::1] skew_mabs, double [::1] kurt_mabs, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    state_l : ndarray[np.float64, ndim=1]
        local forecast state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hx_l : ndarray[np.float64, ndim=2]
        local observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        local observed ens. mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    py__likelihood_l_pdaf : Callable
        Provide likelihood of an ensemble state

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    screen : int
        Verbosity flag
    type_forget : int
        Type of forgetting factor
    eff_dimens : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    type_hyb : int
        Type of hybrid weight
    hyb_g : double
        Prescribed hybrid weight for state transformation
    hyb_k : double
        Scale factor kappa (for type_hyb 3 and 4)
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    skew_mabs : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    kurt_mabs : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        local forecast state
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        on exit: local weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    eff_dimens : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    skew_mabs : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    kurt_mabs : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv_l = ainv_l_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] rndmat_np = np.asarray(rndmat, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] eff_dimens_np = np.asarray(eff_dimens, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gamma_np = np.asarray(gamma, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] skew_mabs_np = np.asarray(skew_mabs, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] kurt_mabs_np = np.asarray(kurt_mabs, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    with nogil:
        c__pdaf_lknetf_analysis_t(&domain_p, &step, &dim_l, &dim_obs_l,
                                  &dim_ens, &state_l[0], &ainv_l[0,0],
                                  &ens_l[0,0], &hx_l[0,0], &hxbar_l[0],
                                  &obs_l[0], &rndmat[0,0], &forget,
                                  pdaf_cb.c__prodrinva_l_pdaf,
                                  pdaf_cb.c__init_obsvar_l_pdaf,
                                  pdaf_cb.c__likelihood_l_pdaf, &screen,
                                  &type_forget, &eff_dimens[0], &type_hyb,
                                  &hyb_g, &hyb_k, &gamma[0], &skew_mabs[0],
                                  &kurt_mabs[0], &flag)

    return state_l_np, ainv_l_np, ens_l_np, rndmat_np, forget, eff_dimens_np, gamma_np, skew_mabs_np, kurt_mabs_np, flag


def get_ensstats():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    cdef CFI_cdesc_rank1 skew_ptr_cfi
    cdef CFI_cdesc_t *skew_ptr_ptr = <CFI_cdesc_t *> &skew_ptr_cfi
    cdef CFI_cdesc_rank1 kurt_ptr_cfi
    cdef CFI_cdesc_t *kurt_ptr_ptr = <CFI_cdesc_t *> &kurt_ptr_cfi
    cdef int  status
    with nogil:
        c__pdaf_get_ensstats(skew_ptr_ptr, kurt_ptr_ptr, &status)

    cdef CFI_index_t skew_ptr_subscripts[1]
    skew_ptr_subscripts[0] = 0
    cdef double * skew_ptr_ptr_np
    skew_ptr_ptr_np = <double *>CFI_address(skew_ptr_ptr, skew_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] skew_ptr_np = np.asarray(<double [:skew_ptr_ptr.dim[0].extent:1]> skew_ptr_ptr_np, order="F")
    cdef CFI_index_t kurt_ptr_subscripts[1]
    kurt_ptr_subscripts[0] = 0
    cdef double * kurt_ptr_ptr_np
    kurt_ptr_ptr_np = <double *>CFI_address(kurt_ptr_ptr, kurt_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] kurt_ptr_np = np.asarray(<double [:kurt_ptr_ptr.dim[0].extent:1]> kurt_ptr_ptr_np, order="F")
    return skew_ptr_np, kurt_ptr_np, status


def estkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_estkf_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def estkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_estkf_alloc(&outflag)

    return outflag


def estkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_estkf_config(&subtype, &verbose)

    return subtype


def estkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_estkf_set_iparam(&id, &value, &flag)

    return flag


def estkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_estkf_set_rparam(&id, &value, &flag)

    return flag


def estkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_estkf_options()



def estkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_estkf_memtime(&printtype)



def gen_obs(int  step, int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_f_pdaf,
    py__obs_op_f_pdaf, py__get_obs_f_pdaf, py__init_obserr_f_pdaf,
    py__prepoststep_pdaf, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_f_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_f_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__get_obs_f_pdaf : Callable
        Provide observation vector to user

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_f : int
                Size of the full observation vector

        Callback Returns
        ----------------
        observation_f : ndarray[np.float64, ndim=1]
                Full vector of synthetic observations (process-local)
                Array shape: (dim_obs_f)

    py__init_obserr_f_pdaf : Callable
        Initialize vector of observation error standard deviations

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_f : int
                Full dimension of observation vector
        obs_f : ndarray[np.float64, ndim=1]
                Full observation vector
                Array shape: (dim_obs_f)

        Callback Returns
        ----------------
        obserr_f : ndarray[np.float64, ndim=1]
                Full observation error stddev
                Array shape: (dim_obs_f)

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

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    pdaf_cb.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaf_gen_obs(&step, &dim_p, &dim_obs_f, &dim_ens, &state_p[0],
                        &ainv[0,0], &ens_p[0,0],
                        pdaf_cb.c__init_dim_obs_f_pdaf,
                        pdaf_cb.c__obs_op_f_pdaf,
                        pdaf_cb.c__get_obs_f_pdaf,
                        pdaf_cb.c__init_obserr_f_pdaf,
                        pdaf_cb.c__prepoststep_pdaf, &screen, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, flag


def obs_init(int  step, int  dim_p, int  dim_ens, int  dim_obs_p,
    double [::1] state_p, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, int  screen, int  debug,
    bint  do_ens_mean, bint  do_init_dim, bint  do_hx, bint  do_hxbar,
    bint  do_init_obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    do_ens_mean : bint
        Whether to compute ensemble mean
    do_init_dim : bint
        Whether to call U_init_dim_obs
    do_hx : bint
        Whether to initialize HX_p
    do_hxbar : bint
        Whether to initialize HXbar
    do_init_obs : bint
        Whether to initialize obs_p

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    with nogil:
        c__pdafobs_init(&step, &dim_p, &dim_ens, &dim_obs_p, &state_p[0],
                        &ens_p[0,0], pdaf_cb.c__init_dim_obs_pdaf,
                        pdaf_cb.c__obs_op_pdaf, pdaf_cb.c__init_obs_pdaf,
                        &screen, &debug, &do_ens_mean, &do_init_dim,
                        &do_hx, &do_hxbar, &do_init_obs)

    return dim_obs_p, state_p_np, ens_p_np


def obs_init_local(int  domain_p, int  step, int  dim_ens,
    py__init_dim_obs_l_pdaf, py__g2l_obs_pdaf, py__init_obs_l_pdaf, int  debug):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_ens : int
        Size of ensemble
    py__init_dim_obs_l_pdaf : Callable
        Init. dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    debug : int
        Flag for writing debug output

    Returns
    -------
    dim_obs_l : int
        Size of local observation vector
    dim_obs_f : int
        PE-local dimension of observation vector
    """
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    cdef int  dim_obs_l
    cdef int  dim_obs_f
    with nogil:
        c__pdafobs_init_local(&domain_p, &step, &dim_obs_l, &dim_obs_f,
                              &dim_ens, pdaf_cb.c__init_dim_obs_l_pdaf,
                              pdaf_cb.c__g2l_obs_pdaf,
                              pdaf_cb.c__init_obs_l_pdaf, &debug)

    return dim_obs_l, dim_obs_f


def obs_init_obsvars(int  step, int  dim_obs_p, py__init_obsvars_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        PE-local dimension of observation vector
    py__init_obsvars_pdaf : Callable
        Initialize vector of observation error variances

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_f : int
                Dimension of full observation vector

        Callback Returns
        ----------------
        var_f : ndarray[np.float64, ndim=1]
                vector of observation error variances
                Array shape: (dim_obs_f)


    Returns
    -------
    """
    pdaf_cb.init_obsvars_pdaf = <void*>py__init_obsvars_pdaf
    with nogil:
        c__pdafobs_init_obsvars(&step, &dim_obs_p, pdaf_cb.c__init_obsvars_pdaf)



def obs_dealloc():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdafobs_dealloc()



def obs_dealloc_local():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdafobs_dealloc_local()



def netf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_netf_init(&subtype, &param_int[0], &dim_pint,
                          &param_real[0], &dim_preal, &ensemblefilter,
                          &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def netf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_netf_alloc(&outflag)

    return outflag


def netf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_netf_config(&subtype, &verbose)

    return subtype


def netf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_netf_set_iparam(&id, &value, &flag)

    return flag


def netf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_netf_set_rparam(&id, &value, &flag)

    return flag


def netf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_netf_options()



def netf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_netf_memtime(&printtype)



def lenkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_lenkf_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def lenkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_lenkf_alloc(&outflag)

    return outflag


def lenkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_lenkf_config(&subtype, &verbose)

    return subtype


def lenkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lenkf_set_iparam(&id, &value, &flag)

    return flag


def lenkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lenkf_set_rparam(&id, &value, &flag)

    return flag


def lenkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_lenkf_options()



def lenkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_lenkf_memtime(&printtype)



def lseik_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_lseik_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def lseik_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_lseik_alloc(&outflag)

    return outflag


def lseik_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_lseik_config(&subtype, &verbose)

    return subtype


def lseik_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lseik_set_iparam(&id, &value, &flag)

    return flag


def lseik_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lseik_set_rparam(&id, &value, &flag)

    return flag


def lseik_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_lseik_options()



def lseik_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_lseik_memtime(&printtype)



def etkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_etkf_init(&subtype, &param_int[0], &dim_pint,
                          &param_real[0], &dim_preal, &ensemblefilter,
                          &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def etkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_etkf_alloc(&outflag)

    return outflag


def etkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_etkf_config(&subtype, &verbose)

    return subtype


def etkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_etkf_set_iparam(&id, &value, &flag)

    return flag


def etkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_etkf_set_rparam(&id, &value, &flag)

    return flag


def etkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_etkf_options()



def etkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_etkf_memtime(&printtype)



def lenkf_update(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__add_obs_err_pdaf, py__init_obs_pdaf,
    py__init_obs_covar_pdaf, py__prepoststep_pdaf, py__localize_covar_pdaf,
    int  screen, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of state ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Dimension of observation vector
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


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

    py__localize_covar_pdaf : Callable
        Apply localization to HP and HPH^T

        Callback Parameters
        -------------------
        dim_p : int
                pe-local state dimension
        dim_obs : int
                number of observations
        hp_p : ndarray[np.float64, ndim=2]
                pe local part of matrix hp
                Array shape: (dim_obs, dim_p)
        hph : ndarray[np.float64, ndim=2]
                matrix hph
                Array shape: (dim_obs, dim_obs)

        Callback Returns
        ----------------
        hp_p : ndarray[np.float64, ndim=2]
                pe local part of matrix hp
                Array shape: (dim_obs, dim_p)
        hph : ndarray[np.float64, ndim=2]
                matrix hph
                Array shape: (dim_obs, dim_obs)

    screen : int
        Verbosity flag
    subtype : int
        Specification of filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdaflenkf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                            &state_p[0], &ens_p[0,0],
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__add_obs_err_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__init_obs_covar_pdaf,
                            pdaf_cb.c__prepoststep_pdaf,
                            pdaf_cb.c__localize_covar_pdaf, &screen,
                            &subtype, &flag)

    return dim_obs_p, state_p_np, ens_p_np, flag


def pf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_pf_init(&subtype, &param_int[0], &dim_pint, &param_real[0],
                        &dim_preal, &ensemblefilter, &fixedbasis, &verbose,
                        &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def pf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_pf_alloc(&outflag)

    return outflag


def pf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_pf_config(&subtype, &verbose)

    return subtype


def pf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_pf_set_iparam(&id, &value, &flag)

    return flag


def pf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_pf_set_rparam(&id, &value, &flag)

    return flag


def pf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_pf_options()



def pf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_pf_memtime(&printtype)



def lknetf_ana_letkft(int  domain_p, int  step, int  dim_l,
    int  dim_obs_l, int  dim_ens, double [::1] state_l,
    double [::1,:] ens_l, double [::1,:] hz_l, double [::1] hxbar_l,
    double [::1] obs_l, double [::1,:] rndmat, double  forget,
    py__prodrinva_hyb_l_pdaf, py__init_obsvar_l_pdaf, double [::1] gamma,
    int  screen, int  type_forget, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    state_l : ndarray[np.float64, ndim=1]
        local forecast state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        PE-local full observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        local observed ens. mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_hyb_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain including hybrid weight

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        dim_ens : int
                Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        gamma : double
                Hybrid weight provided by PDAF
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, dim_ens)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, dim_ens)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, dim_ens)

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    screen : int
        Verbosity flag
    type_forget : int
        Type of forgetting factor
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        local forecast state
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        local weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        PE-local full observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv_l = ainv_l_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_l_np = np.asarray(hz_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] rndmat_np = np.asarray(rndmat, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gamma_np = np.asarray(gamma, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    with nogil:
        c__pdaf_lknetf_ana_letkft(&domain_p, &step, &dim_l, &dim_obs_l,
                                  &dim_ens, &state_l[0], &ainv_l[0,0],
                                  &ens_l[0,0], &hz_l[0,0], &hxbar_l[0],
                                  &obs_l[0], &rndmat[0,0], &forget,
                                  pdaf_cb.c__prodrinva_hyb_l_pdaf,
                                  pdaf_cb.c__init_obsvar_l_pdaf, &gamma[0],
                                  &screen, &type_forget, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hz_l_np, rndmat_np, forget, gamma_np, flag


def lknetf_ana_lnetf(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, double [::1,:] ens_l, double [::1,:] hx_l,
    double [::1,:] rndmat, double [::1] obs_l, py__likelihood_hyb_l_pdaf,
    int  cnt_small_svals, double [::1] n_eff_all, double [::1] gamma,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hx_l : ndarray[np.float64, ndim=2]
        local observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    py__likelihood_hyb_l_pdaf : Callable
        Compute observation likelihood for an ensemble member with hybrid weight

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                Input vector holding the local residual
                Array shape: (dim_obs_l)
        gamma : double
                Hybrid weight provided by PDAF

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                Input vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    cnt_small_svals : int
        Number of small eigen values
    n_eff_all : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    cnt_small_svals : int
        Number of small eigen values
    n_eff_all : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] n_eff_all_np = np.asarray(n_eff_all, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gamma_np = np.asarray(gamma, dtype=np.float64, order="F")
    pdaf_cb.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    with nogil:
        c__pdaf_lknetf_ana_lnetf(&domain_p, &step, &dim_l, &dim_obs_l,
                                 &dim_ens, &ens_l[0,0], &hx_l[0,0],
                                 &rndmat[0,0], &obs_l[0],
                                 pdaf_cb.c__likelihood_hyb_l_pdaf,
                                 &cnt_small_svals, &n_eff_all[0],
                                 &gamma[0], &screen, &flag)

    return ens_l_np, cnt_small_svals, n_eff_all_np, gamma_np, flag


def enkf_ana_rsm(int  step, int  dim_p, int  dim_obs_p, int dim_obs, int  dim_ens,
    int  rank_ana, double [::1] state_p, double [::1,:] ens_p,
    double [::1,:] hx_p, double [::1] hxbar_p, double [::1] obs_p,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf, int  screen, int  debug,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_obs : int
        Global dimension of observation vector
    dim_ens : int
        Size of state ensemble
    rank_ana : int
        Rank to be considered for inversion of HPH
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Dimension of observation vector
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    with nogil:
        c__pdaf_enkf_ana_rsm(&step, &dim_p, &dim_obs_p, &dim_obs, &dim_ens,
                             &rank_ana, &state_p[0], &ens_p[0,0],
                             &hx_p[0,0], &hxbar_p[0], &obs_p[0],
                             pdaf_cb.c__add_obs_err_pdaf,
                             pdaf_cb.c__init_obs_covar_pdaf, &screen,
                             &debug, &flag)

    return state_p_np, ens_p_np, flag


def lknetf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_lknetf_init(&subtype, &param_int[0], &dim_pint,
                            &param_real[0], &dim_preal, &ensemblefilter,
                            &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def lknetf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_lknetf_alloc(&outflag)

    return outflag


def lknetf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_lknetf_config(&subtype, &verbose)

    return subtype


def lknetf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lknetf_set_iparam(&id, &value, &flag)

    return flag


def lknetf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lknetf_set_rparam(&id, &value, &flag)

    return flag


def lknetf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_lknetf_options()



def lknetf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_lknetf_memtime(&printtype)



def lknetf_alpha_neff(int  dim_ens, double [::1] weights, double  hlimit,
    double  alpha):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_ens : int
        Size of ensemble
    weights : ndarray[np.float64, ndim=1]
        Weights
        Array shape: (dim_ens)
    hlimit : double
        Minimum of n_eff / N
    alpha : double
        hybrid weight

    Returns
    -------
    alpha : double
        hybrid weight
    """
    with nogil:
        c__pdaf_lknetf_alpha_neff(&dim_ens, &weights[0], &hlimit, &alpha)

    return alpha


def lknetf_compute_gamma(int  domain_p, int  step, int  dim_obs_l,
    int  dim_ens, double [::1,:] hx_l, double [::1] hxbar_l,
    double [::1] obs_l, int  type_hyb, double  hyb_g, double  hyb_k,
    double [::1] gamma, double [::1] n_eff_out, double [::1] skew_mabs,
    double [::1] kurt_mabs, py__likelihood_l_pdaf, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    hx_l : ndarray[np.float64, ndim=2]
        local observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        local mean observed ensemble
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    type_hyb : int
        Type of hybrid weight
    hyb_g : double
        Prescribed hybrid weight for state transformation
    hyb_k : double
        Hybrid weight for covariance transformation
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    n_eff_out : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    skew_mabs : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    kurt_mabs : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    py__likelihood_l_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    n_eff_out : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    skew_mabs : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    kurt_mabs : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gamma_np = np.asarray(gamma, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] n_eff_out_np = np.asarray(n_eff_out, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] skew_mabs_np = np.asarray(skew_mabs, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] kurt_mabs_np = np.asarray(kurt_mabs, dtype=np.float64, order="F")
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    with nogil:
        c__pdaf_lknetf_compute_gamma(&domain_p, &step, &dim_obs_l,
                                     &dim_ens, &hx_l[0,0], &hxbar_l[0],
                                     &obs_l[0], &type_hyb, &hyb_g, &hyb_k,
                                     &gamma[0], &n_eff_out[0],
                                     &skew_mabs[0], &kurt_mabs[0],
                                     pdaf_cb.c__likelihood_l_pdaf, &screen,
                                     &flag)

    return gamma_np, n_eff_out_np, skew_mabs_np, kurt_mabs_np, flag


def lknetf_set_gamma(int  domain_p, int  dim_obs_l, int  dim_ens,
    double [::1,:] hx_l, double [::1] hxbar_l, double [::1] weights,
    int  type_hyb, double  hyb_g, double  hyb_k, double [::1] gamma,
    double [::1] n_eff_out, double [::1] maskew, double [::1] makurt,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    hx_l : ndarray[np.float64, ndim=2]
        local observed state ens.
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        local mean observed ensemble
        Array shape: (dim_obs_l)
    weights : ndarray[np.float64, ndim=1]
        Weight vector
        Array shape: (dim_ens)
    type_hyb : int
        Type of hybrid weight
    hyb_g : double
        Prescribed hybrid weight for state transformation
    hyb_k : double
        Scale factor kappa (for type_hyb 3 and 4)
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    n_eff_out : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    maskew : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    makurt : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    gamma : ndarray[np.float64, ndim=1]
        Hybrid weight for state transformation
        Array shape: (1)
    n_eff_out : ndarray[np.float64, ndim=1]
        Effective ensemble size
        Array shape: (1)
    maskew : ndarray[np.float64, ndim=1]
        Mean absolute skewness
        Array shape: (1)
    makurt : ndarray[np.float64, ndim=1]
        Mean absolute kurtosis
        Array shape: (1)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] gamma_np = np.asarray(gamma, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] n_eff_out_np = np.asarray(n_eff_out, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] maskew_np = np.asarray(maskew, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] makurt_np = np.asarray(makurt, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_lknetf_set_gamma(&domain_p, &dim_obs_l, &dim_ens,
                                 &hx_l[0,0], &hxbar_l[0], &weights[0],
                                 &type_hyb, &hyb_g, &hyb_k, &gamma[0],
                                 &n_eff_out[0], &maskew[0], &makurt[0],
                                 &screen, &flag)

    return gamma_np, n_eff_out_np, maskew_np, makurt_np, flag


def lknetf_reset_gamma(double  gamma_in):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    gamma_in : double
        Prescribed hybrid weight

    Returns
    -------
    """
    with nogil:
        c__pdaf_lknetf_reset_gamma(&gamma_in)



def hyb3dvar_analysis_cvt(int  step, int  dim_p, int  dim_obs_p,
    int  dim_ens, int  dim_cvec, int  dim_cvec_ens, double  beta_3dvar,
    double [::1] state_p, double [::1,:] ens_p, double [::1] state_inc_p,
    double [::1] hxbar_p, double [::1] obs_p, py__prodrinva_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  screen, int  type_opt,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    dim_cvec : int
        Size of control vector (parameterized part)
    dim_cvec_ens : int
        Size of control vector (ensemble part)
    beta_3dvar : double
        Hybrid weight for hybrid 3D-Var
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    state_inc_p : ndarray[np.float64, ndim=1]
        PE-local state analysis increment
        Array shape: (dim_p)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix (parameterized)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    screen : int
        Verbosity flag
    type_opt : int
        Type of minimizer for 3DVar
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    state_inc_p : ndarray[np.float64, ndim=1]
        PE-local state analysis increment
        Array shape: (dim_p)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_inc_p_np = np.asarray(state_inc_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdafhyb3dvar_analysis_cvt(&step, &dim_p, &dim_obs_p, &dim_ens,
                                     &dim_cvec, &dim_cvec_ens, &beta_3dvar,
                                     &state_p[0], &ens_p[0,0],
                                     &state_inc_p[0], &hxbar_p[0],
                                     &obs_p[0], pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf, &screen,
                                     &type_opt, &debug, &flag)

    return state_p_np, ens_p_np, state_inc_p_np, flag


def _3dvar_analysis_cvt(int  step, int  dim_p, int  dim_obs_p,
    int  dim_cvec, double [::1] hxbar_p, double [::1] obs_p,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, int  screen, int  type_opt,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_cvec : int
        Size of control vector
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    screen : int
        Verbosity flag
    type_opt : int
        Type of minimizer for 3DVar
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    with nogil:
        c__pdaf3dvar_analysis_cvt(&step, &dim_p, &dim_obs_p, &dim_cvec,
                                  &state_p[0], &hxbar_p[0], &obs_p[0],
                                  pdaf_cb.c__prodrinva_pdaf,
                                  pdaf_cb.c__cvt_pdaf,
                                  pdaf_cb.c__cvt_adj_pdaf,
                                  pdaf_cb.c__obs_op_lin_pdaf,
                                  pdaf_cb.c__obs_op_adj_pdaf, &screen,
                                  &type_opt, &debug, &flag)

    return state_p_np, flag


def lestkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_lestkf_init(&subtype, &param_int[0], &dim_pint,
                            &param_real[0], &dim_preal, &ensemblefilter,
                            &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def lestkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_lestkf_alloc(&outflag)

    return outflag


def lestkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_lestkf_config(&subtype, &verbose)

    return subtype


def lestkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lestkf_set_iparam(&id, &value, &flag)

    return flag


def lestkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lestkf_set_rparam(&id, &value, &flag)

    return flag


def lestkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_lestkf_options()



def lestkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_lestkf_memtime(&printtype)



def seik_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    int  rank, double [::1] state_p, double [::1,:] uinv,
    double [::1,:] ens_p, double [::1,:] hl_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  debug,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of eigenvalue matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of eigenvalue matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_p_np = np.asarray(hl_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_seik_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                         &state_p[0], &uinv[0,0], &ens_p[0,0], &hl_p[0,0],
                         &hxbar_p[0], &obs_p[0], &forget,
                         pdaf_cb.c__prodrinva_pdaf, &debug, &flag)

    return state_p_np, uinv_np, ens_p_np, hl_p_np, flag


def seik_resample(int  subtype, int  dim_p, int  dim_ens, int  rank,
    double [::1,:] uinv, double [::1] state_p, double [::1,:] enst_p,
    int  type_sqrt, int  type_trans, int  nm1vsn, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Filter subtype
    dim_p : int
        PE-local state dimension
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    enst_p : ndarray[np.float64, ndim=2]
        PE-local ensemble times T
        Array shape: (dim_p, dim_ens)
    type_sqrt : int
        Type of square-root of A
    type_trans : int
        Type of ensemble transformation
    nm1vsn : int
        Flag which normalization of P ist used in SEIK
    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    enst_p : ndarray[np.float64, ndim=2]
        PE-local ensemble times T
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] enst_p_np = np.asarray(enst_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_seik_resample(&subtype, &dim_p, &dim_ens, &rank,
                              &uinv[0,0], &state_p[0], &enst_p[0,0],
                              &type_sqrt, &type_trans, &nm1vsn, &screen, &flag)

    return uinv_np, state_p_np, enst_p_np, flag


def lseik_ana(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, int  rank, double [::1] state_l, double [::1,:] uinv_l,
    double [::1,:] ens_l, double [::1,:] hl_l, double [::1] hxbar_l,
    double [::1] obs_l, double  forget, py__prodrinva_l_pdaf, int  screen,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_l : ndarray[np.float64, ndim=1]
        State on local analysis domain
        Array shape: (dim_l)
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        State on local analysis domain
        Array shape: (dim_l)
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_l_np = np.asarray(uinv_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_l_np = np.asarray(hl_l, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_lseik_ana(&domain_p, &step, &dim_l, &dim_obs_l, &dim_ens,
                          &rank, &state_l[0], &uinv_l[0,0], &ens_l[0,0],
                          &hl_l[0,0], &hxbar_l[0], &obs_l[0], &forget,
                          pdaf_cb.c__prodrinva_l_pdaf, &screen, &debug, &flag)

    return state_l_np, uinv_l_np, hl_l_np, forget, flag


def lseik_resample(int  domain_p, int  subtype, int  dim_l, int  dim_ens,
    int  rank, double [::1,:] uinv_l, double [::1] state_l,
    double [::1,:] ens_l, double [::1,:] omegat_in, int  type_sqrt,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    subtype : int
        Specification of filter subtype
    dim_l : int
        State dimension on local analysis domain
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_l : ndarray[np.float64, ndim=1]
        Local model state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    omegat_in : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank, dim_ens)
    type_sqrt : int
        Type of square-root of A
    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    uinv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_l : ndarray[np.float64, ndim=1]
        Local model state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    omegat_in : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_l_np = np.asarray(uinv_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] omegat_in_np = np.asarray(omegat_in, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_lseik_resample(&domain_p, &subtype, &dim_l, &dim_ens,
                               &rank, &uinv_l[0,0], &state_l[0],
                               &ens_l[0,0], &omegat_in[0,0], &type_sqrt,
                               &screen, &flag)

    return uinv_l_np, state_l_np, ens_l_np, omegat_in_np, flag


def prepost(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

        Callback Parameters
        -------------------
        dim_p : int
                pe-local state dimension
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                local state vector
                Array shape: (dim_p)

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

    py__next_observation_pdaf : Callable
        Routine to provide time step, time and dimensionof next observation

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_prepost(pdaf_cb.c__collect_state_pdaf,
                        pdaf_cb.c__distribute_state_pdaf,
                        pdaf_cb.c__prepoststep_pdaf,
                        pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def enkf_update(int  step, int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ens_p, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__add_obs_err_pdaf, py__init_obs_pdaf, py__init_obs_covar_pdaf,
    py__prepoststep_pdaf, int  screen, int  subtype, int  dim_lag,
    double [::1,:,:] sens_p, int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of state ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Dimension of observation vector
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


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

    screen : int
        Verbosity flag
    subtype : int
        Specification of filter subtype
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafenkf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                           &state_p[0], &ens_p[0,0],
                           pdaf_cb.c__init_dim_obs_pdaf,
                           pdaf_cb.c__obs_op_pdaf,
                           pdaf_cb.c__add_obs_err_pdaf,
                           pdaf_cb.c__init_obs_pdaf,
                           pdaf_cb.c__init_obs_covar_pdaf,
                           pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                           &dim_lag, &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_p, state_p_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def init_parallel(int  dim_ens, bint  ensemblefilter, bint  fixedbasis,
    int  comm_model, int  in_comm_filter, int  in_comm_couple,
    int  in_n_modeltasks, int  in_task_id, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_ens : int
        Rank of covar matrix/ensemble size
    ensemblefilter : bint
        Is the filter ensemble-based?
    fixedbasis : bint
        Run with fixed error-space basis?
    comm_model : int
        Model communicator (not shared)
    in_comm_filter : int
        Filter communicator
    in_comm_couple : int
        Coupling communicator
    in_n_modeltasks : int
        Number of model tasks
    in_task_id : int
        Task ID of current PE
    screen : int
        Whether screen information is shown
    flag : int
        Status flag

    Returns
    -------
    dim_ens : int
        Rank of covar matrix/ensemble size
    flag : int
        Status flag
    """
    with nogil:
        c__pdaf_init_parallel(&dim_ens, &ensemblefilter, &fixedbasis,
                              &comm_model, &in_comm_filter,
                              &in_comm_couple, &in_n_modeltasks,
                              &in_task_id, &screen, &flag)

    return dim_ens, flag


def seik_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_seik_init(&subtype, &param_int[0], &dim_pint,
                          &param_real[0], &dim_preal, &ensemblefilter,
                          &fixedbasis, &verbose, &outflag)

    return param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def seik_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_seik_alloc(&outflag)

    return outflag


def seik_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_seik_config(&subtype, &verbose)

    return subtype


def seik_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_seik_set_iparam(&id, &value, &flag)

    return flag


def seik_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_seik_set_rparam(&id, &value, &flag)

    return flag


def seik_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_seik_options()



def seik_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_seik_memtime(&printtype)



def netf_update(int  step, int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__likelihood_pdaf,
    py__prepoststep_pdaf, int  screen, int  subtype, int  dim_lag,
    double [::1,:,:] sens_p, int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__likelihood_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        resid : ndarray[np.float64, ndim=1]
                Input vector holding the residual
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        likely : double
                Output value of the likelihood

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafnetf_update(&step, &dim_p, &dim_obs_p, &dim_ens,
                           &state_p[0], &ainv[0,0], &ens_p[0,0],
                           pdaf_cb.c__init_dim_obs_pdaf,
                           pdaf_cb.c__obs_op_pdaf,
                           pdaf_cb.c__init_obs_pdaf,
                           pdaf_cb.c__likelihood_pdaf,
                           pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                           &dim_lag, &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def seik_ana_newt(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    int  rank, double [::1] state_p, double [::1,:] uinv,
    double [::1,:] ens_p, double [::1,:] hl_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of eigenvalue matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of eigenvalue matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_p_np = np.asarray(hl_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_seik_ana_newt(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                              &state_p[0], &uinv[0,0], &ens_p[0,0],
                              &hl_p[0,0], &hxbar_p[0], &obs_p[0], &forget,
                              pdaf_cb.c__prodrinva_pdaf, &screen, &debug, &flag)

    return state_p_np, uinv_np, ens_p_np, hl_p_np, flag


def seik_resample_newt(int  subtype, int  dim_p, int  dim_ens, int  rank,
    double [::1,:] uinv, double [::1] state_p, double [::1,:] ens_p,
    int  type_sqrt, int  type_trans, int  nm1vsn, int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Filter subtype
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    type_sqrt : int
        Type of square-root of A
    type_trans : int
        Type of ensemble transformation
    nm1vsn : int
        Flag which normalization of P ist used in SEIK
    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_seik_resample_newt(&subtype, &dim_p, &dim_ens, &rank,
                                   &uinv[0,0], &state_p[0], &ens_p[0,0],
                                   &type_sqrt, &type_trans, &nm1vsn,
                                   &screen, &flag)

    return uinv_np, state_p_np, ens_p_np, flag


def lenkf_ana_rsm(int  step, int  dim_p, int  dim_obs_p, int dim_obs, int  dim_ens,
    int  rank_ana, double [::1] state_p, double [::1,:] ens_p,
    double [::1,:] hx_p, double [::1] hxbar_p, double [::1] obs_p,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf, py__localize_covar_pdaf,
    int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_obs: int
        Global dimension of observation vector
    dim_ens : int
        Size of state ensemble
    rank_ana : int
        Rank to be considered for inversion of HPH
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Dimension of observation vector
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Matrix to that observation covariance R is added
                Array shape: (dim_obs_p,dim_obs_p)

    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


    py__localize_covar_pdaf : Callable
        Apply localization to HP and HPH^T

        Callback Parameters
        -------------------
        dim_p : int
                pe-local state dimension
        dim_obs : int
                number of observations
        hp_p : ndarray[np.float64, ndim=2]
                pe local part of matrix hp
                Array shape: (dim_obs, dim_p)
        hph : ndarray[np.float64, ndim=2]
                matrix hph
                Array shape: (dim_obs, dim_obs)

        Callback Returns
        ----------------
        hp_p : ndarray[np.float64, ndim=2]
                pe local part of matrix hp
                Array shape: (dim_obs, dim_p)
        hph : ndarray[np.float64, ndim=2]
                matrix hph
                Array shape: (dim_obs, dim_obs)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    with nogil:
        c__pdaf_lenkf_ana_rsm(&step, &dim_p, &dim_obs_p, &dim_obs, &dim_ens,
                              &rank_ana, &state_p[0], &ens_p[0,0],
                              &hx_p[0,0], &hxbar_p[0], &obs_p[0],
                              pdaf_cb.c__add_obs_err_pdaf,
                              pdaf_cb.c__init_obs_covar_pdaf,
                              pdaf_cb.c__localize_covar_pdaf, &screen,
                              &debug, &flag)

    return state_p_np, ens_p_np, flag


def lestkf_ana(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, int  rank, double [::1] state_l, double [::1,:] ainv_l,
    double [::1,:] ens_l, double [::1,:] hl_l, double [::1] hxbar_l,
    double [::1] obs_l, double [::1,:] omegat_in, double  forget,
    py__prodrinva_l_pdaf, int  envar_mode, int  type_sqrt,
    double [::1,:] ta, int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_l : ndarray[np.float64, ndim=1]
        state on local analysis domain
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    omegat_in : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    envar_mode : int
        Flag whether routine is called from 3DVar for special functionality
    type_sqrt : int
        Type of square-root of A
    ta : ndarray[np.float64, ndim=2]
        Ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        state on local analysis domain
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    forget : double
        Forgetting factor
    ta : ndarray[np.float64, ndim=2]
        Ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.asarray(ainv_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_l_np = np.asarray(hl_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ta_np = np.asarray(ta, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_lestkf_ana(&domain_p, &step, &dim_l, &dim_obs_l, &dim_ens,
                           &rank, &state_l[0], &ainv_l[0,0], &ens_l[0,0],
                           &hl_l[0,0], &hxbar_l[0], &obs_l[0],
                           &omegat_in[0,0], &forget,
                           pdaf_cb.c__prodrinva_l_pdaf, &envar_mode,
                           &type_sqrt, &ta[0,0], &screen, &debug, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hl_l_np, forget, ta_np, flag


def lestkf_update(int  step, int  dim_p, int  dim_ens, int  rank,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__prepoststep_pdaf, int  screen,
    int  subtype, int  envar_mode, int  dim_lag, double [::1,:,:] sens_p,
    int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Compute product of R^(-1) with HV

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    envar_mode : int
        Flag whether routine is called from 3DVar for special functionality
    dim_lag : int
        Number of past time instances for smoother
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaflestkf_update(&step, &dim_p, &dim_obs_f, &dim_ens, &rank,
                             &state_p[0], &ainv[0,0], &ens_p[0,0],
                             pdaf_cb.c__init_dim_obs_pdaf,
                             pdaf_cb.c__obs_op_pdaf,
                             pdaf_cb.c__init_obs_pdaf,
                             pdaf_cb.c__init_obs_l_pdaf,
                             pdaf_cb.c__prodrinva_l_pdaf,
                             pdaf_cb.c__init_n_domains_p_pdaf,
                             pdaf_cb.c__init_dim_l_pdaf,
                             pdaf_cb.c__init_dim_obs_l_pdaf,
                             pdaf_cb.c__g2l_state_pdaf,
                             pdaf_cb.c__l2g_state_pdaf,
                             pdaf_cb.c__g2l_obs_pdaf,
                             pdaf_cb.c__init_obsvar_pdaf,
                             pdaf_cb.c__init_obsvar_l_pdaf,
                             pdaf_cb.c__prepoststep_pdaf, &screen,
                             &subtype, &envar_mode, &dim_lag,
                             &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, sens_p_np, cnt_maxlag, flag


def lnetf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_lnetf_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def lnetf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_lnetf_alloc(&outflag)

    return outflag


def lnetf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_lnetf_config(&subtype, &verbose)

    return subtype


def lnetf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lnetf_set_iparam(&id, &value, &flag)

    return flag


def lnetf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_lnetf_set_rparam(&id, &value, &flag)

    return flag


def lnetf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_lnetf_options()



def lnetf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_lnetf_memtime(&printtype)



def enkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_enkf_init(&subtype, &param_int[0], &dim_pint,
                          &param_real[0], &dim_preal, &ensemblefilter,
                          &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def enkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_enkf_alloc(&outflag)

    return outflag


def enkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_enkf_config(&subtype, &verbose)

    return subtype


def enkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_enkf_set_iparam(&id, &value, &flag)

    return flag


def enkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_enkf_set_rparam(&id, &value, &flag)

    return flag


def enkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_enkf_options()



def enkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_enkf_memtime(&printtype)



def enkf_gather_resid(int  dim_obs, int  dim_obs_p, int  dim_ens,
    double [::1,:] resid_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs : int
        Global observation dimension
    dim_obs_p : int
        PE-local observation dimension
    dim_ens : int
        Ensemble size
    resid_p : ndarray[np.float64, ndim=2]
        PE-local residual matrix
        Array shape: (dim_obs_p, dim_ens)

    Returns
    -------
    resid : ndarray[np.float64, ndim=2]
        Global residual matrix
        Array shape: (dim_obs, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] resid_np = np.zeros((dim_obs, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] resid = resid_np
    with nogil:
        c__pdaf_enkf_gather_resid(&dim_obs, &dim_obs_p, &dim_ens,
                                  &resid_p[0,0], &resid[0,0])

    return resid_np


def enkf_obs_ensemble(int  step, int  dim_obs_p, int  dim_obs,
    int  dim_ens, double [::1] obs_p, py__init_obs_covar_pdaf, int  screen,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Local dimension of current observation
    dim_obs : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    py__init_obs_covar_pdaf : Callable
        Initialize observation error covariance matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs : int
                Global size of observation vector
        dim_obs_p : int
                Size of process-local observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Process-local vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        covar : ndarray[np.float64, ndim=2]
                Observation error covariance matrix
                Array shape: (dim_obs_p,dim_obs_p)
        isdiag : bint


    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    obsens_p : ndarray[np.float64, ndim=2]
        PE-local obs. ensemble
        Array shape: (dim_obs_p,dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] obsens_p_np = np.zeros((dim_obs_p,dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] obsens_p = obsens_p_np
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    with nogil:
        c__pdaf_enkf_obs_ensemble(&step, &dim_obs_p, &dim_obs, &dim_ens,
                                  &obsens_p[0,0], &obs_p[0],
                                  pdaf_cb.c__init_obs_covar_pdaf, &screen,
                                  &flag)

    return obsens_p_np, flag


def pf_update(int  step, int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__likelihood_pdaf,
    py__prepoststep_pdaf, int  screen, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__likelihood_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        resid : ndarray[np.float64, ndim=1]
                Input vector holding the residual
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        likely : double
                Output value of the likelihood

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafpf_update(&step, &dim_p, &dim_obs_p, &dim_ens, &state_p[0],
                         &ainv[0,0], &ens_p[0,0],
                         pdaf_cb.c__init_dim_obs_pdaf,
                         pdaf_cb.c__obs_op_pdaf, pdaf_cb.c__init_obs_pdaf,
                         pdaf_cb.c__likelihood_pdaf,
                         pdaf_cb.c__prepoststep_pdaf, &screen, &subtype, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag


def generate_rndmat(int  dim, int  mattype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        Size of matrix rndmat
    mattype : int
        Select type of random matrix:

    Returns
    -------
    rndmat : ndarray[np.float64, ndim=2]
        Matrix
        Array shape: (dim, dim)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] rndmat_np = np.zeros((dim, dim), dtype=np.float64, order="F")
    cdef double [::1,:] rndmat = rndmat_np
    with nogil:
        c__pdaf_generate_rndmat(&dim, &rndmat[0,0], &mattype)

    return rndmat_np


def print_domain_stats(int  n_domains_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    n_domains_p : int
        Number of PE-local analysis domains

    Returns
    -------
    """
    with nogil:
        c__pdaf_print_domain_stats(&n_domains_p)



def init_local_obsstats():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_init_local_obsstats()



def incr_local_obsstats(int  dim_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs_l : int
        Number of locally assimilated observations

    Returns
    -------
    """
    with nogil:
        c__pdaf_incr_local_obsstats(&dim_obs_l)



def print_local_obsstats(int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    n_domains_with_obs : int

    """
    cdef int  n_domains_with_obs
    with nogil:
        c__pdaf_print_local_obsstats(&screen, &n_domains_with_obs)

    return n_domains_with_obs


def seik_matrixt(int  dim, int  dim_ens, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        dimension of states
    dim_ens : int
        Size of ensemble
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)

    Returns
    -------
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_np = np.asarray(a, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_seik_matrixt(&dim, &dim_ens, &a[0,0])

    return a_np


def seik_ttimesa(int  rank, int  dim_col, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    rank : int
        Rank of initial covariance matrix
    dim_col : int
        Number of columns in A and B
    a : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (rank, dim_col)

    Returns
    -------
    b : ndarray[np.float64, ndim=2]
        Output matrix (TA)
        Array shape: (rank+1, dim_col)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] b_np = np.zeros((rank+1, dim_col), dtype=np.float64, order="F")
    cdef double [::1,:] b = b_np
    with nogil:
        c__pdaf_seik_ttimesa(&rank, &dim_col, &a[0,0], &b[0,0])

    return b_np


def seik_omega(int  rank, double [::1,:] omega, int  omegatype, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    rank : int
        Approximated rank of covar matrix
    omega : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank+1, rank)
    omegatype : int
        Select type of Omega:
    screen : int
        Verbosity flag

    Returns
    -------
    omega : ndarray[np.float64, ndim=2]
        Matrix Omega
        Array shape: (rank+1, rank)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] omega_np = np.asarray(omega, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_seik_omega(&rank, &omega[0,0], &omegatype, &screen)

    return omega_np


def seik_uinv(int  rank, double [::1,:] uinv):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    rank : int
        Rank of initial covariance matrix
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)

    Returns
    -------
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (rank, rank)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_seik_uinv(&rank, &uinv[0,0])

    return uinv_np


def ens_omega(int [::1] seed, int  r, int  dim_ens, double [::1,:] omega,
    double  norm, int  otype, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    seed : ndarray[np.intc, ndim=1]
        Seed for random number generation
        Array shape: (4)
    r : int
        Approximated rank of covar matrix
    dim_ens : int
        Ensemble size
    omega : ndarray[np.float64, ndim=2]
        Random matrix
        Array shape: (dim_ens,r)
    norm : double
        Norm for ensemble transformation
    otype : int
        Type of Omega:
    screen : int
        Control verbosity

    Returns
    -------
    omega : ndarray[np.float64, ndim=2]
        Random matrix
        Array shape: (dim_ens,r)
    norm : double
        Norm for ensemble transformation
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] omega_np = np.asarray(omega, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_ens_omega(&seed[0], &r, &dim_ens, &omega[0,0], &norm,
                          &otype, &screen)

    return omega_np, norm


def estkf_omegaa(int  rank, int  dim_col, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    rank : int
        Rank of initial covariance matrix
    dim_col : int
        Number of columns in A and B
    a : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (rank, dim_col)

    Returns
    -------
    b : ndarray[np.float64, ndim=2]
        Output matrix (TA)
        Array shape: (rank+1, dim_col)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] b_np = np.zeros((rank+1, dim_col), dtype=np.float64, order="F")
    cdef double [::1,:] b = b_np
    with nogil:
        c__pdaf_estkf_omegaa(&rank, &dim_col, &a[0,0], &b[0,0])

    return b_np


def estkf_aomega(int  dim, int  dim_ens, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        dimension of states
    dim_ens : int
        Size of ensemble
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)

    Returns
    -------
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_np = np.asarray(a, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_estkf_aomega(&dim, &dim_ens, &a[0,0])

    return a_np


def subtract_rowmean(int  dim, int  dim_ens, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        dimension of states
    dim_ens : int
        Size of ensemble
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)

    Returns
    -------
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_np = np.asarray(a, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_subtract_rowmean(&dim, &dim_ens, &a[0,0])

    return a_np


def subtract_colmean(int  dim_ens, int  dim, double [::1,:] a):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_ens : int
        Rank of initial covariance matrix
    dim : int
        Number of columns in A and B
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim_ens, dim)

    Returns
    -------
    a : ndarray[np.float64, ndim=2]
        Input/output matrix
        Array shape: (dim_ens, dim)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_np = np.asarray(a, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_subtract_colmean(&dim_ens, &dim, &a[0,0])

    return a_np


def add_particle_noise(int  dim_p, int  dim_ens, double [::1] state_p,
    double [::1,:] ens_p, int  type_noise, double  noise_amp, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        State dimension
    dim_ens : int
        Number of particles
    state_p : ndarray[np.float64, ndim=1]
        State vector (not filled)
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        Ensemble array
        Array shape: (dim_p, dim_ens)
    type_noise : int
        Type of noise
    noise_amp : double
        Noise amplitude
    screen : int
        Verbosity flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        State vector (not filled)
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        Ensemble array
        Array shape: (dim_p, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_add_particle_noise(&dim_p, &dim_ens, &state_p[0],
                                   &ens_p[0,0], &type_noise, &noise_amp,
                                   &screen)

    return state_p_np, ens_p_np


def inflate_weights(int  screen, int  dim_ens, double  alpha,
    double [::1] weights):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    screen : int
        verbosity flag
    dim_ens : int
        Ensemble size
    alpha : double
        Minimum limit of n_eff / N
    weights : ndarray[np.float64, ndim=1]
        weights (before and after inflation)
        Array shape: (dim_ens)

    Returns
    -------
    weights : ndarray[np.float64, ndim=1]
        weights (before and after inflation)
        Array shape: (dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] weights_np = np.asarray(weights, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_inflate_weights(&screen, &dim_ens, &alpha, &weights[0])

    return weights_np


def inflate_ens(int  dim, int  dim_ens, double [::1] meanstate,
    double [::1,:] ens, double  forget, bint  do_ensmean):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim : int
        dimension of states
    dim_ens : int
        Size of ensemble
    meanstate : ndarray[np.float64, ndim=1]
        state vector to hold ensemble mean
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        Input/output ensemble matrix
        Array shape: (dim, dim_ens)
    forget : double
        Forgetting factor
    do_ensmean : bint
        Whether to compute the ensemble mean state

    Returns
    -------
    meanstate : ndarray[np.float64, ndim=1]
        state vector to hold ensemble mean
        Array shape: (dim)
    ens : ndarray[np.float64, ndim=2]
        Input/output ensemble matrix
        Array shape: (dim, dim_ens)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] meanstate_np = np.asarray(meanstate, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_np = np.asarray(ens, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_inflate_ens(&dim, &dim_ens, &meanstate[0], &ens[0,0],
                            &forget, &do_ensmean)

    return meanstate_np, ens_np


def alloc(int  dim_p, int  dim_ens, int  dim_ens_task, int  dim_es,
    int  statetask, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        Size of state vector
    dim_ens : int
        Ensemble size
    dim_ens_task : int
        Ensemble size handled by a model task
    dim_es : int
        Dimension of error space (size of Ainv)
    statetask : int
        Task ID forecasting a single state
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_alloc(&dim_p, &dim_ens, &dim_ens_task, &dim_es,
                      &statetask, &outflag)

    return outflag

def alloc_sens(int dim_p, int dim_ens, int dim_lag, int outflag):
    """Allocate PDAF smoother array

    Parameters
    ----------
    dim_p : int
        Dimension of PE-local state vector
    dim_ens: int
        Ensemble size
    dim_lag: int
        Smoothing lag
    outflag : int
        Status flag

    Returns
    -------
    outflag: int
        Status flag
    """
    with nogil:
        c__pdaf_alloc_sens(&dim_p, &dim_ens, &dim_lag, &outflag)
    return outflag

def alloc_bias(int dim_bias_p, int outflag):
    """Allocate PDAF bias array

    Parameters
    ----------
    dim_bias_p : int
        Dimension of PE-local bias vector
    outflag : int
        Status flag

    Returns
    -------
    outflag: int
        Status flag
    """
    with nogil:
        c__pdaf_alloc_bias(&dim_bias_p, &outflag)
    return outflag

def smoothing(int  dim_p, int  dim_ens, int  dim_lag, double [::1,:] ainv,
    double [::1,:,:] sens_p, int  cnt_maxlag, double  forget, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ainv : ndarray[np.float64, ndim=2]
        Weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    forget : double
        Forgetting factor
    screen : int
        Verbosity flag

    Returns
    -------
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_smoothing(&dim_p, &dim_ens, &dim_lag, &ainv[0,0],
                          &sens_p[0,0,0], &cnt_maxlag, &forget, &screen)

    return sens_p_np, cnt_maxlag


def smoothing_local(int  domain_p, int  step, int  dim_p, int  dim_l,
    int  dim_ens, int  dim_lag, double [::1,:] ainv, double [::1,:] ens_l,
    double [::1,:,:] sens_p, int  cnt_maxlag, py__g2l_state_pdaf,
    py__l2g_state_pdaf, double  forget, int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_l : int
        State dimension on local analysis domain
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ainv : ndarray[np.float64, ndim=2]
        Weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        local past ensemble (temporary)
        Array shape: (dim_l, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    forget : double
        Forgetting factor
    screen : int
        Verbosity flag

    Returns
    -------
    ens_l : ndarray[np.float64, ndim=2]
        local past ensemble (temporary)
        Array shape: (dim_l, dim_ens)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    with nogil:
        c__pdaf_smoothing_local(&domain_p, &step, &dim_p, &dim_l, &dim_ens,
                                &dim_lag, &ainv[0,0], &ens_l[0,0],
                                &sens_p[0,0,0], &cnt_maxlag,
                                pdaf_cb.c__g2l_state_pdaf,
                                pdaf_cb.c__l2g_state_pdaf, &forget, &screen)

    return ens_l_np, sens_p_np, cnt_maxlag


def smoother_shift(int  dim_p, int  dim_ens, int  dim_lag,
    double [::1,:,:] ens_p, double [::1,:,:] sens_p, int  cnt_maxlag,
    int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_lag : int
        Number of past time instances for smoother
    ens_p : ndarray[np.float64, ndim=3]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens, 1)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    screen : int
        Verbosity flag

    Returns
    -------
    ens_p : ndarray[np.float64, ndim=3]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens, 1)
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count available number of time steps for smoothing
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaf_smoother_shift(&dim_p, &dim_ens, &dim_lag, &ens_p[0,0,0],
                               &sens_p[0,0,0], &cnt_maxlag, &screen)

    return ens_p_np, sens_p_np, cnt_maxlag


def lknetf_update_sync(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__likelihood_l_pdaf, py__prepoststep_pdaf,
    int  screen, int  subtype, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Compute product of R^(-1) with HV

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    py__likelihood_l_pdaf : Callable
        Compute likelihood

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaflknetf_update_sync(&step, &dim_p, &dim_obs_f, &dim_ens,
                                  &state_p[0], &ainv[0,0], &ens_p[0,0],
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__init_obs_pdaf,
                                  pdaf_cb.c__init_obs_l_pdaf,
                                  pdaf_cb.c__prodrinva_l_pdaf,
                                  pdaf_cb.c__init_n_domains_p_pdaf,
                                  pdaf_cb.c__init_dim_l_pdaf,
                                  pdaf_cb.c__init_dim_obs_l_pdaf,
                                  pdaf_cb.c__g2l_state_pdaf,
                                  pdaf_cb.c__l2g_state_pdaf,
                                  pdaf_cb.c__g2l_obs_pdaf,
                                  pdaf_cb.c__init_obsvar_pdaf,
                                  pdaf_cb.c__init_obsvar_l_pdaf,
                                  pdaf_cb.c__likelihood_l_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf, &screen,
                                  &subtype, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, flag


def etkf_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1,:] ens_p, double [::1,:] hz_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  type_trans, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    type_trans : int
        Type of ensemble transformation
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        on exit: weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv = ainv_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_p_np = np.asarray(hz_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_etkf_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &state_p[0],
                         &ainv[0,0], &ens_p[0,0], &hz_p[0,0], &hxbar_p[0],
                         &obs_p[0], &forget, pdaf_cb.c__prodrinva_pdaf,
                         &screen, &type_trans, &debug, &flag)

    return state_p_np, ainv_np, ens_p_np, hz_p_np, flag


def letkf_ana(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, double [::1] state_l, double [::1,:] ens_l,
    double [::1,:] hz_l, double [::1] hxbar_l, double [::1] obs_l,
    double [::1,:] rndmat, double  forget, py__prodrinva_l_pdaf,
    int  type_trans, int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    type_trans : int
        Type of ensemble transformation
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        on exit: local weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    rndmat : ndarray[np.float64, ndim=2]
        Global random rotation matrix
        Array shape: (dim_ens, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv_l = ainv_l_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_l_np = np.asarray(hz_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] rndmat_np = np.asarray(rndmat, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_letkf_ana(&domain_p, &step, &dim_l, &dim_obs_l, &dim_ens,
                          &state_l[0], &ainv_l[0,0], &ens_l[0,0],
                          &hz_l[0,0], &hxbar_l[0], &obs_l[0], &rndmat[0,0],
                          &forget, pdaf_cb.c__prodrinva_l_pdaf,
                          &type_trans, &screen, &debug, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hz_l_np, rndmat_np, forget, flag


def letkf_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_letkf_init(&subtype, &param_int[0], &dim_pint,
                           &param_real[0], &dim_preal, &ensemblefilter,
                           &fixedbasis, &verbose, &outflag)

    return subtype, param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def letkf_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_letkf_alloc(&outflag)

    return outflag


def letkf_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_letkf_config(&subtype, &verbose)

    return subtype


def letkf_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_letkf_set_iparam(&id, &value, &flag)

    return flag


def letkf_set_rparam(int  id, double  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : double
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_letkf_set_rparam(&id, &value, &flag)

    return flag


def letkf_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_letkf_options()



def letkf_memtime(int  printtype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    printtype : int
        Type of screen output:

    Returns
    -------
    """
    with nogil:
        c__pdaf_letkf_memtime(&printtype)



def estkf_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    int  rank, double [::1] state_p, double [::1,:] ainv,
    double [::1,:] ens_p, double [::1,:] hl_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  envar_mode, int  type_sqrt, int  type_trans, double [::1,:] ta,
    int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast mean state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix A - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 with some matrix

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    envar_mode : int
        Flag whether routine is called from 3DVar for special functionality
    type_sqrt : int
        Type of square-root of A
    type_trans : int
        Type of ensemble transformation
    ta : ndarray[np.float64, ndim=2]
        Ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast mean state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix A - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    ta : ndarray[np.float64, ndim=2]
        Ensemble transformation matrix
        Array shape: (dim_ens, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_p_np = np.asarray(hl_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ta_np = np.asarray(ta, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_estkf_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                          &state_p[0], &ainv[0,0], &ens_p[0,0], &hl_p[0,0],
                          &hxbar_p[0], &obs_p[0], &forget,
                          pdaf_cb.c__prodrinva_pdaf, &screen, &envar_mode,
                          &type_sqrt, &type_trans, &ta[0,0], &debug, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, hl_p_np, ta_np, flag


def ensrf_ana(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ens_p, double [::1,:] hx_p,
    double [::1] hxbar_p, double [::1] obs_p, double [::1] var_obs_p,
    py__localize_covar_serial_pdaf, int  screen, int  debug):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of state ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    var_obs_p : ndarray[np.float64, ndim=1]
        PE-local vector of observation eror variances
        Array shape: (dim_obs_p)
    py__localize_covar_serial_pdaf : Callable
        Apply localization for single-observation vectors

        Callback Parameters
        -------------------
        iobs : int
                Index of current observation
        dim_p : int
                Process-local state dimension
        dim_obs : int
                Number of observations
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

        Callback Returns
        ----------------
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hx_p_np = np.asarray(hx_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxbar_p_np = np.asarray(hxbar_p, dtype=np.float64, order="F")
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    with nogil:
        c__pdaf_ensrf_ana(&step, &dim_p, &dim_obs_p, &dim_ens, &state_p[0],
                          &ens_p[0,0], &hx_p[0,0], &hxbar_p[0], &obs_p[0],
                          &var_obs_p[0],
                          pdaf_cb.c__localize_covar_serial_pdaf, &screen,
                          &debug)

    return state_p_np, ens_p_np, hx_p_np, hxbar_p_np


def ensrf_ana_2step(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ens_p, double [::1,:] hx_p,
    double [::1] hxbar_p, double [::1] obs_p, double [::1] var_obs_p,
    py__localize_covar_serial_pdaf, int  screen, int  debug):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of state ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    var_obs_p : ndarray[np.float64, ndim=1]
        PE-local vector of observation eror variances
        Array shape: (dim_obs_p)
    py__localize_covar_serial_pdaf : Callable
        Apply localization for single-observation vectors

        Callback Parameters
        -------------------
        iobs : int
                Index of current observation
        dim_p : int
                Process-local state dimension
        dim_obs : int
                Number of observations
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

        Callback Returns
        ----------------
        hp_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HP for observation iobs
                Array shape: (dim_p)
        hxy_p : ndarray[np.float64, ndim=1]
                Process-local part of matrix HX(HX_all) for full observations
                Array shape: (dim_obs)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local ensemble mean state
        Array shape: (dim_p)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hx_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hx_p_np = np.asarray(hx_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxbar_p_np = np.asarray(hxbar_p, dtype=np.float64, order="F")
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    with nogil:
        c__pdaf_ensrf_ana_2step(&step, &dim_p, &dim_obs_p, &dim_ens,
                                &state_p[0], &ens_p[0,0], &hx_p[0,0],
                                &hxbar_p[0], &obs_p[0], &var_obs_p[0],
                                pdaf_cb.c__localize_covar_serial_pdaf,
                                &screen, &debug)

    return state_p_np, ens_p_np, hx_p_np, hxbar_p_np


def lnetf_update(int  step, int  dim_p, int  dim_ens,
    double [::1] state_p, double [::1,:] ainv, double [::1,:] ens_p,
    py__obs_op_pdaf, py__init_dim_obs_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__likelihood_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__prepoststep_pdaf,
    int  screen, int  subtype, int  dim_lag, double [::1,:,:] sens_p,
    int  cnt_maxlag, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__init_obs_pdaf : Callable
        Initialize PE-local observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__likelihood_l_pdaf : Callable
        Compute observation likelihood for an ensemble member

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)

        Callback Returns
        ----------------
        resid_l : ndarray[np.float64, ndim=1]
                nput vector holding the local residual
                Array shape: (dim_obs_l)
        likely_l : double
                Output value of the local likelihood

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from global state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

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

    screen : int
        Verbosity flag
    subtype : int
        Filter subtype
    dim_lag : int
        Status flag
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag

    Returns
    -------
    dim_obs_f : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Inverse of matrix U
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    dim_lag : int
        Status flag
    sens_p : ndarray[np.float64, ndim=3]
        PE-local smoother ensemble
        Array shape: (dim_p, dim_ens, dim_lag)
    cnt_maxlag : int
        Count number of past time steps for smoothing
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=3, mode="fortran", negative_indices=False, cast=False] sens_p_np = np.asarray(sens_p, dtype=np.float64, order="F")
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  dim_obs_f
    with nogil:
        c__pdaflnetf_update(&step, &dim_p, &dim_obs_f, &dim_ens,
                            &state_p[0], &ainv[0,0], &ens_p[0,0],
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__init_obs_pdaf,
                            pdaf_cb.c__init_obs_l_pdaf,
                            pdaf_cb.c__likelihood_l_pdaf,
                            pdaf_cb.c__init_n_domains_p_pdaf,
                            pdaf_cb.c__init_dim_l_pdaf,
                            pdaf_cb.c__init_dim_obs_l_pdaf,
                            pdaf_cb.c__g2l_state_pdaf,
                            pdaf_cb.c__l2g_state_pdaf,
                            pdaf_cb.c__g2l_obs_pdaf,
                            pdaf_cb.c__prepoststep_pdaf, &screen, &subtype,
                            &dim_lag, &sens_p[0,0,0], &cnt_maxlag, &flag)

    return dim_obs_f, state_p_np, ainv_np, ens_p_np, dim_lag, sens_p_np, cnt_maxlag, flag


def seik_ana_trans(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    int  rank, double [::1] state_p, double [::1,:] uinv,
    double [::1,:] ens_p, double [::1,:] hl_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  type_sqrt, int  type_trans, int  nm1vsn, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_p : ndarray[np.float64, ndim=1]
        PE-local forecast mean state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    type_sqrt : int
        Type of square-root of A
    type_trans : int
        Type of ensemble transformation
    nm1vsn : int
        Type of normalization in covariance matrix computation
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local forecast mean state
        Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hl_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble (perturbations)
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] uinv_np = np.asarray(uinv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_p_np = np.asarray(hl_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_seik_ana_trans(&step, &dim_p, &dim_obs_p, &dim_ens, &rank,
                               &state_p[0], &uinv[0,0], &ens_p[0,0],
                               &hl_p[0,0], &hxbar_p[0], &obs_p[0], &forget,
                               pdaf_cb.c__prodrinva_pdaf, &screen,
                               &type_sqrt, &type_trans, &nm1vsn, &debug, &flag)

    return state_p_np, uinv_np, ens_p_np, hl_p_np, flag


def hyb3dvar_update_estkf(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec, int  dim_cvec_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__prodrinva_pdaf,
    py__prepoststep_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_obsvar_pdaf, int  screen,
    int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_cvec : int
        Size of control vector (parameterized part)
    dim_cvec_ens : int
        Size of control vector (ensemble part)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for 3DVAR analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

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

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble var)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafhyb3dvar_update_estkf(&step, &dim_p, &dim_obs_p, &dim_ens,
                                     &dim_cvec, &dim_cvec_ens, &state_p[0],
                                     &ainv[0,0], &ens_p[0,0],
                                     pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__init_obs_pdaf,
                                     pdaf_cb.c__prodrinva_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     pdaf_cb.c__init_obsvar_pdaf, &screen,
                                     &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag


def hyb3dvar_update_lestkf(int  step, int  dim_p, int  dim_ens,
    int  dim_cvec, int  dim_cvec_ens, double [::1] state_p,
    double [::1,:] ainv, double [::1,:] ens_p, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_obs_pdaf, py__prodrinva_pdaf,
    py__prepoststep_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf,
    py__init_obs_f_pdaf, py__init_obs_l_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    int  screen, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_ens : int
        Size of ensemble
    dim_cvec : int
        Size of control vector (parameterized part)
    dim_cvec_ens : int
        Size of control vector (ensemble part)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix for LESKTF
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector
                (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of PE-local observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector
                (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_pdaf : Callable
        Initialize observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of the observation vector

        Callback Returns
        ----------------
        observation_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A for 3DVAR analysis

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

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

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix (ensemble)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local dimension of state
        dim_ens : int
                Ensemble size
        dim_cvec_ens : int
                Dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        v_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec_ens)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local state increment
                Array shape: (dim_p)

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix (ensemble var)

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_ens : int
                Ensemble size
        dim_cv_ens_p : int
                PE-local dimension of control vector
        ens_p : ndarray[np.float64, ndim=2]
                PE-local ensemble
                Array shape: (dim_p, dim_ens)
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local input vector
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local result vector
                Array shape: (dim_cv_ens_p)

    py__cvt_pdaf : Callable
        Apply control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

        Callback Returns
        ----------------
        vv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

        Callback Parameters
        -------------------
        iter : int
                Iteration of optimization
        dim_p : int
                PE-local observation dimension
        dim_cvec : int
                Dimension of control vector
        vcv_p : ndarray[np.float64, ndim=1]
                PE-local result vector (state vector increment)
                Array shape: (dim_p)
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

        Callback Returns
        ----------------
        cv_p : ndarray[np.float64, ndim=1]
                PE-local control vector
                Array shape: (dim_cvec)

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                PE-local dimension of state
        dim_obs_p : int
                Dimension of observed state
        m_state_p : ndarray[np.float64, ndim=1]
                PE-local observed state
                Array shape: (dim_obs_p)
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                PE-local model state
                Array shape: (dim_p)

    py__init_dim_obs_f_pdaf : Callable
        Initialize dimension of observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_f_pdaf : Callable
        Observation operator

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_p : int
                Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p : int
                Size of observation vector
        state_p : ndarray[np.float64, ndim=1]
                Model state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        m_state_p : ndarray[np.float64, ndim=1]
                Observed state vector (i.e. the result after applying the observation operator to state_p)
                Array shape: (dim_obs_p)

    py__init_obs_f_pdaf : Callable
        Initialize PE-local observation vector

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_f : int
                Size of the full observation vector

        Callback Returns
        ----------------
        observation_f : ndarray[np.float64, ndim=1]
                Full vector of observations
                Array shape: (dim_obs_f)

    py__init_obs_l_pdaf : Callable
        Init. observation vector on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local size of the observation vector

        Callback Returns
        ----------------
        observation_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A on local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        n_domains_p : int
                pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain

        Callback Returns
        ----------------
        dim_l : int
                local state dimension

    py__init_dim_obs_l_pdaf : Callable
        Initialize dim. of obs. vector for local ana. domain

        Callback Parameters
        -------------------
        domain_p : int
                index of current local analysis domain
        step : int
                current time step
        dim_obs_f : int
                full dimension of observation vector

        Callback Returns
        ----------------
        dim_obs_l : int
                local dimension of observation vector

    py__g2l_state_pdaf : Callable
        Get state on local ana. domain from full state

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)
        dim_l : int
                local state dimension

        Callback Returns
        ----------------
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)

    py__l2g_state_pdaf : Callable
        Init full state from state on local analysis domain

        Callback Parameters
        -------------------
        step : int
                current time step
        domain_p : int
                current local analysis domain
        dim_l : int
                local state dimension
        state_l : ndarray[np.float64, ndim=1]
                state vector on local analysis domain
                Array shape: (dim_l)
        dim_p : int
                pe-local full state dimension
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

        Callback Returns
        ----------------
        state_p : ndarray[np.float64, ndim=1]
                pe-local full state vector
                Array shape: (dim_p)

    py__g2l_obs_pdaf : Callable
        Restrict full obs. vector to local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_f : int
                Size of full observation vector for model sub-domain
        dim_obs_l : int
                Size of observation vector for local analysis domain
        mstate_f : ndarray[np.intc, ndim=1]
                Full observation vector for model sub-domain
                Array shape: (dim_p)
        dim_p : int
                Size of full observation vector for model sub-domain
        dim_l : int
                Size of observation vector for local analysis domain

        Callback Returns
        ----------------
        mstate_l : ndarray[np.intc, ndim=1]
                Observation vector for local analysis domain
                Array shape: (dim_l)

    py__init_obsvar_pdaf : Callable
        Initialize mean observation error variance

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Size of observation vector
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)

        Callback Returns
        ----------------
        meanvar : double
                Mean observation error variance

    py__init_obsvar_l_pdaf : Callable
        Initialize local mean observation error variance

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Local dimension of observation vector
        obs_l : ndarray[np.float64, ndim=1]
                Local observation vector
                Array shape: (dim_obs_p)
        dim_obs_p : int
                Dimension of local observation vector

        Callback Returns
        ----------------
        meanvar_l : double
                Mean local observation error variance

    screen : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    dim_obs_p : int
        PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        Transform matrix for LESKTF
        Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local ensemble matrix
        Array shape: (dim_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.asarray(ainv, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    pdaf_cb.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    pdaf_cb.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef int  dim_obs_p
    with nogil:
        c__pdafhyb3dvar_update_lestkf(&step, &dim_p, &dim_obs_p, &dim_ens,
                                      &dim_cvec, &dim_cvec_ens,
                                      &state_p[0], &ainv[0,0], &ens_p[0,0],
                                      pdaf_cb.c__init_dim_obs_pdaf,
                                      pdaf_cb.c__obs_op_pdaf,
                                      pdaf_cb.c__init_obs_pdaf,
                                      pdaf_cb.c__prodrinva_pdaf,
                                      pdaf_cb.c__prepoststep_pdaf,
                                      pdaf_cb.c__cvt_ens_pdaf,
                                      pdaf_cb.c__cvt_adj_ens_pdaf,
                                      pdaf_cb.c__cvt_pdaf,
                                      pdaf_cb.c__cvt_adj_pdaf,
                                      pdaf_cb.c__obs_op_lin_pdaf,
                                      pdaf_cb.c__obs_op_adj_pdaf,
                                      pdaf_cb.c__init_dim_obs_f_pdaf,
                                      pdaf_cb.c__obs_op_f_pdaf,
                                      pdaf_cb.c__init_obs_f_pdaf,
                                      pdaf_cb.c__init_obs_l_pdaf,
                                      pdaf_cb.c__prodrinva_l_pdaf,
                                      pdaf_cb.c__init_n_domains_p_pdaf,
                                      pdaf_cb.c__init_dim_l_pdaf,
                                      pdaf_cb.c__init_dim_obs_l_pdaf,
                                      pdaf_cb.c__g2l_state_pdaf,
                                      pdaf_cb.c__l2g_state_pdaf,
                                      pdaf_cb.c__g2l_obs_pdaf,
                                      pdaf_cb.c__init_obsvar_pdaf,
                                      pdaf_cb.c__init_obsvar_l_pdaf,
                                      &screen, &flag)

    return dim_obs_p, state_p_np, ainv_np, ens_p_np, flag

def lestkf_ana_fixed(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, int  rank, double [::1] state_l, double [::1,:] ainv_l,
    double [::1,:] ens_l, double [::1,:] hl_l, double [::1] hxbar_l,
    double [::1] obs_l, double  forget, py__prodrinva_l_pdaf,
    int  type_sqrt, int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    rank : int
        Rank of initial covariance matrix
    state_l : ndarray[np.float64, ndim=1]
        state on local analysis domain
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    type_sqrt : int
        Type of square-root of A
    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        state on local analysis domain
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        Inverse of matrix U - temporary use only
        Array shape: (rank, rank)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hl_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.asarray(ainv_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hl_l_np = np.asarray(hl_l, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_lestkf_ana_fixed(&domain_p, &step, &dim_l, &dim_obs_l,
                                 &dim_ens, &rank, &state_l[0],
                                 &ainv_l[0,0], &ens_l[0,0], &hl_l[0,0],
                                 &hxbar_l[0], &obs_l[0], &forget,
                                 pdaf_cb.c__prodrinva_l_pdaf, &type_sqrt,
                                 &screen, &debug, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hl_l_np, forget, flag


def genobs_init(int  subtype, int [::1] param_int, int  dim_pint,
    double [::1] param_real, int  dim_preal, int  verbose, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    dim_pint : int
        Number of integer parameters
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    dim_preal : int
        Number of real parameters
    verbose : int
        Control screen output
    outflag : int
        Status flag

    Returns
    -------
    param_int : ndarray[np.intc, ndim=1]
        Integer parameter array
        Array shape: (dim_pint)
    param_real : ndarray[np.float64, ndim=1]
        Real parameter array
        Array shape: (dim_preal)
    ensemblefilter : bint
        Is the chosen filter ensemble-based?
    fixedbasis : bint
        Does the filter run with fixed error-space basis?
    outflag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_int_np = np.asarray(param_int, dtype=np.intc, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] param_real_np = np.asarray(param_real, dtype=np.float64, order="F")
    cdef bint  ensemblefilter
    cdef bint  fixedbasis
    with nogil:
        c__pdaf_genobs_init(&subtype, &param_int[0], &dim_pint,
                            &param_real[0], &dim_preal, &ensemblefilter,
                            &fixedbasis, &verbose, &outflag)

    return param_int_np, param_real_np, ensemblefilter, fixedbasis, outflag


def genobs_alloc(int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    with nogil:
        c__pdaf_genobs_alloc(&outflag)

    return outflag


def genobs_config(int  subtype, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    subtype : int
        Sub-type of filter
    verbose : int
        Control screen output

    Returns
    -------
    subtype : int
        Sub-type of filter
    """
    with nogil:
        c__pdaf_genobs_config(&subtype, &verbose)

    return subtype


def genobs_set_iparam(int  id, int  value):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    id : int
        Index of parameter
    value : int
        Parameter value

    Returns
    -------
    flag : int
        Status flag: 0 for no error
    """
    cdef int  flag
    with nogil:
        c__pdaf_genobs_set_iparam(&id, &value, &flag)

    return flag


def genobs_options():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdaf_genobs_options()



def etkf_ana_t(int  step, int  dim_p, int  dim_obs_p, int  dim_ens,
    double [::1,:] ens_p, double [::1,:] hz_p, double [::1] hxbar_p,
    double [::1] obs_p, double  forget, py__prodrinva_pdaf, int  screen,
    int  type_trans, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_p : int
        PE-local dimension of model state
    dim_obs_p : int
        PE-local dimension of observation vector
    dim_ens : int
        Size of ensemble
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    hxbar_p : ndarray[np.float64, ndim=1]
        PE-local observed state
        Array shape: (dim_obs_p)
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)
    forget : double
        Forgetting factor
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

        Callback Parameters
        -------------------
        step : int
                Current time step
        dim_obs_p : int
                Number of observations at current time step (i.e. the size of the observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one
                (or the rank of the initial covariance matrix)
        obs_p : ndarray[np.float64, ndim=1]
                Vector of observations
                Array shape: (dim_obs_p)
        a_p : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_p, rank)

        Callback Returns
        ----------------
        c_p : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_p, rank)

    screen : int
        Verbosity flag
    type_trans : int
        Type of ensemble transformation
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        on exit: PE-local forecast state
        Array shape: (dim_p)
    ainv : ndarray[np.float64, ndim=2]
        on exit: weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_p : ndarray[np.float64, ndim=2]
        PE-local state ensemble
        Array shape: (dim_p, dim_ens)
    hz_p : ndarray[np.float64, ndim=2]
        PE-local observed ensemble
        Array shape: (dim_obs_p, dim_ens)
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.zeros((dim_p), dtype=np.float64, order="F")
    cdef double [::1] state_p = state_p_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv = ainv_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_p_np = np.asarray(ens_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_p_np = np.asarray(hz_p, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    with nogil:
        c__pdaf_etkf_ana_t(&step, &dim_p, &dim_obs_p, &dim_ens,
                           &state_p[0], &ainv[0,0], &ens_p[0,0],
                           &hz_p[0,0], &hxbar_p[0], &obs_p[0], &forget,
                           pdaf_cb.c__prodrinva_pdaf, &screen, &type_trans,
                           &debug, &flag)

    return state_p_np, ainv_np, ens_p_np, hz_p_np, flag


def letkf_ana_fixed(int  domain_p, int  step, int  dim_l, int  dim_obs_l,
    int  dim_ens, double [::1] state_l, double [::1,:] ens_l,
    double [::1,:] hz_l, double [::1] hxbar_l, double [::1] obs_l,
    double  forget, py__prodrinva_l_pdaf, int  screen, int  debug, int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_l : int
        State dimension on local analysis domain
    dim_obs_l : int
        Size of obs. vector on local ana. domain
    dim_ens : int
        Size of ensemble
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    hxbar_l : ndarray[np.float64, ndim=1]
        Local observed ensemble mean
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    forget : double
        Forgetting factor
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A for local analysis domain

        Callback Parameters
        -------------------
        domain_p : int
                Index of current local analysis domain
        step : int
                Current time step
        dim_obs_l : int
                Number of local observations at current time step (i.e. the size of the local observation vector)
        rank : int
                Number of the columns in the matrix processes here.
                This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[np.float64, ndim=1]
                Local vector of observations
                Array shape: (dim_obs_l)
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)

        Callback Returns
        ----------------
        a_l : ndarray[np.float64, ndim=2]
                Input matrix provided by PDAF
                Array shape: (dim_obs_l, rank)
        c_l : ndarray[np.float64, ndim=2]
                Output matrix
                Array shape: (dim_obs_l, rank)

    screen : int
        Verbosity flag
    debug : int
        Flag for writing debug output
    flag : int
        Status flag

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        Local forecast state
        Array shape: (dim_l)
    ainv_l : ndarray[np.float64, ndim=2]
        on exit: local weight matrix for ensemble transformation
        Array shape: (dim_ens, dim_ens)
    ens_l : ndarray[np.float64, ndim=2]
        Local state ensemble
        Array shape: (dim_l, dim_ens)
    hz_l : ndarray[np.float64, ndim=2]
        Local observed state ensemble (perturbation)
        Array shape: (dim_obs_l, dim_ens)
    forget : double
        Forgetting factor
    flag : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.asarray(state_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ainv_l_np = np.zeros((dim_ens, dim_ens), dtype=np.float64, order="F")
    cdef double [::1,:] ainv_l = ainv_l_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ens_l_np = np.asarray(ens_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hz_l_np = np.asarray(hz_l, dtype=np.float64, order="F")
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    with nogil:
        c__pdaf_letkf_ana_fixed(&domain_p, &step, &dim_l, &dim_obs_l,
                                &dim_ens, &state_l[0], &ainv_l[0,0],
                                &ens_l[0,0], &hz_l[0,0], &hxbar_l[0],
                                &obs_l[0], &forget,
                                pdaf_cb.c__prodrinva_l_pdaf, &screen,
                                &debug, &flag)

    return state_l_np, ainv_l_np, ens_l_np, hz_l_np, forget, flag


