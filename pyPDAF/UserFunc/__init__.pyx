import numpy as np
import mpi4py.MPI as MPI
import sys

from traceback import print_exception
# Global error handler
def global_except_hook(exctype, value, traceback):
    
    try:
        
        sys.stderr.write("\n*****************************************************\n")
        sys.stderr.write("Uncaught exception was detected on rank {}. \n".format(
            MPI.COMM_WORLD.Get_rank()))
        
        print_exception(exctype, value, traceback)
        sys.stderr.write("*****************************************************\n\n\n")
        sys.stderr.write("\n")
        sys.stderr.write("Calling MPI_Abort() to shut down MPI processes...\n")
        sys.stderr.flush()
    finally:
        try:
            MPI.COMM_WORLD.Abort(1)
        except Exception as e:
            sys.stderr.write("*****************************************************\n")
            sys.stderr.write("Sorry, we failed to stop MPI, this process will hang.\n")
            sys.stderr.write("*****************************************************\n")
            sys.stderr.flush()
            raise e

sys.excepthook = global_except_hook


def py__add_obs_err_pdaf(step, dim_obs_p, c_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of observation vector
    c_p : ndarray[float]
        matrix to that observation covariance r is added
        shape is (dim_obs_p,dim_obs_p)

    Returns
    -------
    c_p : ndarray[float]
        matrix to that observation covariance r is added

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__add_obs_err_pdaf is called!!!...')


def py__init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, flag):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    filtertype : int
        type of filter to initialize
    dim_p : int
        pe-local state dimension
    dim_ens : int
        size of ensemble
    state_p : ndarray[float]
        pe-local model state
        shape is (dim_p)
    uinv : ndarray[float]
        array not referenced for ensemble filters
        shape is (dim_ens-1,dim_ens-1)
    ens_p : ndarray[float]
        pe-local state ensemble
        shape is (dim_p,dim_ens)
    flag : int
        pdaf status flag

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state
    uinv : ndarray[float]
        array not referenced for ensemble filters
    ens_p : ndarray[float]
        pe-local state ensemble
    flag : int
        pdaf status flag

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_ens_pdaf is called!!!...')


def py__next_observation_pdaf(stepnow, nsteps, doexit, time):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    stepnow : int
        number of the current time step
    nsteps : int
        number of time steps until next obs
    doexit : int
        whether to exit forecasting (1 for exit)
    time : float
        current model (physical) time

    Returns
    -------
    nsteps : int
        number of time steps until next obs
    doexit : int
        whether to exit forecasting (1 for exit)
    time : float
        current model (physical) time

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__next_observation_pdaf is called!!!...')


def py__collect_state_pdaf(dim_p, state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_p : int
        pe-local state dimension
    state_p : ndarray[float]
        local state vector
        shape is (dim_p)

    Returns
    -------
    state_p : ndarray[float]
        local state vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__collect_state_pdaf is called!!!...')


def py__distribute_state_pdaf(dim_p, state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_p : int
        
    state_p : ndarray[float]
        
        shape is (dim_p)

    Returns
    -------
    state_p : ndarray[float]
        

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__distribute_state_pdaf is called!!!...')


def py__prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, state_p, uinv, ens_p, flag):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step (negative for call after forecast)
    dim_p : int
        pe-local state dimension
    dim_ens : int
        size of state ensemble
    dim_ens_p : int
        pe-local size of ensemble
    dim_obs_p : int
        pe-local dimension of observation vector
    state_p : ndarray[float]
        pe-local forecast/analysis state
         (the array 'state_p' is not generally not
         initialized in the case of seik.
         it can be used freely here.)
        shape is (dim_p)
    uinv : ndarray[float]
        inverse of matrix u
        shape is (dim_ens-1,dim_ens-1)
    ens_p : ndarray[float]
        pe-local state ensemble
        shape is (dim_p,dim_ens)
    flag : int
        pdaf status flag

    Returns
    -------
    state_p : ndarray[float]
        pe-local forecast/analysis state
         (the array 'state_p' is not generally not
         initialized in the case of seik.
         it can be used freely here.)
    uinv : ndarray[float]
        inverse of matrix u
    ens_p : ndarray[float]
        pe-local state ensemble

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__prepoststep_pdaf is called!!!...')


def py__init_dim_obs_pdaf(step, dim_obs_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of observation vector

    Returns
    -------
    dim_obs_p : int
        dimension of observation vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_dim_obs_pdaf is called!!!...')


def py__init_obs_pdaf(step, dim_obs_p, observation_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        size of the observation vector
    observation_p : ndarray[float]
        vector of observations
        shape is (dim_obs_p)

    Returns
    -------
    observation_p : ndarray[float]
        vector of observations

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obs_pdaf is called!!!...')


def py__init_obs_covar_pdaf(step, dim_obs, dim_obs_p, covar, obs_p, isdiag):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        global size of observation vector
    dim_obs_p : int
        size of process-local observation vector
    covar : float
        observation error covariance matrix
    obs_p : ndarray[float]
        process-local vector of observations
        shape is (dim_obs_p)
    isdiag : bool
        

    Returns
    -------
    covar : float
        observation error covariance matrix
    isdiag : bool
        

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obs_covar_pdaf is called!!!...')


def py__init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        size of observation vector
    obs_p : ndarray[float]
        vector of observations
        shape is (dim_obs_p)
    meanvar : float
        mean observation error variance

    Returns
    -------
    meanvar : float
        mean observation error variance

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obsvar_pdaf is called!!!...')


def py__prodrinva_pdaf(step, dim_obs_p, rank, obs_p, a_p, c_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        number of observations at current time step (i.e. the size of the observation vector)
    rank : int
        number of the columns in the matrix processes here.
         this is usually the ensemble size minus one
         (or the rank of the initial covariance matrix)
    obs_p : ndarray[float]
        vector of observations
        shape is (dim_obs_p)
    a_p : ndarray[float]
        input matrix provided by pdaf
        shape is (dim_obs_p,rank)
    c_p : ndarray[float]
        output matrix
        shape is (dim_obs_p,rank)

    Returns
    -------
    c_p : ndarray[float]
        output matrix

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__prodrinva_pdaf is called!!!...')


def py__obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_p : int
        size of state vector (local part in case of parallel decomposed state)
    dim_obs_p : int
        size of observation vector
    state_p : ndarray[float]
        model state vector
        shape is (dim_p)
    m_state_p : ndarray[float]
        observed state vector (i.e. the result after applying the observation operator to state_p)
        shape is (dim_obs_p)

    Returns
    -------
    m_state_p : ndarray[float]
        observed state vector (i.e. the result after applying the observation operator to state_p)

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__obs_op_pdaf is called!!!...')


def py__g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, dim_p, mstate_l, dim_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_f : int
        size of full observation vector for model sub-domain
    dim_obs_l : int
        size of observation vector for local analysis domain
    mstate_f : ndarray[int]
        full observation vector for model sub-domain
        shape is (dim_p)
    dim_p : int
        size of full observation vector for model sub-domain
    mstate_l : ndarray[int]
        observation vector for local analysis domain
        shape is (dim_l)
    dim_l : int
        size of observation vector for local analysis domain

    Returns
    -------
    mstate_l : ndarray[int]
        observation vector for local analysis domain

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__g2l_obs_pdaf is called!!!...')


def py__g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    domain_p : int
        current local analysis domain
    dim_p : int
        pe-local full state dimension
    state_p : ndarray[float]
        pe-local full state vector
        shape is (dim_p)
    dim_l : int
        local state dimension
    state_l : ndarray[float]
        state vector on local analysis domain
        shape is (dim_l)

    Returns
    -------
    state_l : ndarray[float]
        state vector on local analysis domain

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__g2l_state_pdaf is called!!!...')


def py__init_dim_l_pdaf(step, domain_p, dim_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    domain_p : int
        current local analysis domain
    dim_l : int
        local state dimension

    Returns
    -------
    dim_l : int
        local state dimension

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_dim_l_pdaf is called!!!...')


def py__init_dim_obs_f_pdaf(step, dim_obs_f):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        size of the full observation vector

    Returns
    -------
    dim_obs_f : int
        size of the full observation vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_dim_obs_f_pdaf is called!!!...')


def py__init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_f : int
        full dimension of observation vector
    dim_obs_l : int
        local dimension of observation vector

    Returns
    -------
    dim_obs_l : int
        local dimension of observation vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_dim_obs_l_pdaf is called!!!...')


def py__init_n_domains_p_pdaf(step, n_domains_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    n_domains_p : int
        pe-local number of analysis domains

    Returns
    -------
    n_domains_p : int
        pe-local number of analysis domains

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_n_domains_p_pdaf is called!!!...')


def py__init_obs_f_pdaf(step, dim_obs_f, observation_f):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        size of the full observation vector
    observation_f : ndarray[float]
        full vector of observations
        shape is (dim_obs_f)

    Returns
    -------
    observation_f : ndarray[float]
        full vector of observations

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obs_f_pdaf is called!!!...')


def py__init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        local size of the observation vector
    observation_l : ndarray[float]
        local vector of observations
        shape is (dim_obs_l)

    Returns
    -------
    observation_l : ndarray[float]
        local vector of observations

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obs_l_pdaf is called!!!...')


def py__init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, dim_obs_p, meanvar_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        local dimension of observation vector
    obs_l : ndarray[float]
        local observation vector
        shape is (dim_obs_p)
    dim_obs_p : int
        dimension of local observation vector
    meanvar_l : float
        mean local observation error variance

    Returns
    -------
    meanvar_l : float
        mean local observation error variance

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obsvar_l_pdaf is called!!!...')


def py__init_obserr_f_pdaf(step, dim_obs_f, obs_f, obserr_f):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        full dimension of observation vector
    obs_f : ndarray[float]
        full observation vector
        shape is (dim_obs_f)
    obserr_f : ndarray[float]
        full observation error stddev
        shape is (dim_obs_f)

    Returns
    -------
    obserr_f : ndarray[float]
        full observation error stddev

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__init_obserr_f_pdaf is called!!!...')


def py__l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    domain_p : int
        current local analysis domain
    dim_l : int
        local state dimension
    state_l : ndarray[float]
        state vector on local analysis domain
        shape is (dim_l)
    dim_p : int
        pe-local full state dimension
    state_p : ndarray[float]
        pe-local full state vector
        shape is (dim_p)

    Returns
    -------
    state_p : ndarray[float]
        pe-local full state vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__l2g_state_pdaf is called!!!...')


def py__obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_p : int
        size of state vector (local part in case of parallel decomposed state)
    dim_obs_f : int
        size of full observation vector
    state_p : ndarray[float]
        model state vector
        shape is (dim_p)
    m_state_f : ndarray[float]
        full observed state (i.e. the result after applying the observation operator to state_p)
        shape is (dim_obs_f)

    Returns
    -------
    m_state_f : ndarray[float]
        full observed state (i.e. the result after applying the observation operator to state_p)

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__obs_op_f_pdaf is called!!!...')


def py__prodrinva_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, a_l, c_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        number of local observations at current time step (i.e. the size of the local observation vector)
    rank : int
        number of the columns in the matrix processes here.
         this is usually the ensemble size minus one (or the rank of the initial covariance matrix)
    obs_l : ndarray[float]
        local vector of observations
        shape is (dim_obs_l)
    a_l : ndarray[float]
        input matrix provided by pdaf
        shape is (dim_obs_l,rank)
    c_l : ndarray[float]
        output matrix
        shape is (dim_obs_l,rank)

    Returns
    -------
    c_l : ndarray[float]
        output matrix

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__prodrinva_l_pdaf is called!!!...')


def py__localize_covar_pdaf(dim_p, dim_obs, hp_p, hph):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_p : int
        pe-local state dimension
    dim_obs : int
        number of observations
    hp_p : ndarray[float]
        pe local part of matrix hp
        shape is (dim_obs,dim_p)
    hph : ndarray[float]
        matrix hph
        shape is (dim_obs,dim_obs)

    Returns
    -------
    hp_p : ndarray[float]
        pe local part of matrix hp
    hph : ndarray[float]
        matrix hph

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__localize_covar_pdaf is called!!!...')


def py__likelihood_pdaf(step, dim_obs_p, obs_p, resid, likely):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        number of observations at current time step (i.e. the size of the observation vector)
    obs_p : ndarray[float]
        vector of observations
        shape is (dim_obs_p)
    resid : ndarray[float]
        input vector holding the residual
        shape is (dim_obs_p)
    likely : float
        output value of the likelihood

    Returns
    -------
    likely : float
        output value of the likelihood

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__likelihood_pdaf is called!!!...')


def py__likelihood_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        number of local observations at current time step (i.e. the size of the local observation vector)
    obs_l : ndarray[float]
        local vector of observations
        shape is (dim_obs_l)
    resid_l : ndarray[float]
        nput vector holding the local residual
        shape is (dim_obs_l)
    likely_l : float
        output value of the local likelihood

    Returns
    -------
    likely_l : float
        output value of the local likelihood

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__likelihood_l_pdaf is called!!!...')


def py__get_obs_f_pdaf(step, dim_obs_f, observation_f):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        size of the full observation vector
    observation_f : ndarray[float]
        full vector of synthetic observations (process-local)
        shape is (dim_obs_f)

    Returns
    -------
    observation_f : ndarray[float]
        full vector of synthetic observations (process-local)

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__get_obs_f_pdaf is called!!!...')


def py__cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, vcv_p, cv_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    iter : int
        iteration of optimization
    dim_p : int
        pe-local observation dimension
    dim_ens : int
        ensemble size
    dim_cv_ens_p : int
        pe-local dimension of control vector
    ens_p : ndarray[float]
        pe-local ensemble
        shape is (dim_p,dim_ens)
    vcv_p : ndarray[float]
        pe-local input vector
        shape is (dim_p)
    cv_p : ndarray[float]
        pe-local result vector
        shape is (dim_cv_ens_p)

    Returns
    -------
    cv_p : ndarray[float]
        pe-local result vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__cvt_adj_ens_pdaf is called!!!...')


def py__cvt_adj_pdaf(iter, dim_p, dim_cvec, vcv_p, cv_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    iter : int
        iteration of optimization
    dim_p : int
        pe-local observation dimension
    dim_cvec : int
        dimension of control vector
    vcv_p : ndarray[float]
        pe-local result vector (state vector increment)
        shape is (dim_p)
    cv_p : ndarray[float]
        pe-local control vector
        shape is (dim_cvec)

    Returns
    -------
    cv_p : ndarray[float]
        pe-local control vector

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__cvt_adj_pdaf is called!!!...')


def py__cvt_pdaf(iter, dim_p, dim_cvec, cv_p, vv_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    iter : int
        iteration of optimization
    dim_p : int
        pe-local observation dimension
    dim_cvec : int
        dimension of control vector
    cv_p : ndarray[float]
        pe-local control vector
        shape is (dim_cvec)
    vv_p : ndarray[float]
        pe-local result vector (state vector increment)
        shape is (dim_p)

    Returns
    -------
    vv_p : ndarray[float]
        pe-local result vector (state vector increment)

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__cvt_pdaf is called!!!...')


def py__cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, vv_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    iter : int
        iteration of optimization
    dim_p : int
        pe-local dimension of state
    dim_ens : int
        ensemble size
    dim_cvec_ens : int
        dimension of control vector
    ens_p : ndarray[float]
        pe-local ensemble
        shape is (dim_p,dim_ens)
    v_p : ndarray[float]
        pe-local control vector
        shape is (dim_cvec_ens)
    vv_p : ndarray[float]
        pe-local state increment
        shape is (dim_p)

    Returns
    -------
    vv_p : ndarray[float]
        pe-local state increment

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__cvt_ens_pdaf is called!!!...')


def py__obs_op_adj_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_p : int
        pe-local dimension of state
    dim_obs_p : int
        dimension of observed state
    state_p : ndarray[float]
        pe-local model state
        shape is (dim_p)
    m_state_p : ndarray[float]
        pe-local observed state
        shape is (dim_obs_p)

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__obs_op_adj_pdaf is called!!!...')


def py__obs_op_lin_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_p : int
        pe-local dimension of state
    dim_obs_p : int
        dimension of observed state
    state_p : ndarray[float]
        pe-local model state
        shape is (dim_p)
    m_state_p : ndarray[float]
        pe-local observed state
        shape is (dim_obs_p)

    Returns
    -------
    m_state_p : ndarray[float]
        pe-local observed state

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__obs_op_lin_pdaf is called!!!...')


def py__dist_stateinc_pdaf(dim_p, state_inc_p, first, steps):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_p : int
        dimension of pe-local state
    state_inc_p : ndarray[float]
        pe-local state vector
        shape is (dim_p)
    first : int
        flag for first call of each forecast
    steps : int
        number of time steps in forecast
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__dist_stateinc_pdaf is called!!!...')


def py__prodrinva_hyb_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, gamma, likely_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        number of local observations at current time step (i.e. the size of the local observation vector)
    obs_l : ndarray[float]
        local vector of observations
        shape is (dim_obs_l)
    resid_l : ndarray[float]
        input vector holding the local residual
        shape is (dim_obs_l)
    gamma : float
        hybrid weight provided by pdaf
    likely_l : float
        output value of the local likelihood

    Returns
    -------
    likely_l : float
        output value of the local likelihood

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__prodrinva_hyb_l_pdaf is called!!!...')


def py__likelihood_hyb_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, gamma, a_l, c_l):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        number of local observations at current time step (i.e. the size of the local observation vector)
    rank : int
        number of the columns in the matrix processes here. this is usually the ensemble size minus one (or the rank of the initial covariance matrix)
    obs_l : ndarray[float]
        local vector of observations
        shape is (dim_obs_l)
    gamma : float
        hybrid weight provided by pdaf
    a_l : ndarray[float]
        input matrix provided by pdaf
        shape is (dim_obs_l,rank)
    c_l : ndarray[float]
        output matrix
        shape is (dim_obs_l,rank)

    Returns
    -------
    c_l : ndarray[float]
        output matrix

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong py__likelihood_hyb_l_pdaf is called!!!...')


cdef void c__add_obs_err_pdaf (int* step, int* dim_obs_p, double* c_p):
    c_p_np = np.asarray(<double[:np.prod((dim_obs_p[0], dim_obs_p[0]))]> 
                        c_p).reshape((dim_obs_p[0], dim_obs_p[0]), order='F')

    c_p_np_tmp = py__add_obs_err_pdaf(step[0], dim_obs_p[0], c_p_np)

    c_p_np[:] = c_p_np_tmp[:]
    cdef double[::1] c_p_view = c_p_np.ravel(order='F')
    assert c_p == &c_p_view[0], 'reference (memory address) of c_p has changed in c__add_obs_err_pdaf.'


cdef void c__init_ens_pdaf (int* filtertype, int* dim_p, int* dim_ens, double* state_p, double* uinv, double* ens_p, int* flag):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    if dim_ens[0] > 1:
        uinv_np = np.asarray(<double[:np.prod((dim_ens[0]-1, dim_ens[0]-1))]> 
                            uinv).reshape((dim_ens[0]-1, dim_ens[0]-1), order='F')
    else:
        uinv_np = None
    ens_p_np = np.asarray(<double[:np.prod((dim_p[0], dim_ens[0]))]> 
                          ens_p).reshape((dim_p[0], dim_ens[0]), order='F')

    state_p_np_tmp, uinv_np_tmp, ens_p_np_tmp, flag[0] = py__init_ens_pdaf(filtertype[0], dim_p[0], dim_ens[0], state_p_np, uinv_np, ens_p_np, flag[0])

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__init_ens_pdaf.'
    cdef double[::1] uinv_view
    if dim_ens[0] > 1:
        uinv_np[:] = uinv_np_tmp[:]
        uinv_view = uinv_np.ravel(order='F')
        assert uinv == &uinv_view[0], 'reference (memory address) of uinv has changed in c__init_ens_pdaf.'
    ens_p_np[:] = ens_p_np_tmp[:]
    cdef double[::1] ens_p_view = ens_p_np.ravel(order='F')
    assert ens_p == &ens_p_view[0], 'reference (memory address) of ens_p has changed in c__init_ens_pdaf.'


cdef void c__next_observation_pdaf (int* stepnow, int* nsteps, int* doexit, double* time):
    nsteps[0], doexit[0], time[0] = py__next_observation_pdaf(stepnow[0], nsteps[0], doexit[0], time[0])



cdef void c__collect_state_pdaf (int* dim_p, double* state_p):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')

    state_p_np_tmp = py__collect_state_pdaf(dim_p[0], state_p_np)

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__collect_state_pdaf.'


cdef void c__distribute_state_pdaf (int* dim_p, double* state_p):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')

    state_p_np_tmp = py__distribute_state_pdaf(dim_p[0], state_p_np)

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__distribute_state_pdaf.'


cdef void c__prepoststep_pdaf (int* step, int* dim_p, int* dim_ens, int* dim_ens_p, int* dim_obs_p, double* state_p, double* uinv, double* ens_p, int* flag):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    if dim_ens[0] > 1:
        uinv_np = np.asarray(<double[:np.prod((dim_ens[0]-1, dim_ens[0]-1))]> 
                            uinv).reshape((dim_ens[0]-1, dim_ens[0]-1), order='F')
    else:
        uinv_np = None
    ens_p_np = np.asarray(<double[:np.prod((dim_p[0], dim_ens[0]))]> 
                          ens_p).reshape((dim_p[0], dim_ens[0]), order='F')

    state_p_np_tmp, uinv_np_tmp, ens_p_np_tmp = py__prepoststep_pdaf(step[0], dim_p[0], dim_ens[0], dim_ens_p[0], dim_obs_p[0], state_p_np, uinv_np, ens_p_np, flag[0])

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__prepoststep_pdaf.'
    cdef double[::1] uinv_view
    if dim_ens[0] > 1:
        uinv_np[:] = uinv_np_tmp[:]
        uinv_view = uinv_np.ravel(order='F')
        assert uinv == &uinv_view[0], 'reference (memory address) of uinv has changed in c__prepoststep_pdaf.'
    ens_p_np[:] = ens_p_np_tmp[:]
    cdef double[::1] ens_p_view = ens_p_np.ravel(order='F')
    assert ens_p == &ens_p_view[0], 'reference (memory address) of ens_p has changed in c__prepoststep_pdaf.'


cdef void c__init_dim_obs_pdaf (int* step, int* dim_obs_p):
    dim_obs_p[0] = py__init_dim_obs_pdaf(step[0], dim_obs_p[0])



cdef void c__init_obs_pdaf (int* step, int* dim_obs_p, double* observation_p):
    observation_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                                  observation_p).reshape((dim_obs_p[0]), order='F')

    observation_p_np_tmp = py__init_obs_pdaf(step[0], dim_obs_p[0], observation_p_np)

    observation_p_np[:] = observation_p_np_tmp[:]
    cdef double[::1] observation_p_view = observation_p_np.ravel(order='F')
    assert observation_p == &observation_p_view[0], 'reference (memory address) of observation_p has changed in c__init_obs_pdaf.'


cdef void c__init_obs_covar_pdaf (int* step, int* dim_obs, int* dim_obs_p, double* covar, double* obs_p, bint* isdiag):
    obs_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          obs_p).reshape((dim_obs_p[0]), order='F')

    covar[0], isdiag[0] = py__init_obs_covar_pdaf(step[0], dim_obs[0], dim_obs_p[0], covar[0], obs_p_np, isdiag[0])



cdef void c__init_obsvar_pdaf (int* step, int* dim_obs_p, double* obs_p, double* meanvar):
    obs_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          obs_p).reshape((dim_obs_p[0]), order='F')

    meanvar[0] = py__init_obsvar_pdaf(step[0], dim_obs_p[0], obs_p_np, meanvar[0])



cdef void c__prodrinva_pdaf (int* step, int* dim_obs_p, int* rank, double* obs_p, double* a_p, double* c_p):
    obs_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          obs_p).reshape((dim_obs_p[0]), order='F')
    a_p_np = np.asarray(<double[:np.prod((dim_obs_p[0], rank[0]))]> 
                        a_p).reshape((dim_obs_p[0], rank[0]), order='F')
    c_p_np = np.asarray(<double[:np.prod((dim_obs_p[0], rank[0]))]> 
                        c_p).reshape((dim_obs_p[0], rank[0]), order='F')

    c_p_np_tmp = py__prodrinva_pdaf(step[0], dim_obs_p[0], rank[0], obs_p_np, a_p_np, c_p_np)

    c_p_np[:] = c_p_np_tmp[:]
    cdef double[::1] c_p_view = c_p_np.ravel(order='F')
    assert c_p == &c_p_view[0], 'reference (memory address) of c_p has changed in c__prodrinva_pdaf.'


cdef void c__obs_op_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    m_state_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                              m_state_p).reshape((dim_obs_p[0]), order='F')

    m_state_p_np_tmp = py__obs_op_pdaf(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)

    m_state_p_np[:] = m_state_p_np_tmp[:]
    cdef double[::1] m_state_p_view = m_state_p_np.ravel(order='F')
    assert m_state_p == &m_state_p_view[0], 'reference (memory address) of m_state_p has changed in c__obs_op_pdaf.'


cdef void c__g2l_obs_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l, int* mstate_f, int* dim_p, int* mstate_l, int* dim_l):
    mstate_f_np = np.asarray(<int[:np.prod((dim_p[0]))]> 
                             mstate_f).reshape((dim_p[0]), order='F')
    mstate_l_np = np.asarray(<int[:np.prod((dim_l[0]))]> 
                             mstate_l).reshape((dim_l[0]), order='F')

    mstate_l_np_tmp = py__g2l_obs_pdaf(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0], mstate_f_np, dim_p[0], mstate_l_np, dim_l[0])

    mstate_l_np[:] = mstate_l_np_tmp[:]
    cdef int[::1] mstate_l_view = mstate_l_np.ravel(order='F')
    assert mstate_l == &mstate_l_view[0], 'reference (memory address) of mstate_l has changed in c__g2l_obs_pdaf.'


cdef void c__g2l_state_pdaf (int* step, int* domain_p, int* dim_p, double* state_p, int* dim_l, double* state_l):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    state_l_np = np.asarray(<double[:np.prod((dim_l[0]))]> 
                            state_l).reshape((dim_l[0]), order='F')

    state_l_np_tmp = py__g2l_state_pdaf(step[0], domain_p[0], dim_p[0], state_p_np, dim_l[0], state_l_np)

    state_l_np[:] = state_l_np_tmp[:]
    cdef double[::1] state_l_view = state_l_np.ravel(order='F')
    assert state_l == &state_l_view[0], 'reference (memory address) of state_l has changed in c__g2l_state_pdaf.'


cdef void c__init_dim_l_pdaf (int* step, int* domain_p, int* dim_l):
    dim_l[0] = py__init_dim_l_pdaf(step[0], domain_p[0], dim_l[0])



cdef void c__init_dim_obs_f_pdaf (int* step, int* dim_obs_f):
    dim_obs_f[0] = py__init_dim_obs_f_pdaf(step[0], dim_obs_f[0])



cdef void c__init_dim_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l):
    dim_obs_l[0] = py__init_dim_obs_l_pdaf(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0])



cdef void c__init_n_domains_p_pdaf (int* step, int* n_domains_p):
    n_domains_p[0] = py__init_n_domains_p_pdaf(step[0], n_domains_p[0])



cdef void c__init_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f):
    observation_f_np = np.asarray(<double[:np.prod((dim_obs_f[0]))]> 
                                  observation_f).reshape((dim_obs_f[0]), order='F')

    observation_f_np_tmp = py__init_obs_f_pdaf(step[0], dim_obs_f[0], observation_f_np)

    observation_f_np[:] = observation_f_np_tmp[:]
    cdef double[::1] observation_f_view = observation_f_np.ravel(order='F')
    assert observation_f == &observation_f_view[0], 'reference (memory address) of observation_f has changed in c__init_obs_f_pdaf.'


cdef void c__init_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* observation_l):
    observation_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                                  observation_l).reshape((dim_obs_l[0]), order='F')

    observation_l_np_tmp = py__init_obs_l_pdaf(domain_p[0], step[0], dim_obs_l[0], observation_l_np)

    observation_l_np[:] = observation_l_np_tmp[:]
    cdef double[::1] observation_l_view = observation_l_np.ravel(order='F')
    assert observation_l == &observation_l_view[0], 'reference (memory address) of observation_l has changed in c__init_obs_l_pdaf.'


cdef void c__init_obsvar_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, int* dim_obs_p, double* meanvar_l):
    obs_l_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          obs_l).reshape((dim_obs_p[0]), order='F')

    meanvar_l[0] = py__init_obsvar_l_pdaf(domain_p[0], step[0], dim_obs_l[0], obs_l_np, dim_obs_p[0], meanvar_l[0])



cdef void c__init_obserr_f_pdaf (int* step, int* dim_obs_f, double* obs_f, double* obserr_f):
    obs_f_np = np.asarray(<double[:np.prod((dim_obs_f[0]))]> 
                          obs_f).reshape((dim_obs_f[0]), order='F')
    obserr_f_np = np.asarray(<double[:np.prod((dim_obs_f[0]))]> 
                             obserr_f).reshape((dim_obs_f[0]), order='F')

    obserr_f_np_tmp = py__init_obserr_f_pdaf(step[0], dim_obs_f[0], obs_f_np, obserr_f_np)

    obserr_f_np[:] = obserr_f_np_tmp[:]
    cdef double[::1] obserr_f_view = obserr_f_np.ravel(order='F')
    assert obserr_f == &obserr_f_view[0], 'reference (memory address) of obserr_f has changed in c__init_obserr_f_pdaf.'


cdef void c__l2g_state_pdaf (int* step, int* domain_p, int* dim_l, double* state_l, int* dim_p, double* state_p):
    state_l_np = np.asarray(<double[:np.prod((dim_l[0]))]> 
                            state_l).reshape((dim_l[0]), order='F')
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')

    state_p_np_tmp = py__l2g_state_pdaf(step[0], domain_p[0], dim_l[0], state_l_np, dim_p[0], state_p_np)

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__l2g_state_pdaf.'


cdef void c__obs_op_f_pdaf (int* step, int* dim_p, int* dim_obs_f, double* state_p, double* m_state_f):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    m_state_f_np = np.asarray(<double[:np.prod((dim_obs_f[0]))]> 
                              m_state_f).reshape((dim_obs_f[0]), order='F')

    m_state_f_np_tmp = py__obs_op_f_pdaf(step[0], dim_p[0], dim_obs_f[0], state_p_np, m_state_f_np)

    m_state_f_np[:] = m_state_f_np_tmp[:]
    cdef double[::1] m_state_f_view = m_state_f_np.ravel(order='F')
    assert m_state_f == &m_state_f_view[0], 'reference (memory address) of m_state_f has changed in c__obs_op_f_pdaf.'


cdef void c__prodrinva_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* a_l, double* c_l):
    obs_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                          obs_l).reshape((dim_obs_l[0]), order='F')
    a_l_np = np.asarray(<double[:np.prod((dim_obs_l[0], rank[0]))]> 
                        a_l).reshape((dim_obs_l[0], rank[0]), order='F')
    c_l_np = np.asarray(<double[:np.prod((dim_obs_l[0], rank[0]))]> 
                        c_l).reshape((dim_obs_l[0], rank[0]), order='F')

    c_l_np_tmp = py__prodrinva_l_pdaf(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np, a_l_np, c_l_np)

    c_l_np[:] = c_l_np_tmp[:]
    cdef double[::1] c_l_view = c_l_np.ravel(order='F')
    assert c_l == &c_l_view[0], 'reference (memory address) of c_l has changed in c__prodrinva_l_pdaf.'


cdef void c__localize_covar_pdaf (int* dim_p, int* dim_obs, double* hp_p, double* hph):
    hp_p_np = np.asarray(<double[:np.prod((dim_obs[0], dim_p[0]))]> 
                         hp_p).reshape((dim_obs[0], dim_p[0]), order='F')
    hph_np = np.asarray(<double[:np.prod((dim_obs[0], dim_obs[0]))]> 
                        hph).reshape((dim_obs[0], dim_obs[0]), order='F')

    hp_p_np_tmp, hph_np_tmp = py__localize_covar_pdaf(dim_p[0], dim_obs[0], hp_p_np, hph_np)

    hp_p_np[:] = hp_p_np_tmp[:]
    cdef double[::1] hp_p_view = hp_p_np.ravel(order='F')
    assert hp_p == &hp_p_view[0], 'reference (memory address) of hp_p has changed in c__localize_covar_pdaf.'
    hph_np[:] = hph_np_tmp[:]
    cdef double[::1] hph_view = hph_np.ravel(order='F')
    assert hph == &hph_view[0], 'reference (memory address) of hph has changed in c__localize_covar_pdaf.'


cdef void c__likelihood_pdaf (int* step, int* dim_obs_p, double* obs_p, double* resid, double* likely):
    obs_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          obs_p).reshape((dim_obs_p[0]), order='F')
    resid_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                          resid).reshape((dim_obs_p[0]), order='F')

    likely[0] = py__likelihood_pdaf(step[0], dim_obs_p[0], obs_p_np, resid_np, likely[0])



cdef void c__likelihood_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* likely_l):
    obs_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                          obs_l).reshape((dim_obs_l[0]), order='F')
    resid_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                            resid_l).reshape((dim_obs_l[0]), order='F')

    likely_l[0] = py__likelihood_l_pdaf(domain_p[0], step[0], dim_obs_l[0], obs_l_np, resid_l_np, likely_l[0])



cdef void c__get_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f):
    observation_f_np = np.asarray(<double[:np.prod((dim_obs_f[0]))]> 
                                  observation_f).reshape((dim_obs_f[0]), order='F')

    observation_f_np_tmp = py__get_obs_f_pdaf(step[0], dim_obs_f[0], observation_f_np)

    observation_f_np[:] = observation_f_np_tmp[:]
    cdef double[::1] observation_f_view = observation_f_np.ravel(order='F')
    assert observation_f == &observation_f_view[0], 'reference (memory address) of observation_f has changed in c__get_obs_f_pdaf.'


cdef void c__cvt_adj_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cv_ens_p, double* ens_p, double* vcv_p, double* cv_p):
    ens_p_np = np.asarray(<double[:np.prod((dim_p[0], dim_ens[0]))]> 
                          ens_p).reshape((dim_p[0], dim_ens[0]), order='F')
    vcv_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                          vcv_p).reshape((dim_p[0]), order='F')
    cv_p_np = np.asarray(<double[:np.prod((dim_cv_ens_p[0]))]> 
                         cv_p).reshape((dim_cv_ens_p[0]), order='F')

    cv_p_np_tmp = py__cvt_adj_ens_pdaf(iter[0], dim_p[0], dim_ens[0], dim_cv_ens_p[0], ens_p_np, vcv_p_np, cv_p_np)

    cv_p_np[:] = cv_p_np_tmp[:]
    cdef double[::1] cv_p_view = cv_p_np.ravel(order='F')
    assert cv_p == &cv_p_view[0], 'reference (memory address) of cv_p has changed in c__cvt_adj_ens_pdaf.'


cdef void c__cvt_adj_pdaf (int* iter, int* dim_p, int* dim_cvec, double* vcv_p, double* cv_p):
    vcv_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                          vcv_p).reshape((dim_p[0]), order='F')
    cv_p_np = np.asarray(<double[:np.prod((dim_cvec[0]))]> 
                         cv_p).reshape((dim_cvec[0]), order='F')

    cv_p_np_tmp = py__cvt_adj_pdaf(iter[0], dim_p[0], dim_cvec[0], vcv_p_np, cv_p_np)

    cv_p_np[:] = cv_p_np_tmp[:]
    cdef double[::1] cv_p_view = cv_p_np.ravel(order='F')
    assert cv_p == &cv_p_view[0], 'reference (memory address) of cv_p has changed in c__cvt_adj_pdaf.'


cdef void c__cvt_pdaf (int* iter, int* dim_p, int* dim_cvec, double* cv_p, double* vv_p):
    cv_p_np = np.asarray(<double[:np.prod((dim_cvec[0]))]> 
                         cv_p).reshape((dim_cvec[0]), order='F')
    vv_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                         vv_p).reshape((dim_p[0]), order='F')

    vv_p_np_tmp = py__cvt_pdaf(iter[0], dim_p[0], dim_cvec[0], cv_p_np, vv_p_np)

    vv_p_np[:] = vv_p_np_tmp[:]
    cdef double[::1] vv_p_view = vv_p_np.ravel(order='F')
    assert vv_p == &vv_p_view[0], 'reference (memory address) of vv_p has changed in c__cvt_pdaf.'


cdef void c__cvt_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cvec_ens, double* ens_p, double* v_p, double* vv_p):
    ens_p_np = np.asarray(<double[:np.prod((dim_p[0], dim_ens[0]))]> 
                          ens_p).reshape((dim_p[0], dim_ens[0]), order='F')
    v_p_np = np.asarray(<double[:np.prod((dim_cvec_ens[0]))]> 
                        v_p).reshape((dim_cvec_ens[0]), order='F')
    vv_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                         vv_p).reshape((dim_p[0]), order='F')

    vv_p_np_tmp = py__cvt_ens_pdaf(iter[0], dim_p[0], dim_ens[0], dim_cvec_ens[0], ens_p_np, v_p_np, vv_p_np)

    vv_p_np[:] = vv_p_np_tmp[:]
    cdef double[::1] vv_p_view = vv_p_np.ravel(order='F')
    assert vv_p == &vv_p_view[0], 'reference (memory address) of vv_p has changed in c__cvt_ens_pdaf.'


cdef void c__obs_op_adj_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    m_state_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                              m_state_p).reshape((dim_obs_p[0]), order='F')

    state_p_np_tmp = py__obs_op_adj_pdaf(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)

    state_p_np[:] = state_p_np_tmp[:]
    cdef double[::1] state_p_view = state_p_np.ravel(order='F')
    assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__obs_op_adj_pdaf.'


cdef void c__obs_op_lin_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p):
    state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                            state_p).reshape((dim_p[0]), order='F')
    m_state_p_np = np.asarray(<double[:np.prod((dim_obs_p[0]))]> 
                              m_state_p).reshape((dim_obs_p[0]), order='F')

    m_state_p_np_tmp = py__obs_op_lin_pdaf(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)

    m_state_p_np[:] = m_state_p_np_tmp[:]
    cdef double[::1] m_state_p_view = m_state_p_np.ravel(order='F')
    assert m_state_p == &m_state_p_view[0], 'reference (memory address) of m_state_p has changed in c__obs_op_lin_pdaf.'


cdef void c__dist_stateinc_pdaf (int* dim_p, double* state_inc_p, int* first, int* steps):
    state_inc_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> 
                                state_inc_p).reshape((dim_p[0]), order='F')

    py__dist_stateinc_pdaf(dim_p[0], state_inc_p_np, first[0], steps[0])



cdef void c__prodrinva_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* gamma, double* likely_l):
    obs_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                          obs_l).reshape((dim_obs_l[0]), order='F')
    resid_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                            resid_l).reshape((dim_obs_l[0]), order='F')

    likely_l[0] = py__prodrinva_hyb_l_pdaf(domain_p[0], step[0], dim_obs_l[0], obs_l_np, resid_l_np, gamma[0], likely_l[0])



cdef void c__likelihood_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* gamma, double* a_l, double* c_l):
    obs_l_np = np.asarray(<double[:np.prod((dim_obs_l[0]))]> 
                          obs_l).reshape((dim_obs_l[0]), order='F')
    a_l_np = np.asarray(<double[:np.prod((dim_obs_l[0], rank[0]))]> 
                        a_l).reshape((dim_obs_l[0], rank[0]), order='F')
    c_l_np = np.asarray(<double[:np.prod((dim_obs_l[0], rank[0]))]> 
                        c_l).reshape((dim_obs_l[0], rank[0]), order='F')

    c_l_np_tmp = py__likelihood_hyb_l_pdaf(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np, gamma[0], a_l_np, c_l_np)

    c_l_np[:] = c_l_np_tmp[:]
    cdef double[::1] c_l_view = c_l_np.ravel(order='F')
    assert c_l == &c_l_view[0], 'reference (memory address) of c_l has changed in c__likelihood_hyb_l_pdaf.'


