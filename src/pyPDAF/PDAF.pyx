import sys
import numpy as np
cimport pyPDAF.UserFunc as c__PDAFcython

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
                sys.stderr.write("Uncaught exception was detected on rank {}. \n".format(
                    mpi4py.MPI.COMM_WORLD.Get_rank()))

                print_exception(exctype, value, traceback)
                sys.stderr.write("\n")
                sys.stderr.flush()
            finally:
                try:
                    mpi4py.MPI.COMM_WORLD.Abort(1)
                except Exception as e:
                    sys.stderr.write("MPI Abort failed, this process will hang.\n")
                    sys.stderr.flush()
                    raise e
        else:
            sys.__excepthook__(exctype, value, traceback)
    except ImportError:
        sys.__excepthook__(exctype, value, traceback)

sys.excepthook = global_except_hook


def assimilate_3dvar (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__prodRinvA_pdaf,
                      py__cvt_pdaf,
                      py__cvt_adj_pdaf,
                      py__obs_op_lin_pdaf,
                      py__obs_op_adj_pdaf,
                      py__prepoststep_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_3dvar or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_3dvar (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__prodRinvA_pdaf,
                              c__PDAFcython.c__cvt_pdaf,
                              c__PDAFcython.c__cvt_adj_pdaf,
                              c__PDAFcython.c__obs_op_lin_pdaf,
                              c__PDAFcython.c__obs_op_adj_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &outflag
                             )

    return outflag

def assimilate_en3dvar_estkf (py__collect_state_pdaf,
                              py__distribute_state_pdaf,
                              py__init_dim_obs_pdaf,
                              py__obs_op_pdaf,
                              py__init_obs_pdaf,
                              py__prodRinvA_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_obsvar_pdaf,
                              py__prepoststep_pdaf,
                              py__next_observation_pdaf
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_en3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__distribute_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodRinvA_pdaf,
                                      c__PDAFcython.c__cvt_ens_pdaf,
                                      c__PDAFcython.c__cvt_adj_ens_pdaf,
                                      c__PDAFcython.c__obs_op_lin_pdaf,
                                      c__PDAFcython.c__obs_op_adj_pdaf,
                                      c__PDAFcython.c__init_obsvar_pdaf,
                                      c__PDAFcython.c__prepoststep_pdaf,
                                      c__PDAFcython.c__next_observation_pdaf,
                                      &outflag
                                     )

    return outflag

def assimilate_en3dvar_lestkf (py__collect_state_pdaf,
                               py__distribute_state_pdaf,
                               py__init_dim_obs_pdaf,
                               py__obs_op_pdaf,
                               py__init_obs_pdaf,
                               py__prodRinvA_pdaf,
                               py__cvt_ens_pdaf,
                               py__cvt_adj_ens_pdaf,
                               py__obs_op_lin_pdaf,
                               py__obs_op_adj_pdaf,
                               py__init_dim_obs_f_pdaf,
                               py__obs_op_f_pdaf,
                               py__init_obs_f_pdaf,
                               py__init_obs_l_pdaf,
                               py__prodRinvA_l_pdaf,
                               py__init_n_domains_p_pdaf,
                               py__init_dim_l_pdaf,
                               py__init_dim_obs_l_pdaf,
                               py__g2l_state_pdaf,
                               py__l2g_state_pdaf,
                               py__g2l_obs_pdaf,
                               py__init_obsvar_pdaf,
                               py__init_obsvar_l_pdaf,
                               py__prepoststep_pdaf,
                               py__next_observation_pdaf,
                               int outflag
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_en3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__init_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdaf_assimilate_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__distribute_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodRinvA_pdaf,
                                       c__PDAFcython.c__cvt_ens_pdaf,
                                       c__PDAFcython.c__cvt_adj_ens_pdaf,
                                       c__PDAFcython.c__obs_op_lin_pdaf,
                                       c__PDAFcython.c__obs_op_adj_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__init_obs_f_pdaf,
                                       c__PDAFcython.c__init_obs_l_pdaf,
                                       c__PDAFcython.c__prodRinvA_l_pdaf,
                                       c__PDAFcython.c__init_n_domains_p_pdaf,
                                       c__PDAFcython.c__init_dim_l_pdaf,
                                       c__PDAFcython.c__init_dim_obs_l_pdaf,
                                       c__PDAFcython.c__g2l_state_pdaf,
                                       c__PDAFcython.c__l2g_state_pdaf,
                                       c__PDAFcython.c__g2l_obs_pdaf,
                                       c__PDAFcython.c__init_obsvar_pdaf,
                                       c__PDAFcython.c__init_obsvar_l_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       c__PDAFcython.c__next_observation_pdaf,
                                       &outflag
                                      )

    return outflag

def assimilate_enkf (py__collect_state_pdaf,
                     py__distribute_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__add_obs_err_pdaf,
                     py__init_obs_covar_pdaf,
                     py__next_observation_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_enkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__add_obs_err_pdaf : Callable[step:int, dim_obs_p:int, C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]]
        Add obs error covariance R to HPH in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Dimension of observation vector
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

    py__init_obs_covar_pdaf : Callable[step:int, dim_obs:int, dim_obs_p:int, covar:float, obs_p : ndarray[tuple[dim_obs_p], np.float64], isdiag:bool]
        Initialize obs. error cov. matrix R in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs:int
            Global size of observation vector
        dim_obs_p:int
            Size of process-local observation vector
        covar:float
            Observation error covariance matrix
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Process-local vector of observations
        isdiag:bool

        Returns
        -------
        covar:float
            Observation error covariance matrix
        isdiag:bool

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_enkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__add_obs_err_pdaf,
                             c__PDAFcython.c__init_obs_covar_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_estkf (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__prepoststep_pdaf,
                      py__prodRinvA_pdaf,
                      py__init_obsvar_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 HV

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_estkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodRinvA_pdaf,
                              c__PDAFcython.c__init_obsvar_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &flag
                             )

    return flag

def assimilate_etkf (py__collect_state_pdaf,
                     py__distribute_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_pdaf,
                     py__init_obsvar_pdaf,
                     py__next_observation_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_etkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 HV

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_etkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_pdaf,
                             c__PDAFcython.c__init_obsvar_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_hyb3dvar_estkf (py__collect_state_pdaf,
                               py__distribute_state_pdaf,
                               py__init_dim_obs_pdaf,
                               py__obs_op_pdaf,
                               py__init_obs_pdaf,
                               py__prodRinvA_pdaf,
                               py__cvt_ens_pdaf,
                               py__cvt_adj_ens_pdaf,
                               py__cvt_pdaf,
                               py__cvt_adj_pdaf,
                               py__obs_op_lin_pdaf,
                               py__obs_op_adj_pdaf,
                               py__init_obsvar_pdaf,
                               py__prepoststep_pdaf,
                               py__next_observation_pdaf
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_hyb3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__distribute_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodRinvA_pdaf,
                                       c__PDAFcython.c__cvt_ens_pdaf,
                                       c__PDAFcython.c__cvt_adj_ens_pdaf,
                                       c__PDAFcython.c__cvt_pdaf,
                                       c__PDAFcython.c__cvt_adj_pdaf,
                                       c__PDAFcython.c__obs_op_lin_pdaf,
                                       c__PDAFcython.c__obs_op_adj_pdaf,
                                       c__PDAFcython.c__init_obsvar_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       c__PDAFcython.c__next_observation_pdaf,
                                       &outflag
                                      )

    return outflag

def assimilate_hyb3dvar_lestkf (py__collect_state_pdaf,
                                py__distribute_state_pdaf,
                                py__init_dim_obs_pdaf,
                                py__obs_op_pdaf,
                                py__init_obs_pdaf,
                                py__prodRinvA_pdaf,
                                py__cvt_ens_pdaf,
                                py__cvt_adj_ens_pdaf,
                                py__cvt_pdaf,
                                py__cvt_adj_pdaf,
                                py__obs_op_lin_pdaf,
                                py__obs_op_adj_pdaf,
                                py__init_dim_obs_f_pdaf,
                                py__obs_op_f_pdaf,
                                py__init_obs_f_pdaf,
                                py__init_obs_l_pdaf,
                                py__prodRinvA_l_pdaf,
                                py__init_n_domains_p_pdaf,
                                py__init_dim_l_pdaf,
                                py__init_dim_obs_l_pdaf,
                                py__g2l_state_pdaf,
                                py__l2g_state_pdaf,
                                py__g2l_obs_pdaf,
                                py__init_obsvar_pdaf,
                                py__init_obsvar_l_pdaf,
                                py__prepoststep_pdaf,
                                py__next_observation_pdaf
                               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_hyb3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__init_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                        c__PDAFcython.c__distribute_state_pdaf,
                                        c__PDAFcython.c__init_dim_obs_pdaf,
                                        c__PDAFcython.c__obs_op_pdaf,
                                        c__PDAFcython.c__init_obs_pdaf,
                                        c__PDAFcython.c__prodRinvA_pdaf,
                                        c__PDAFcython.c__cvt_ens_pdaf,
                                        c__PDAFcython.c__cvt_adj_ens_pdaf,
                                        c__PDAFcython.c__cvt_pdaf,
                                        c__PDAFcython.c__cvt_adj_pdaf,
                                        c__PDAFcython.c__obs_op_lin_pdaf,
                                        c__PDAFcython.c__obs_op_adj_pdaf,
                                        c__PDAFcython.c__init_dim_obs_f_pdaf,
                                        c__PDAFcython.c__obs_op_f_pdaf,
                                        c__PDAFcython.c__init_obs_f_pdaf,
                                        c__PDAFcython.c__init_obs_l_pdaf,
                                        c__PDAFcython.c__prodRinvA_l_pdaf,
                                        c__PDAFcython.c__init_n_domains_p_pdaf,
                                        c__PDAFcython.c__init_dim_l_pdaf,
                                        c__PDAFcython.c__init_dim_obs_l_pdaf,
                                        c__PDAFcython.c__g2l_state_pdaf,
                                        c__PDAFcython.c__l2g_state_pdaf,
                                        c__PDAFcython.c__g2l_obs_pdaf,
                                        c__PDAFcython.c__init_obsvar_pdaf,
                                        c__PDAFcython.c__init_obsvar_l_pdaf,
                                        c__PDAFcython.c__prepoststep_pdaf,
                                        c__PDAFcython.c__next_observation_pdaf,
                                        &outflag
                                       )

    return outflag

def assimilate_lenkf (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__prepoststep_pdaf,
                      py__localize_covar_pdaf,
                      py__add_obs_err_pdaf,
                      py__init_obs_covar_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lenkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__localize_covar_pdaf : Callable[dim_p:int, dim_obs:int, hp_p : ndarray[tuple[dim_obs, dim_p], np.float64], hph : ndarray[tuple[dim_obs, dim_obs], np.float64]]
        Apply localization to HP and HPH^T

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        dim_obs:int
            number of observations
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

        Returns
        -------
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

    py__add_obs_err_pdaf : Callable[step:int, dim_obs_p:int, C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]]
        Add obs error covariance R to HPH in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Dimension of observation vector
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

    py__init_obs_covar_pdaf : Callable[step:int, dim_obs:int, dim_obs_p:int, covar:float, obs_p : ndarray[tuple[dim_obs_p], np.float64], isdiag:bool]
        Initialize obs. error cov. matrix R in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs:int
            Global size of observation vector
        dim_obs_p:int
            Size of process-local observation vector
        covar:float
            Observation error covariance matrix
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Process-local vector of observations
        isdiag:bool

        Returns
        -------
        covar:float
            Observation error covariance matrix
        isdiag:bool

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lenkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__localize_covar_pdaf,
                              c__PDAFcython.c__add_obs_err_pdaf,
                              c__PDAFcython.c__init_obs_covar_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &flag
                             )

    return flag

def assimilate_lestkf (py__collect_state_pdaf,
                       py__distribute_state_pdaf,
                       py__init_dim_obs_pdaf,
                       py__obs_op_pdaf,
                       py__init_obs_pdaf,
                       py__init_obs_l_pdaf,
                       py__prepoststep_pdaf,
                       py__prodRinvA_l_pdaf,
                       py__init_n_domains_p_pdaf,
                       py__init_dim_l_pdaf,
                       py__init_dim_obs_l_pdaf,
                       py__g2l_state_pdaf,
                       py__l2g_state_pdaf,
                       py__g2l_obs_pdaf,
                       py__init_obsvar_pdaf,
                       py__init_obsvar_l_pdaf,
                       py__next_observation_pdaf
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lestkf (c__PDAFcython.c__collect_state_pdaf,
                               c__PDAFcython.c__distribute_state_pdaf,
                               c__PDAFcython.c__init_dim_obs_pdaf,
                               c__PDAFcython.c__obs_op_pdaf,
                               c__PDAFcython.c__init_obs_pdaf,
                               c__PDAFcython.c__init_obs_l_pdaf,
                               c__PDAFcython.c__prepoststep_pdaf,
                               c__PDAFcython.c__prodRinvA_l_pdaf,
                               c__PDAFcython.c__init_n_domains_p_pdaf,
                               c__PDAFcython.c__init_dim_l_pdaf,
                               c__PDAFcython.c__init_dim_obs_l_pdaf,
                               c__PDAFcython.c__g2l_state_pdaf,
                               c__PDAFcython.c__l2g_state_pdaf,
                               c__PDAFcython.c__g2l_obs_pdaf,
                               c__PDAFcython.c__init_obsvar_pdaf,
                               c__PDAFcython.c__init_obsvar_l_pdaf,
                               c__PDAFcython.c__next_observation_pdaf,
                               &flag
                              )

    return flag

def assimilate_letkf (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__init_obs_l_pdaf,
                      py__prepoststep_pdaf,
                      py__prodRinvA_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__init_obsvar_pdaf,
                      py__init_obsvar_l_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_letkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_letkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodRinvA_l_pdaf,
                              c__PDAFcython.c__init_n_domains_p_pdaf,
                              c__PDAFcython.c__init_dim_l_pdaf,
                              c__PDAFcython.c__init_dim_obs_l_pdaf,
                              c__PDAFcython.c__g2l_state_pdaf,
                              c__PDAFcython.c__l2g_state_pdaf,
                              c__PDAFcython.c__g2l_obs_pdaf,
                              c__PDAFcython.c__init_obsvar_pdaf,
                              c__PDAFcython.c__init_obsvar_l_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &flag
                             )

    return flag

def assimilate_lnetf (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_l_pdaf,
                      py__prepoststep_pdaf,
                      py__likelihood_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lnetf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], likely_l:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            nput vector holding the local residual
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lnetf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__likelihood_l_pdaf,
                              c__PDAFcython.c__init_n_domains_p_pdaf,
                              c__PDAFcython.c__init_dim_l_pdaf,
                              c__PDAFcython.c__init_dim_obs_l_pdaf,
                              c__PDAFcython.c__g2l_state_pdaf,
                              c__PDAFcython.c__l2g_state_pdaf,
                              c__PDAFcython.c__g2l_obs_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &flag
                             )

    return flag

def assimilate_lknetf (py__collect_state_pdaf,
                       py__distribute_state_pdaf,
                       py__init_dim_obs_pdaf,
                       py__obs_op_pdaf,
                       py__init_obs_pdaf,
                       py__init_obs_l_pdaf,
                       py__prepoststep_pdaf,
                       py__prodRinvA_l_pdaf,
                       py__prodRinvA_hyb_l_pdaf,
                       py__init_n_domains_p_pdaf,
                       py__init_dim_l_pdaf,
                       py__init_dim_obs_l_pdaf,
                       py__g2l_state_pdaf,
                       py__l2g_state_pdaf,
                       py__g2l_obs_pdaf,
                       py__init_obsvar_pdaf,
                       py__init_obsvar_l_pdaf,
                       py__likelihood_l_pdaf,
                       py__likelihood_hyb_l_pdaf,
                       py__next_observation_pdaf
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lknetf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__prodRinvA_hyb_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], gamma:float, likely_l:float]
        Provide product R^-1 A on local analysis domain with hybrid weight

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            Input vector holding the local residual
        gamma:float
            Hybrid weight provided by PDAF
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__likelihood_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], likely_l:float]
        Compute likelihood

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            nput vector holding the local residual
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__likelihood_hyb_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], gamma:float, A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Compute likelihood with hybrid weight

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        gamma:float
            Hybrid weight provided by PDAF
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.prodRinvA_hyb_l_pdaf = <void*>py__prodRinvA_hyb_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lknetf (c__PDAFcython.c__collect_state_pdaf,
                               c__PDAFcython.c__distribute_state_pdaf,
                               c__PDAFcython.c__init_dim_obs_pdaf,
                               c__PDAFcython.c__obs_op_pdaf,
                               c__PDAFcython.c__init_obs_pdaf,
                               c__PDAFcython.c__init_obs_l_pdaf,
                               c__PDAFcython.c__prepoststep_pdaf,
                               c__PDAFcython.c__prodRinvA_l_pdaf,
                               c__PDAFcython.c__prodRinvA_hyb_l_pdaf,
                               c__PDAFcython.c__init_n_domains_p_pdaf,
                               c__PDAFcython.c__init_dim_l_pdaf,
                               c__PDAFcython.c__init_dim_obs_l_pdaf,
                               c__PDAFcython.c__g2l_state_pdaf,
                               c__PDAFcython.c__l2g_state_pdaf,
                               c__PDAFcython.c__g2l_obs_pdaf,
                               c__PDAFcython.c__init_obsvar_pdaf,
                               c__PDAFcython.c__init_obsvar_l_pdaf,
                               c__PDAFcython.c__likelihood_l_pdaf,
                               c__PDAFcython.c__likelihood_hyb_l_pdaf,
                               c__PDAFcython.c__next_observation_pdaf,
                               &flag
                              )

    return flag

def assimilate_lseik (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__init_obs_l_pdaf,
                      py__prepoststep_pdaf,
                      py__prodRinvA_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__init_obsvar_pdaf,
                      py__init_obsvar_l_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lseik or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lseik (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodRinvA_l_pdaf,
                              c__PDAFcython.c__init_n_domains_p_pdaf,
                              c__PDAFcython.c__init_dim_l_pdaf,
                              c__PDAFcython.c__init_dim_obs_l_pdaf,
                              c__PDAFcython.c__g2l_state_pdaf,
                              c__PDAFcython.c__l2g_state_pdaf,
                              c__PDAFcython.c__g2l_obs_pdaf,
                              c__PDAFcython.c__init_obsvar_pdaf,
                              c__PDAFcython.c__init_obsvar_l_pdaf,
                              c__PDAFcython.c__next_observation_pdaf,
                              &flag
                             )

    return flag

def assimilate_netf (py__collect_state_pdaf,
                     py__distribute_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__likelihood_pdaf,
                     py__next_observation_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_netf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], resid : ndarray[tuple[dim_obs_p], np.float64], likely:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        resid : ndarray[tuple[dim_obs_p], np.float64]
            Input vector holding the residual
        likely:float
            Output value of the likelihood

        Returns
        -------
        likely:float
            Output value of the likelihood

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_netf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__likelihood_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_pf (py__collect_state_pdaf,
                   py__distribute_state_pdaf,
                   py__init_dim_obs_pdaf,
                   py__obs_op_pdaf,
                   py__init_obs_pdaf,
                   py__prepoststep_pdaf,
                   py__likelihood_pdaf,
                   py__next_observation_pdaf
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_pf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], resid : ndarray[tuple[dim_obs_p], np.float64], likely:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        resid : ndarray[tuple[dim_obs_p], np.float64]
            Input vector holding the residual
        likely:float
            Output value of the likelihood

        Returns
        -------
        likely:float
            Output value of the likelihood

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_pf (c__PDAFcython.c__collect_state_pdaf,
                           c__PDAFcython.c__distribute_state_pdaf,
                           c__PDAFcython.c__init_dim_obs_pdaf,
                           c__PDAFcython.c__obs_op_pdaf,
                           c__PDAFcython.c__init_obs_pdaf,
                           c__PDAFcython.c__prepoststep_pdaf,
                           c__PDAFcython.c__likelihood_pdaf,
                           c__PDAFcython.c__next_observation_pdaf,
                           &flag
                          )

    return flag

def assimilate_seek (py__collect_state_pdaf,
                     py__distribute_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_pdaf,
                     py__next_observation_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_seek or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 HV

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_seek (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_seik (py__collect_state_pdaf,
                     py__distribute_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_pdaf,
                     py__next_observation_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_seik or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 HV

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_seik (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_prepost (py__collect_state_pdaf,
                        py__distribute_state_pdaf,
                        py__prepoststep_pdaf,
                        py__next_observation_pdaf
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_prepost or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_prepost (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__distribute_state_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__next_observation_pdaf,
                                &flag
                               )

    return flag

def deallocate ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_deallocate or PDAF source files 

    """
    c__pdaf_deallocate ()

def diag_effsample (double[::1] weights
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_effsample or PDAF source files 

    Parameters
    ----------
    weights : ndarray[tuple[dim_sample], np.float64]
        weights of the samples

    Returns
    -------
    effSample : float
        effecfive sample size
    """
    cdef int dim_sample
    dim_sample = weights.shape[0]


    cdef double effSample

    c__pdaf_diag_effsample (&dim_sample,
                            &weights[0],
                            &effSample
                           )

    return effSample

def diag_ensstats (int element,
                   double[::1] state,
                   double[:,:] ens
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_ensstats or PDAF source files 

    Parameters
    ----------
    element : int
        ID of element to be used. If element=0, mean values over all elements are computed
    state : ndarray[tuple[dim], np.float64]
        State vector
    ens : ndarray[tuple[dim, dim_ens], np.float64]
        State ensemble

    Returns
    -------
    skewness : float
        Skewness of ensemble
    kurtosis : float
        Kurtosis of ensemble
    status : int
        Status flag (0=success)
    """
    cdef double[::1] ens_f = np.asfortranarray(ens).ravel(order="F")
    cdef int dim, dim_ens
    dim = ens.shape[0]
    dim_ens = ens.shape[1]


    cdef double skewness
    cdef double kurtosis
    cdef int status

    c__pdaf_diag_ensstats (&dim,
                           &dim_ens,
                           &element,
                           &state[0],
                           &ens_f[0],
                           &skewness,
                           &kurtosis,
                           &status
                          )

    return skewness, kurtosis, status

def diag_histogram (int ncall,
                    int element,
                    double[::1] state,
                    double[:,:] ens,
                    int[::1] hist
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_histogram or PDAF source files 

    Parameters
    ----------
    ncall : int
        Number of calls to routine
    element : int
        Element of vector used for histogram
    state : ndarray[tuple[dim], np.float64]
        If element=0, all elements are usedState vector
    ens : ndarray[tuple[dim, dim_ens], np.float64]
        State ensemble
    hist : ndarray[tuple[dim_ens+1], np.intc]
        Histogram about the state

    Returns
    -------
    hist : ndarray[tuple[dim_ens+1], np.intc]
         Histogram about the state
    delta : float
        deviation measure from flat histogram
    status : int
        Status flag (0=success)
    """
    cdef double[::1] ens_f = np.asfortranarray(ens).ravel(order="F")
    cdef int dim, dim_ens
    dim = ens.shape[0]
    dim_ens = ens.shape[1]


    cdef double delta
    cdef int status

    c__pdaf_diag_histogram (&ncall,
                            &dim,
                            &dim_ens,
                            &element,
                            &state[0],
                            &ens_f[0],
                            &hist[0],
                            &delta,
                            &status
                           )

    return np.asarray(hist).reshape((dim_ens+1), order='F'), delta, status

def eofcovar (int[::1] dim_fields,
              int[::1] offsets,
              int remove_mstate,
              int do_mv,
              double[:,:] states,
              double[::1] meanstate,
              int verbose
             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_eofcovar or PDAF source files 

    Parameters
    ----------
    dim_fields : ndarray[tuple[nfields], np.intc]
        Size of each field
    offsets : ndarray[tuple[nfields], np.intc]
        Start position of each field
    remove_mstate : int
        1: subtract mean state from statesbefore computing EOFs; 0: don't remove
    do_mv : int
        1: Do multivariate scaling; 0: no scalingnfields, dim_fields and offsets are only used if do_mv=1
    states : ndarray[tuple[dim_state, nstates], np.float64]
        State perturbations
    meanstate : ndarray[tuple[dim_state], np.float64]
        Mean state (only changed if remove_mstate=1)
    verbose : int
        Verbosity flag

    Returns
    -------
    states : ndarray[tuple[dim_state, nstates], np.float64]
         State perturbations
    stddev : ndarray[tuple[nfields], np.float64]
         Standard deviation of field variabilityWithout multivariate scaling (do_mv=0), it is stddev = 1.0
    svals : ndarray[tuple[nstates], np.float64]
         Singular values divided by sqrt(nstates-1)
    svec : ndarray[tuple[dim_state, nstates], np.float64]
         Singular vectors
    meanstate : ndarray[tuple[dim_state], np.float64]
         Mean state (only changed if remove_mstate=1)
    status : int
        Status flag
    """
    cdef double[::1] states_f = np.asfortranarray(states).ravel(order="F")
    cdef int nfields, dim_state, nstates
    dim_state = states.shape[0]
    nstates = states.shape[1]
    nfields = dim_fields.shape[0]


    cdef double [::1] stddev = np.zeros((nfields), dtype=np.float64).ravel()
    cdef double [::1] svals = np.zeros((nstates), dtype=np.float64).ravel()
    cdef double [::1] svec = np.zeros((dim_state, nstates), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_eofcovar (&dim_state,
                      &nstates,
                      &nfields,
                      &dim_fields[0],
                      &offsets[0],
                      &remove_mstate,
                      &do_mv,
                      &states_f[0],
                      &stddev[0],
                      &svals[0],
                      &svec[0],
                      &meanstate[0],
                      &verbose,
                      &status
                     )

    return np.asarray(states).reshape((dim_state, nstates), order='F'), np.asarray(stddev).reshape((nfields), order='F'), np.asarray(svals).reshape((nstates), order='F'), np.asarray(svec).reshape((dim_state, nstates), order='F'), np.asarray(meanstate).reshape((dim_state), order='F'), status

def gather_dim_obs_f (int dim_obs_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_dim_obs_f or PDAF source files 

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension

    Returns
    -------
    dim_obs_f : int
        Full observation dimension
    """

    cdef int dim_obs_f

    c__pdaf_gather_dim_obs_f (&dim_obs_p,
                              &dim_obs_f
                             )

    return dim_obs_f

def gather_obs_f (double[::1] obs_p,
                  int dimobs_f
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f or PDAF source files 

    Parameters
    ----------
    obs_p : ndarray[tuple[dimobs_p], np.float64]
        PE-local vector
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    obs_f : ndarray[tuple[dimobs_f], np.float64]
         Full gathered vector
    status : int
        Status flag:(0) no error(1) when PDAF_gather_dim_obs_f not executed before
    """
    cdef int dimobs_p
    dimobs_p = obs_p.shape[0]


    cdef double [::1] obs_f = np.zeros((dimobs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f (&obs_p[0],
                          &dimobs_p,
                          &obs_f[0],
                          &dimobs_f,
                          &status
                         )

    return np.asarray(obs_f).reshape((dimobs_f), order='F'), status

def gather_obs_f2 (double[:,:] coords_p,
                   int dimobs_f
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f2 or PDAF source files 

    Parameters
    ----------
    coords_p : ndarray[tuple[nrows, dimobs_p], np.float64]
        PE-local array
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    coords_f : ndarray[tuple[nrows, dimobs_f], np.float64]
         Full gathered array
    status : int
        Status flag:(0) no error(1) when PDAF_gather dim_obs_f not executed before
    """
    cdef double[::1] coords_p_f = np.asfortranarray(coords_p).ravel(order="F")
    cdef int nrows, dimobs_p
    nrows = coords_p.shape[0]
    dimobs_p = coords_p.shape[1]


    cdef double [::1] coords_f = np.zeros((nrows, dimobs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f2 (&coords_p_f[0],
                           &dimobs_p,
                           &coords_f[0],
                           &dimobs_f,
                           &nrows,
                           &status
                          )

    return np.asarray(coords_f).reshape((nrows, dimobs_f), order='F'), status

def generate_obs (py__collect_state_pdaf,
                  py__distribute_state_pdaf,
                  py__init_dim_obs_f_pdaf,
                  py__obs_op_f_pdaf,
                  py__get_obs_f_pdaf,
                  py__init_obserr_f_pdaf,
                  py__prepoststep_pdaf,
                  py__next_observation_pdaf
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_generate_obs or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__get_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Provide observation vector to user

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

    py__init_obserr_f_pdaf : Callable[step:int, dim_obs_f:int, obs_f : ndarray[tuple[dim_obs_f], np.float64], obserr_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize vector of observation error standard deviations

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Full dimension of observation vector
        obs_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation vector
        obserr_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation error stddev

        Returns
        -------
        obserr_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation error stddev

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                          c__PDAFcython.c__distribute_state_pdaf,
                          c__PDAFcython.c__init_dim_obs_f_pdaf,
                          c__PDAFcython.c__obs_op_f_pdaf,
                          c__PDAFcython.c__get_obs_f_pdaf,
                          c__PDAFcython.c__init_obserr_f_pdaf,
                          c__PDAFcython.c__prepoststep_pdaf,
                          c__PDAFcython.c__next_observation_pdaf,
                          &flag
                         )

    return flag

def get_assim_flag ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_assim_flag or PDAF source files 


    Returns
    -------
    did_assim : int
        Flag: (1) for assimilation; (0) else
    """

    cdef int did_assim

    c__pdaf_get_assim_flag (&did_assim
                           )

    return did_assim

def get_ensstats ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_ensstats or PDAF source files 


    Returns
    -------
    dims : ndarray[tuple[1], np.intc]
         dimension of pointer
    c_skew_ptr : ndarray[float]
        Pointer to skewness array
    c_kurt_ptr : ndarray[float]
        Pointer to kurtosis array
    status : int
        Status flag
    """

    cdef int [::1] dims = np.zeros((1), dtype=np.intc).ravel()
    cdef double* c_skew_ptr
    cdef double* c_kurt_ptr
    cdef int status

    c__pdaf_get_ensstats (&dims[0],
                          &c_skew_ptr,
                          &c_kurt_ptr,
                          &status
                         )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_skew_ptr).reshape(dims, order='F'), \
           np.asarray(<double[:np.prod(dims)]> c_kurt_ptr).reshape(dims, order='F'), \
           status

def get_localfilter ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_localfilter or PDAF source files 


    Returns
    -------
    lfilter : int
        Whether the filter is domain-localized
    """

    cdef int lfilter

    c__pdaf_get_localfilter (&lfilter
                            )

    return lfilter

def get_memberid (int memberid
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_memberid or PDAF source files 

    Parameters
    ----------
    memberid : int
        Index in the local ensemble

    Returns
    -------
    memberid : int
        Index in the local ensemble
    """

    c__pdaf_get_memberid (&memberid
                         )

    return memberid

def get_obsmemberid (int memberid
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_obsmemberid or PDAF source files 

    Parameters
    ----------
    memberid : int
        Index in the local observed ensemble

    Returns
    -------
    memberid : int
        Index in the local observed ensemble
    """

    c__pdaf_get_obsmemberid (&memberid
                            )

    return memberid

def get_smootherens ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_smootherens or PDAF source files 


    Returns
    -------
    c_sens_point : ndarray[float]
        Pointer to smoother array
    maxlag : int
        Number of past timesteps processed in sens
    dims : ndarray[tuple[3], np.intc]
         dimension of pointer
    status : int
        Status flag
    """

    cdef double* c_sens_point
    cdef int maxlag
    cdef int [::1] dims = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_get_smootherens (&c_sens_point,
                             &maxlag,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           maxlag, status

def get_state (int steps,
               int doexit,
               py__next_observation_pdaf,
               py__distribute_state_pdaf,
               py__prepoststep_pdaf,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_state or PDAF source files 

    Parameters
    ----------
    steps : int
        Flag and number of time steps
    doexit : int
        Whether to exit from forecasts
    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    flag : int
        Status flag

    Returns
    -------
    steps : int
        Flag and number of time steps
    time : float
        current model time
    doexit : int
        Whether to exit from forecasts
    flag : int
        Status flag
    """
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef double time

    c__pdaf_get_state (&steps,
                       &time,
                       &doexit,
                       c__PDAFcython.c__next_observation_pdaf,
                       c__PDAFcython.c__distribute_state_pdaf,
                       c__PDAFcython.c__prepoststep_pdaf,
                       &flag
                      )

    return steps, time, doexit, flag

def init (int filtertype,
          int subtype,
          int stepnull,
          int[::1] param_int,
          double[::1] param_real,
          int COMM_model,
          int COMM_filter,
          int COMM_couple,
          int task_id,
          int n_modeltasks,
          bint in_filterpe,
          py__init_ens_pdaf,
          int in_screen
         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_init or PDAF source files 

    Parameters
    ----------
    filtertype : int
        Type of filter
    subtype : int
        Sub-type of filter
    stepnull : int
        Initial time step of assimilation
    param_int : ndarray[tuple[dim_pint], np.intc]
        Integer parameter array
    param_real : ndarray[tuple[dim_preal], np.float64]
        Real parameter array
    COMM_model : int
        Model communicator
    COMM_filter : int
        Filter communicator
    COMM_couple : int
        Coupling communicator
    task_id : int
        Id of my ensemble task
    n_modeltasks : int
        Number of parallel model tasks
    in_filterpe : bool
        Is my PE a filter-PE?
    py__init_ens_pdaf : Callable[filtertype:int, dim_p:int, dim_ens:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User-supplied routine for ensemble initialization

        Parameters
        ----------
        filtertype:int
            type of filter to initialize
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of ensemble
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local model state
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            array not referenced for ensemble filters
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local model state
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            array not referenced for ensemble filters
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

    in_screen : int
        Control screen output:

    Returns
    -------
    param_int : ndarray[tuple[dim_pint], np.intc]
         Integer parameter array
    param_real : ndarray[tuple[dim_preal], np.float64]
         Real parameter array
    flag : int
        Status flag, 0: no error, error codes:
    """
    cdef int dim_pint, dim_preal
    dim_pint = param_int.shape[0]
    dim_preal = param_real.shape[0]

    c__PDAFcython.init_ens_pdaf = <void*>py__init_ens_pdaf
    c__PDAFcython.init_ens_pdaf_single_member = <void*>py__init_ens_pdaf

    cdef int flag

    if (filtertype == 0) or (filtertype == 200 and subtype == 0):
        c__pdaf_init (&filtertype, &subtype, &stepnull,
                      &param_int[0], &dim_pint,
                      &param_real[0], &dim_preal,
                      &COMM_model, &COMM_filter, &COMM_couple,
                      &task_id, &n_modeltasks, &in_filterpe,
                      c__PDAFcython.c__init_ens_pdaf_single_member,
                      &in_screen, &flag)
    else:
        c__pdaf_init (&filtertype,
                      &subtype,
                      &stepnull,
                      &param_int[0],
                      &dim_pint,
                      &param_real[0],
                      &dim_preal,
                      &COMM_model,
                      &COMM_filter,
                      &COMM_couple,
                      &task_id,
                      &n_modeltasks,
                      &in_filterpe,
                      c__PDAFcython.c__init_ens_pdaf,
                      &in_screen,
                      &flag
                     )

    return np.asarray(param_int).reshape((dim_pint), order='F'), np.asarray(param_real).reshape((dim_preal), order='F'), flag

def local_weight (int wtype,
                  int rtype,
                  double cradius,
                  double sradius,
                  double distance,
                  double[:,:] A,
                  double var_obs,
                  int verbose
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_local_weight or PDAF source files 

    Parameters
    ----------
    wtype : int
        Type of weight function
    rtype : int
        Type of regulated weighting
    cradius : float
        Cut-off radius
    sradius : float
        Support radius
    distance : float
        Distance to observation
    A : ndarray[tuple[nrows, ncols], np.float64]
        Input matrix
    var_obs : float
        Observation variance
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : float
        Weights
    """
    cdef double[::1] A_f = np.asfortranarray(A).ravel(order="F")
    cdef int nrows, ncols
    nrows = A.shape[0]
    ncols = A.shape[1]


    cdef double weight

    c__pdaf_local_weight (&wtype,
                          &rtype,
                          &cradius,
                          &sradius,
                          &distance,
                          &nrows,
                          &ncols,
                          &A_f[0],
                          &var_obs,
                          &weight,
                          &verbose
                         )

    return weight

def print_info (int printtype
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_info or PDAF source files 

    Parameters
    ----------
    printtype : int
        Type of screen output
    """

    c__pdaf_print_info (&printtype
                       )

def put_state_3dvar (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prodRinvA_pdaf,
                     py__cvt_pdaf,
                     py__cvt_adj_pdaf,
                     py__obs_op_lin_pdaf,
                     py__obs_op_adj_pdaf,
                     py__prepoststep_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_3dvar or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_3dvar (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prodRinvA_pdaf,
                             c__PDAFcython.c__cvt_pdaf,
                             c__PDAFcython.c__cvt_adj_pdaf,
                             c__PDAFcython.c__obs_op_lin_pdaf,
                             c__PDAFcython.c__obs_op_adj_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             &outflag
                            )

    return outflag

def put_state_en3dvar_estkf (py__collect_state_pdaf,
                             py__init_dim_obs_pdaf,
                             py__obs_op_pdaf,
                             py__init_obs_pdaf,
                             py__prodRinvA_pdaf,
                             py__cvt_ens_pdaf,
                             py__cvt_adj_ens_pdaf,
                             py__obs_op_lin_pdaf,
                             py__obs_op_adj_pdaf,
                             py__init_obsvar_pdaf,
                             py__prepoststep_pdaf
                            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_en3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                     c__PDAFcython.c__init_dim_obs_pdaf,
                                     c__PDAFcython.c__obs_op_pdaf,
                                     c__PDAFcython.c__init_obs_pdaf,
                                     c__PDAFcython.c__prodRinvA_pdaf,
                                     c__PDAFcython.c__cvt_ens_pdaf,
                                     c__PDAFcython.c__cvt_adj_ens_pdaf,
                                     c__PDAFcython.c__obs_op_lin_pdaf,
                                     c__PDAFcython.c__obs_op_adj_pdaf,
                                     c__PDAFcython.c__init_obsvar_pdaf,
                                     c__PDAFcython.c__prepoststep_pdaf,
                                     &outflag
                                    )

    return outflag

def put_state_en3dvar_lestkf (py__collect_state_pdaf,
                              py__init_dim_obs_pdaf,
                              py__obs_op_pdaf,
                              py__init_obs_pdaf,
                              py__prodRinvA_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_dim_obs_f_pdaf,
                              py__obs_op_f_pdaf,
                              py__init_obs_f_pdaf,
                              py__init_obs_l_pdaf,
                              py__prodRinvA_l_pdaf,
                              py__init_n_domains_p_pdaf,
                              py__init_dim_l_pdaf,
                              py__init_dim_obs_l_pdaf,
                              py__g2l_state_pdaf,
                              py__l2g_state_pdaf,
                              py__g2l_obs_pdaf,
                              py__init_obsvar_pdaf,
                              py__init_obsvar_l_pdaf,
                              py__prepoststep_pdaf
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_en3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__init_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodRinvA_pdaf,
                                      c__PDAFcython.c__cvt_ens_pdaf,
                                      c__PDAFcython.c__cvt_adj_ens_pdaf,
                                      c__PDAFcython.c__obs_op_lin_pdaf,
                                      c__PDAFcython.c__obs_op_adj_pdaf,
                                      c__PDAFcython.c__init_dim_obs_f_pdaf,
                                      c__PDAFcython.c__obs_op_f_pdaf,
                                      c__PDAFcython.c__init_obs_f_pdaf,
                                      c__PDAFcython.c__init_obs_l_pdaf,
                                      c__PDAFcython.c__prodRinvA_l_pdaf,
                                      c__PDAFcython.c__init_n_domains_p_pdaf,
                                      c__PDAFcython.c__init_dim_l_pdaf,
                                      c__PDAFcython.c__init_dim_obs_l_pdaf,
                                      c__PDAFcython.c__g2l_state_pdaf,
                                      c__PDAFcython.c__l2g_state_pdaf,
                                      c__PDAFcython.c__g2l_obs_pdaf,
                                      c__PDAFcython.c__init_obsvar_pdaf,
                                      c__PDAFcython.c__init_obsvar_l_pdaf,
                                      c__PDAFcython.c__prepoststep_pdaf,
                                      &outflag
                                     )

    return outflag

def put_state_enkf (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__add_obs_err_pdaf,
                    py__init_obs_covar_pdaf
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_enkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__add_obs_err_pdaf : Callable[step:int, dim_obs_p:int, C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]]
        Add obs error covariance R to HPH in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Dimension of observation vector
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

    py__init_obs_covar_pdaf : Callable[step:int, dim_obs:int, dim_obs_p:int, covar:float, obs_p : ndarray[tuple[dim_obs_p], np.float64], isdiag:bool]
        Initialize obs. error cov. matrix R in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs:int
            Global size of observation vector
        dim_obs_p:int
            Size of process-local observation vector
        covar:float
            Observation error covariance matrix
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Process-local vector of observations
        isdiag:bool

        Returns
        -------
        covar:float
            Observation error covariance matrix
        isdiag:bool


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf

    cdef int flag

    c__pdaf_put_state_enkf (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__add_obs_err_pdaf,
                            c__PDAFcython.c__init_obs_covar_pdaf,
                            &flag
                           )

    return flag

def put_state_estkf (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_pdaf,
                     py__init_obsvar_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_estkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_pdaf,
                             c__PDAFcython.c__init_obsvar_pdaf,
                             &flag
                            )

    return flag

def put_state_etkf (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__prodRinvA_pdaf,
                    py__init_obsvar_pdaf
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_etkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_etkf (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodRinvA_pdaf,
                            c__PDAFcython.c__init_obsvar_pdaf,
                            &flag
                           )

    return flag

def put_state_generate_obs (py__collect_state_pdaf,
                            py__init_dim_obs_f_pdaf,
                            py__obs_op_f_pdaf,
                            py__get_obs_f_pdaf,
                            py__init_obserr_f_pdaf,
                            py__prepoststep_pdaf
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_generate_obs or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__get_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Provide observation vector to user

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

    py__init_obserr_f_pdaf : Callable[step:int, dim_obs_f:int, obs_f : ndarray[tuple[dim_obs_f], np.float64], obserr_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize vector of observation errors

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Full dimension of observation vector
        obs_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation vector
        obserr_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation error stddev

        Returns
        -------
        obserr_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observation error stddev

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdaf_put_state_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                                    c__PDAFcython.c__init_dim_obs_f_pdaf,
                                    c__PDAFcython.c__obs_op_f_pdaf,
                                    c__PDAFcython.c__get_obs_f_pdaf,
                                    c__PDAFcython.c__init_obserr_f_pdaf,
                                    c__PDAFcython.c__prepoststep_pdaf,
                                    &flag
                                   )

    return flag

def put_state_hyb3dvar_estkf (py__collect_state_pdaf,
                              py__init_dim_obs_pdaf,
                              py__obs_op_pdaf,
                              py__init_obs_pdaf,
                              py__prodRinvA_pdaf,
                              py__cvt_pdaf,
                              py__cvt_adj_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_obsvar_pdaf,
                              py__prepoststep_pdaf
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_hyb3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodRinvA_pdaf,
                                      c__PDAFcython.c__cvt_pdaf,
                                      c__PDAFcython.c__cvt_adj_pdaf,
                                      c__PDAFcython.c__cvt_ens_pdaf,
                                      c__PDAFcython.c__cvt_adj_ens_pdaf,
                                      c__PDAFcython.c__obs_op_lin_pdaf,
                                      c__PDAFcython.c__obs_op_adj_pdaf,
                                      c__PDAFcython.c__init_obsvar_pdaf,
                                      c__PDAFcython.c__prepoststep_pdaf,
                                      &outflag
                                     )

    return outflag

def put_state_hyb3dvar_lestkf (py__collect_state_pdaf,
                               py__init_dim_obs_pdaf,
                               py__obs_op_pdaf,
                               py__init_obs_pdaf,
                               py__prodRinvA_pdaf,
                               py__cvt_ens_pdaf,
                               py__cvt_adj_ens_pdaf,
                               py__cvt_pdaf,
                               py__cvt_adj_pdaf,
                               py__obs_op_lin_pdaf,
                               py__obs_op_adj_pdaf,
                               py__init_dim_obs_f_pdaf,
                               py__obs_op_f_pdaf,
                               py__init_obs_f_pdaf,
                               py__init_obs_l_pdaf,
                               py__prodRinvA_l_pdaf,
                               py__init_n_domains_p_pdaf,
                               py__init_dim_l_pdaf,
                               py__init_dim_obs_l_pdaf,
                               py__g2l_state_pdaf,
                               py__l2g_state_pdaf,
                               py__g2l_obs_pdaf,
                               py__init_obsvar_pdaf,
                               py__init_obsvar_l_pdaf,
                               py__prepoststep_pdaf
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_hyb3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix (ensemble)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix (ensemble var)

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__init_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodRinvA_pdaf,
                                       c__PDAFcython.c__cvt_ens_pdaf,
                                       c__PDAFcython.c__cvt_adj_ens_pdaf,
                                       c__PDAFcython.c__cvt_pdaf,
                                       c__PDAFcython.c__cvt_adj_pdaf,
                                       c__PDAFcython.c__obs_op_lin_pdaf,
                                       c__PDAFcython.c__obs_op_adj_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__init_obs_f_pdaf,
                                       c__PDAFcython.c__init_obs_l_pdaf,
                                       c__PDAFcython.c__prodRinvA_l_pdaf,
                                       c__PDAFcython.c__init_n_domains_p_pdaf,
                                       c__PDAFcython.c__init_dim_l_pdaf,
                                       c__PDAFcython.c__init_dim_obs_l_pdaf,
                                       c__PDAFcython.c__g2l_state_pdaf,
                                       c__PDAFcython.c__l2g_state_pdaf,
                                       c__PDAFcython.c__g2l_obs_pdaf,
                                       c__PDAFcython.c__init_obsvar_pdaf,
                                       c__PDAFcython.c__init_obsvar_l_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       &outflag
                                      )

    return outflag

def put_state_lenkf (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prepoststep_pdaf,
                     py__localize_covar_pdaf,
                     py__add_obs_err_pdaf,
                     py__init_obs_covar_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lenkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__localize_covar_pdaf : Callable[dim_p:int, dim_obs:int, hp_p : ndarray[tuple[dim_obs, dim_p], np.float64], hph : ndarray[tuple[dim_obs, dim_obs], np.float64]]
        Apply localization to HP and HPH^T

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        dim_obs:int
            number of observations
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

        Returns
        -------
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

    py__add_obs_err_pdaf : Callable[step:int, dim_obs_p:int, C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]]
        Add obs error covariance R to HPH in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Dimension of observation vector
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
            Matrix to that observation covariance R is added

    py__init_obs_covar_pdaf : Callable[step:int, dim_obs:int, dim_obs_p:int, covar:float, obs_p : ndarray[tuple[dim_obs_p], np.float64], isdiag:bool]
        Initialize obs. error cov. matrix R in EnKF

        Parameters
        ----------
        step:int
            Current time step
        dim_obs:int
            Global size of observation vector
        dim_obs_p:int
            Size of process-local observation vector
        covar:float
            Observation error covariance matrix
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Process-local vector of observations
        isdiag:bool

        Returns
        -------
        covar:float
            Observation error covariance matrix
        isdiag:bool


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf

    cdef int flag

    c__pdaf_put_state_lenkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__localize_covar_pdaf,
                             c__PDAFcython.c__add_obs_err_pdaf,
                             c__PDAFcython.c__init_obs_covar_pdaf,
                             &flag
                            )

    return flag

def put_state_lestkf (py__collect_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__init_obs_l_pdaf,
                      py__prepoststep_pdaf,
                      py__prodRinvA_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__init_obsvar_pdaf,
                      py__init_obsvar_l_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_lestkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodRinvA_l_pdaf,
                              c__PDAFcython.c__init_n_domains_p_pdaf,
                              c__PDAFcython.c__init_dim_l_pdaf,
                              c__PDAFcython.c__init_dim_obs_l_pdaf,
                              c__PDAFcython.c__g2l_state_pdaf,
                              c__PDAFcython.c__l2g_state_pdaf,
                              c__PDAFcython.c__g2l_obs_pdaf,
                              c__PDAFcython.c__init_obsvar_pdaf,
                              c__PDAFcython.c__init_obsvar_l_pdaf,
                              &flag
                             )

    return flag

def put_state_letkf (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__init_obs_l_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_l_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf,
                     py__g2l_obs_pdaf,
                     py__init_obsvar_pdaf,
                     py__init_obsvar_l_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_letkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_letkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__init_obs_l_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_l_pdaf,
                             c__PDAFcython.c__init_n_domains_p_pdaf,
                             c__PDAFcython.c__init_dim_l_pdaf,
                             c__PDAFcython.c__init_dim_obs_l_pdaf,
                             c__PDAFcython.c__g2l_state_pdaf,
                             c__PDAFcython.c__l2g_state_pdaf,
                             c__PDAFcython.c__g2l_obs_pdaf,
                             c__PDAFcython.c__init_obsvar_pdaf,
                             c__PDAFcython.c__init_obsvar_l_pdaf,
                             &flag
                            )

    return flag

def put_state_lnetf (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_l_pdaf,
                     py__prepoststep_pdaf,
                     py__likelihood_l_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf,
                     py__g2l_obs_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lnetf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], likely_l:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            nput vector holding the local residual
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf

    cdef int outflag

    c__pdaf_put_state_lnetf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_l_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__likelihood_l_pdaf,
                             c__PDAFcython.c__init_n_domains_p_pdaf,
                             c__PDAFcython.c__init_dim_l_pdaf,
                             c__PDAFcython.c__init_dim_obs_l_pdaf,
                             c__PDAFcython.c__g2l_state_pdaf,
                             c__PDAFcython.c__l2g_state_pdaf,
                             c__PDAFcython.c__g2l_obs_pdaf,
                             &outflag
                            )

    return outflag

def put_state_lknetf (py__collect_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__init_obs_l_pdaf,
                      py__prepoststep_pdaf,
                      py__prodRinvA_l_pdaf,
                      py__prodRinvA_hyb_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__init_obsvar_pdaf,
                      py__init_obsvar_l_pdaf,
                      py__likelihood_l_pdaf,
                      py__likelihood_hyb_l_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lknetf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__prodRinvA_hyb_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], gamma:float, likely_l:float]
        Provide product R^-1 A on local analysis domain with hybrid weight

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            Input vector holding the local residual
        gamma:float
            Hybrid weight provided by PDAF
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance

    py__likelihood_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], resid_l : ndarray[tuple[dim_obs_l], np.float64], likely_l:float]
        Compute likelihood

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        resid_l : ndarray[tuple[dim_obs_l], np.float64]
            nput vector holding the local residual
        likely_l:float
            Output value of the local likelihood

        Returns
        -------
        likely_l:float
            Output value of the local likelihood

    py__likelihood_hyb_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], gamma:float, A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Compute likelihood with hybrid weight

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        gamma:float
            Hybrid weight provided by PDAF
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.prodRinvA_hyb_l_pdaf = <void*>py__prodRinvA_hyb_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf

    cdef int outflag

    c__pdaf_put_state_lknetf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodRinvA_l_pdaf,
                              c__PDAFcython.c__prodRinvA_hyb_l_pdaf,
                              c__PDAFcython.c__init_n_domains_p_pdaf,
                              c__PDAFcython.c__init_dim_l_pdaf,
                              c__PDAFcython.c__init_dim_obs_l_pdaf,
                              c__PDAFcython.c__g2l_state_pdaf,
                              c__PDAFcython.c__l2g_state_pdaf,
                              c__PDAFcython.c__g2l_obs_pdaf,
                              c__PDAFcython.c__init_obsvar_pdaf,
                              c__PDAFcython.c__init_obsvar_l_pdaf,
                              c__PDAFcython.c__likelihood_l_pdaf,
                              c__PDAFcython.c__likelihood_hyb_l_pdaf,
                              &outflag
                             )

    return outflag

def put_state_lseik (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__init_obs_l_pdaf,
                     py__prepoststep_pdaf,
                     py__prodRinvA_l_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf,
                     py__g2l_obs_pdaf,
                     py__init_obsvar_pdaf,
                     py__init_obsvar_l_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lseik or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize PE-local observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__init_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, observation_l : ndarray[tuple[dim_obs_l], np.float64]]
        Init. observation vector on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local size of the observation vector
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

        Returns
        -------
        observation_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, rank:int, obs_l : ndarray[tuple[dim_obs_l], np.float64], A_l : ndarray[tuple[dim_obs_l, rank], np.float64], C_l : ndarray[tuple[dim_obs_l, rank], np.float64]]
        Provide product R^-1 A on local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Number of local observations at current time step (i.e. the size of the local observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
        obs_l : ndarray[tuple[dim_obs_l], np.float64]
            Local vector of observations
        A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Input matrix provided by PDAF
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

        Returns
        -------
        C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
            Output matrix

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__g2l_obs_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int, mstate_f : ndarray[tuple[dim_p], np.intc], dim_p:int, mstate_l : ndarray[tuple[dim_l], np.intc], dim_l:int]
        Restrict full obs. vector to local analysis domain

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_f:int
            Size of full observation vector for model sub-domain
        dim_obs_l:int
            Size of observation vector for local analysis domain
        mstate_f : ndarray[tuple[dim_p], np.intc]
            Full observation vector for model sub-domain
        dim_p:int
            Size of full observation vector for model sub-domain
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain
        dim_l:int
            Size of observation vector for local analysis domain

        Returns
        -------
        mstate_l : ndarray[tuple[dim_l], np.intc]
            Observation vector for local analysis domain

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance

    py__init_obsvar_l_pdaf : Callable[domain_p:int, step:int, dim_obs_l:int, obs_l : ndarray[tuple[dim_obs_p], np.float64], dim_obs_p:int, meanvar_l:float]
        Initialize local mean observation error variance

        Parameters
        ----------
        domain_p:int
            Index of current local analysis domain
        step:int
            Current time step
        dim_obs_l:int
            Local dimension of observation vector
        obs_l : ndarray[tuple[dim_obs_p], np.float64]
            Local observation vector
        dim_obs_p:int
            Dimension of local observation vector
        meanvar_l:float
            Mean local observation error variance

        Returns
        -------
        meanvar_l:float
            Mean local observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_l_pdaf = <void*>py__prodRinvA_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_lseik (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__init_obs_l_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodRinvA_l_pdaf,
                             c__PDAFcython.c__init_n_domains_p_pdaf,
                             c__PDAFcython.c__init_dim_l_pdaf,
                             c__PDAFcython.c__init_dim_obs_l_pdaf,
                             c__PDAFcython.c__g2l_state_pdaf,
                             c__PDAFcython.c__l2g_state_pdaf,
                             c__PDAFcython.c__g2l_obs_pdaf,
                             c__PDAFcython.c__init_obsvar_pdaf,
                             c__PDAFcython.c__init_obsvar_l_pdaf,
                             &flag
                            )

    return flag

def put_state_netf (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__likelihood_pdaf
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_netf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], resid : ndarray[tuple[dim_obs_p], np.float64], likely:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        resid : ndarray[tuple[dim_obs_p], np.float64]
            Input vector holding the residual
        likely:float
            Output value of the likelihood

        Returns
        -------
        likely:float
            Output value of the likelihood


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf

    cdef int flag

    c__pdaf_put_state_netf (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__likelihood_pdaf,
                            &flag
                           )

    return flag

def put_state_pf (py__collect_state_pdaf,
                  py__init_dim_obs_pdaf,
                  py__obs_op_pdaf,
                  py__init_obs_pdaf,
                  py__prepoststep_pdaf,
                  py__likelihood_pdaf
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_pf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__likelihood_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], resid : ndarray[tuple[dim_obs_p], np.float64], likely:float]
        Compute observation likelihood for an ensemble member

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        resid : ndarray[tuple[dim_obs_p], np.float64]
            Input vector holding the residual
        likely:float
            Output value of the likelihood

        Returns
        -------
        likely:float
            Output value of the likelihood


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf

    cdef int flag

    c__pdaf_put_state_pf (c__PDAFcython.c__collect_state_pdaf,
                          c__PDAFcython.c__init_dim_obs_pdaf,
                          c__PDAFcython.c__obs_op_pdaf,
                          c__PDAFcython.c__init_obs_pdaf,
                          c__PDAFcython.c__prepoststep_pdaf,
                          c__PDAFcython.c__likelihood_pdaf,
                          &flag
                         )

    return flag

def put_state_prepost (py__collect_state_pdaf,
                       py__prepoststep_pdaf
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_prepost or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdaf_put_state_prepost (c__PDAFcython.c__collect_state_pdaf,
                               c__PDAFcython.c__prepoststep_pdaf,
                               &flag
                              )

    return flag

def put_state_seek (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__prodRinvA_pdaf
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_seek or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 HV

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf

    cdef int flag

    c__pdaf_put_state_seek (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodRinvA_pdaf,
                            &flag
                           )

    return flag

def put_state_seik (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__prodRinvA_pdaf,
                    py__init_obsvar_pdaf
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_seik or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__prodRinvA_pdaf : Callable[step:int, dim_obs_p:int, rank:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], A_p : ndarray[tuple[dim_obs_p, rank], np.float64], C_p : ndarray[tuple[dim_obs_p, rank], np.float64]]
        Provide product R^-1 A

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Number of observations at current time step (i.e. the size of the observation vector)
        rank:int
            Number of the columns in the matrix processes here.This is usually the ensemble size minus one(or the rank of the initial covariance matrix)
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        A_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Input matrix provided by PDAF
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

        Returns
        -------
        C_p : ndarray[tuple[dim_obs_p, rank], np.float64]
            Output matrix

    py__init_obsvar_pdaf : Callable[step:int, dim_obs_p:int, obs_p : ndarray[tuple[dim_obs_p], np.float64], meanvar:float]
        Initialize mean observation error variance

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of observation vector
        obs_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations
        meanvar:float
            Mean observation error variance

        Returns
        -------
        meanvar:float
            Mean observation error variance


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodRinvA_pdaf = <void*>py__prodRinvA_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_seik (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodRinvA_pdaf,
                            c__PDAFcython.c__init_obsvar_pdaf,
                            &flag
                           )

    return flag

def reset_forget (double forget_in
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_reset_forget or PDAF source files 

    Parameters
    ----------
    forget_in : float
        New value of forgetting factor
    """

    c__pdaf_reset_forget (&forget_in
                         )

def SampleEns (double[:,:] modes,
               double[::1] svals,
               double[::1] state,
               int verbose,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_SampleEns or PDAF source files 

    Parameters
    ----------
    modes : ndarray[tuple[dim, dim_ens-1], np.float64]
        Array of EOF modes
    svals : ndarray[tuple[dim_ens-1], np.float64]
        Vector of singular values
    state : ndarray[tuple[dim], np.float64]
        PE-local model state
    verbose : int
        Verbosity flag
    flag : int
        Status flag

    Returns
    -------
    modes : ndarray[tuple[dim, dim_ens-1], np.float64]
         Array of EOF modes
    state : ndarray[tuple[dim], np.float64]
         PE-local model state
    ens : ndarray[tuple[dim, dim_ens], np.float64]
         State ensemble
    flag : int
        Status flag
    """
    cdef double[::1] modes_f = np.asfortranarray(modes).ravel(order="F")
    cdef int dim, dim_ens
    dim = modes.shape[0]
    dim_ens = modes.shape[1]
    dim_ens = dim_ens + 1


    cdef double [::1] ens = np.zeros((dim, dim_ens), dtype=np.float64).ravel()

    c__pdaf_sampleens (&dim,
                       &dim_ens,
                       &modes_f[0],
                       &svals[0],
                       &state[0],
                       &ens[0],
                       &verbose,
                       &flag
                      )

    return np.asarray(modes).reshape((dim, dim_ens-1), order='F'), np.asarray(state).reshape((dim), order='F'), np.asarray(ens).reshape((dim, dim_ens), order='F'), flag

def set_debug_flag (int debugval
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_debug_flag or PDAF source files 

    Parameters
    ----------
    debugval : int
        Value of debugging flag; print debug information for >0
    """

    c__pdaf_set_debug_flag (&debugval
                           )

def set_ens_pointer ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_ens_pointer or PDAF source files 


    Returns
    -------
    c_ens_point : ndarray[float]
        Pointer to smoother array
    dims : ndarray[tuple[2], np.intc]
         dimension of the pointer
    status : int
        Status flag
    """

    cdef double* c_ens_point
    cdef int [::1] dims = np.zeros((2), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_ens_pointer (&c_ens_point,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_ens_point).reshape(dims, order='F'), \
           status

def set_smootherens (int maxlag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_smootherens or PDAF source files 

    Parameters
    ----------
    maxlag : int
        Number of past timesteps processed in sens

    Returns
    -------
    c_sens_point : ndarray[float]
        Pointer to smoother array
    dims : ndarray[tuple[3], np.intc]
         dimension of the pointer
    status : int
        Status flag
    """

    cdef double* c_sens_point
    cdef int [::1] dims = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_smootherens (&c_sens_point,
                             &maxlag,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           status

def seik_TtimesA (double[:,:] A
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_seik_TtimesA or PDAF source files 

    Parameters
    ----------
    A : ndarray[tuple[rank, dim_col], np.float64]
        Input matrix

    Returns
    -------
    B : ndarray[tuple[rank+1, dim_col], np.float64]
         Output matrix (TA)
    """
    cdef double[::1] A_f = np.asfortranarray(A).ravel(order="F")
    cdef int rank, dim_col
    rank = A.shape[0]
    dim_col = A.shape[1]


    cdef double [::1] B = np.zeros((rank+1, dim_col), dtype=np.float64).ravel()

    c__pdaf_seik_ttimesa (&rank,
                          &dim_col,
                          &A_f[0],
                          &B[0]
                         )

    return np.asarray(B).reshape((rank+1, dim_col), order='F')

def etkf_Tleft (double[:,:] A
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_etkf_Tleft or PDAF source files 

    Parameters
    ----------
    A : ndarray[tuple[dim_ens, dim], np.float64]
        Input/output matrix

    Returns
    -------
    A : ndarray[tuple[dim_ens, dim], np.float64]
         Input/output matrix
    """
    cdef double[::1] A_f = np.asfortranarray(A).ravel(order="F")
    cdef int dim_ens, dim
    dim_ens = A.shape[0]
    dim = A.shape[1]


    c__pdaf_etkf_tleft (&dim_ens,
                        &dim,
                        &A_f[0]
                       )

    return np.asarray(A).reshape((dim_ens, dim), order='F')

def estkf_OmegaA (double[:,:] A
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_estkf_OmegaA or PDAF source files 

    Parameters
    ----------
    A : ndarray[tuple[rank, dim_col], np.float64]
        Input matrix

    Returns
    -------
    B : ndarray[tuple[rank+1, dim_col], np.float64]
         Output matrix (TA)
    """
    cdef double[::1] A_f = np.asfortranarray(A).ravel(order="F")
    cdef int rank, dim_col
    rank = A.shape[0]
    dim_col = A.shape[1]


    cdef double [::1] B = np.zeros((rank+1, dim_col), dtype=np.float64).ravel()

    c__pdaf_estkf_omegaa (&rank,
                          &dim_col,
                          &A_f[0],
                          &B[0]
                         )

    return np.asarray(B).reshape((rank+1, dim_col), order='F')

def enkf_omega (int[::1] seed,
                double[:,:] omega,
                double norm,
                int otype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_enkf_omega or PDAF source files 

    Parameters
    ----------
    seed : ndarray[tuple[4], np.intc]
        Seed for random number generation
    omega : ndarray[tuple[dim_ens, r], np.float64]
        Random matrix
    norm : float
        Norm for ensemble transformation
    otype : int
        Type of omega
    screen : int
        Verbosity flag

    Returns
    -------
    omega : ndarray[tuple[dim_ens, r], np.float64]
         Random matrix
    norm : float
        Norm for ensemble transformation
    """
    cdef double[::1] omega_f = np.asfortranarray(omega).ravel(order="F")
    cdef int dim_ens, r
    dim_ens = omega.shape[0]
    r = omega.shape[1]


    c__pdaf_enkf_omega (&seed[0],
                        &r,
                        &dim_ens,
                        &omega_f[0],
                        &norm,
                        &otype,
                        &screen
                       )

    return np.asarray(omega).reshape((dim_ens, r), order='F'), norm

def seik_omega (double[:,:] omega,
                int omegatype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_seik_omega or PDAF source files 

    Parameters
    ----------
    omega : ndarray[tuple[rank+1, rank], np.float64]
        Matrix Omega
    omegatype : int
        Select type of omega
    screen : int
        Verbosity flag

    Returns
    -------
    omega : ndarray[tuple[rank+1, rank], np.float64]
         Matrix Omega
    """
    cdef double[::1] omega_f = np.asfortranarray(omega).ravel(order="F")
    cdef int rank
    rank = omega.shape[0]
    _ = omega.shape[1]
    rank = rank - 1


    c__pdaf_seik_omega (&rank,
                        &omega_f[0],
                        &omegatype,
                        &screen
                       )

    return np.asarray(omega).reshape((rank+1, rank), order='F')

def incremental (int steps,
                 py__dist_stateinc_pdaf
                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_incremental or PDAF source files 

    Parameters
    ----------
    steps : int
        Time steps over which increment is distributed
    py__dist_stateinc_pdaf : Callable[dim_p:int, state_inc_p : ndarray[tuple[dim_p], np.float64], first:int, steps:int]
        Add state increment during integration

        Parameters
        ----------
        dim_p:int
            Dimension of PE-local state
        state_inc_p : ndarray[tuple[dim_p], np.float64]
            PE-local state vector
        first:int
            Flag for first call of each forecast
        steps:int
            number of time steps in forecast

        Returns
        -------

    """
    c__PDAFcython.dist_stateinc_pdaf = <void*>py__dist_stateinc_pdaf

    c__pdaf_incremental (&steps,
                         c__PDAFcython.c__dist_stateinc_pdaf
                        )

def add_increment (double[::1] state_p
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_add_increment or PDAF source files 

    Parameters
    ----------
    state_p : ndarray[tuple[dim_p], np.float64]
        State vector

    Returns
    -------
    state_p : ndarray[tuple[dim_p], np.float64]
         State vector
    """
    cdef int dim_p
    dim_p = state_p.shape[0]


    c__pdaf_add_increment (&dim_p,
                           &state_p[0]
                          )

    return np.asarray(state_p).reshape((dim_p), order='F')

def local_weights (int wtype,
                   double cradius,
                   double sradius,
                   double[::1] distance,
                   int verbose
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_local_weights or PDAF source files 

    Parameters
    ----------
    wtype : int
        Type of weight function(0): unit weight (=1 up to distance=cradius)(1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)(2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
    cradius : float
        Parameter for cut-off
    sradius : float
        Support radius
    distance : ndarray[tuple[dim], np.float64]
        Array holding distances
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : ndarray[tuple[dim], np.float64]
         Array for weights
    """
    cdef int dim
    dim = distance.shape[0]


    cdef double [::1] weight = np.zeros((dim), dtype=np.float64).ravel()

    c__pdaf_local_weights (&wtype,
                           &cradius,
                           &sradius,
                           &dim,
                           &distance[0],
                           &weight[0],
                           &verbose
                          )

    return np.asarray(weight).reshape((dim), order='F')

def diag_CRPS (int element,
               double[:,:] oens,
               double[::1] obs
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_CRPS or PDAF source files 

    Parameters
    ----------
    element : int
        ID of element to be used. If element=0, mean values over all elements are computed
    oens : ndarray[tuple[dim, dim_ens], np.float64]
        State ensemble
    obs : ndarray[tuple[dim], np.float64]
        State ensemble

    Returns
    -------
    CRPS : float
        CRPS
    reli : float
        Reliability
    resol : float
        resolution
    uncert : float
        uncertainty
    status : int
        Status flag (0=success)
    """
    cdef double[::1] oens_f = np.asfortranarray(oens).ravel(order="F")
    cdef int dim, dim_ens
    dim = oens.shape[0]
    dim_ens = oens.shape[1]


    cdef double CRPS
    cdef double reli
    cdef double resol
    cdef double uncert
    cdef int status

    c__pdaf_diag_crps (&dim,
                       &dim_ens,
                       &element,
                       &oens_f[0],
                       &obs[0],
                       &CRPS,
                       &reli,
                       &resol,
                       &uncert,
                       &status
                      )

    return CRPS, reli, resol, uncert, status

def force_analysis ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_force_analysis or PDAF source files 

    """
    c__pdaf_force_analysis ()

def gather_obs_f2_flex (int dim_obs_f,
                        double[:,:] coords_p
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f2_flex or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        Full observation dimension
    coords_p : ndarray[tuple[nrows, dim_obs_p], np.float64]
        PE-local array

    Returns
    -------
    coords_f : ndarray[tuple[nrows, dim_obs_f], np.float64]
         Full gathered array
    status : int
        Status flag: (0) no error
    """
    cdef double[::1] coords_p_f = np.asfortranarray(coords_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = coords_p.shape[0]
    dim_obs_p = coords_p.shape[1]


    cdef double [::1] coords_f = np.zeros((nrows, dim_obs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f2_flex (&dim_obs_p,
                                &dim_obs_f,
                                &coords_p_f[0],
                                &coords_f[0],
                                &nrows,
                                &status
                               )

    return np.asarray(coords_f).reshape((nrows, dim_obs_f), order='F'), status

def gather_obs_f_flex (int dim_obs_f,
                       double[::1] obs_p
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f_flex or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        Full observation dimension
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
        PE-local vector

    Returns
    -------
    obs_f : ndarray[tuple[dim_obs_f], np.float64]
         Full gathered vector
    status : int
        Status flag: (0) no error
    """
    cdef int dim_obs_p
    dim_obs_p = obs_p.shape[0]


    cdef double [::1] obs_f = np.zeros((dim_obs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f_flex (&dim_obs_p,
                               &dim_obs_f,
                               &obs_p[0],
                               &obs_f[0],
                               &status
                              )

    return np.asarray(obs_f).reshape((dim_obs_f), order='F'), status

def prepost (py__collect_state_pdaf,
             py__distribute_state_pdaf,
             py__prepoststep_pdaf,
             py__next_observation_pdaf
            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_prepost or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Routine to provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_prepost (c__PDAFcython.c__collect_state_pdaf,
                     c__PDAFcython.c__distribute_state_pdaf,
                     c__PDAFcython.c__prepoststep_pdaf,
                     c__PDAFcython.c__next_observation_pdaf,
                     &outflag
                    )

    return outflag

def set_memberid (int memberid
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_memberid or PDAF source files 

    Parameters
    ----------
    memberid : int
        Index in the local ensemble

    Returns
    -------
    memberid : int
        Index in the local ensemble
    """

    c__pdaf_set_memberid (&memberid
                         )

    return memberid

def set_comm_pdaf (int in_COMM_pdaf
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_comm_pdaf or PDAF source files 

    Parameters
    ----------
    in_COMM_pdaf : int
        MPI communicator for PDAF
    """

    c__pdaf_set_comm_pdaf (&in_COMM_pdaf
                          )

def set_offline_mode (int screen
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_offline_mode or PDAF source files 

    Parameters
    ----------
    screen : int
        Verbosity flag
    """

    c__pdaf_set_offline_mode (&screen
                             )

def print_domain_stats (int n_domains_p
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_domain_stats or PDAF source files 

    Parameters
    ----------
    n_domains_p : int
        Number of PE-local analysis domains
    """

    c__pdaf_print_domain_stats (&n_domains_p
                               )

def init_local_obsstats ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_init_local_obsstats or PDAF source files 

    """
    c__pdaf_init_local_obsstats ()

def incr_local_obsstats (int dim_obs_l
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_incr_local_obsstats or PDAF source files 

    Parameters
    ----------
    dim_obs_l : int
        Number of locally assimilated observations
    """

    c__pdaf_incr_local_obsstats (&dim_obs_l
                                )

def print_local_obsstats (int screen
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_local_obsstats or PDAF source files 

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    n_domains_with_obs : int
        number of domains with observations
    """

    cdef int n_domains_with_obs

    c__pdaf_print_local_obsstats (&screen,
                                  &n_domains_with_obs
                                 )

    return n_domains_with_obs

def omit_obs_omi (double[::1] state_p,
                  double[:,:] ens_p,
                  double[::1] obs_p,
                  py__init_obs_pdaf,
                  py__obs_op_pdaf,
                  int compute_mean,
                  int screen
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_omit_obs_omi or PDAF source files 

    Parameters
    ----------
    state_p : ndarray[tuple[dim_p], np.float64]
        on exit: PE-local forecast mean state
    ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
        PE-local state ensemble
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
        PE-local observation vector
    py__init_obs_pdaf : Callable[step:int, dim_obs_p:int, observation_p : ndarray[tuple[dim_obs_p], np.float64]]
        Initialize observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_p:int
            Size of the observation vector
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

        Returns
        -------
        observation_p : ndarray[tuple[dim_obs_p], np.float64]
            Vector of observations

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    compute_mean : int
        (1) compute mean; (0) state_p holds mean
    screen : int
        Verbosity flag

    Returns
    -------
    state_p : ndarray[tuple[dim_p], np.float64]
         on exit: PE-local forecast mean state
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
         PE-local observation vector
    """
    cdef double[::1] ens_p_f = np.asfortranarray(ens_p).ravel(order="F")
    cdef int dim_p, dim_ens, dim_obs_p
    dim_p = ens_p.shape[0]
    dim_ens = ens_p.shape[1]
    dim_obs_p = obs_p.shape[0]

    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf

    c__pdaf_omit_obs_omi (&dim_p,
                          &dim_obs_p,
                          &dim_ens,
                          &state_p[0],
                          &ens_p_f[0],
                          &obs_p[0],
                          c__PDAFcython.c__init_obs_pdaf,
                          c__PDAFcython.c__obs_op_pdaf,
                          &compute_mean,
                          &screen
                         )

    return np.asarray(state_p).reshape((dim_p), order='F'), np.asarray(obs_p).reshape((dim_obs_p), order='F')

def omi_init (int n_obs
             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init or PDAF source files 

    Parameters
    ----------
    n_obs : int
        number of observations
    """

    c__pdafomi_init (&n_obs
                    )

def omi_set_doassim (int i_obs,
                     int doassim
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_doassim or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    doassim : int
        setter value
    """

    c__pdafomi_set_doassim (&i_obs,
                            &doassim
                           )

def omi_set_disttype (int i_obs,
                      int disttype
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_disttype or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    disttype : int
        setter value
    """

    c__pdafomi_set_disttype (&i_obs,
                             &disttype
                            )

def omi_set_ncoord (int i_obs,
                    int ncoord
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_ncoord or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    ncoord : int
        setter value
    """

    c__pdafomi_set_ncoord (&i_obs,
                           &ncoord
                          )

def omi_set_id_obs_p (int i_obs,
                      int[:,:] id_obs_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_id_obs_p or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    id_obs_p : ndarray[tuple[nrows, dim_obs_p], np.intc]
        setter value
    """
    cdef int[::1] id_obs_p_f = np.asfortranarray(id_obs_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = id_obs_p.shape[0]
    dim_obs_p = id_obs_p.shape[1]


    c__pdafomi_set_id_obs_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &id_obs_p_f[0]
                            )

def omi_set_icoeff_p (int i_obs,
                      double[:,:] icoeff_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_icoeff_p or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    icoeff_p : ndarray[tuple[nrows, dim_obs_p], np.float64]
        setter value
    """
    cdef double[::1] icoeff_p_f = np.asfortranarray(icoeff_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = icoeff_p.shape[0]
    dim_obs_p = icoeff_p.shape[1]


    c__pdafomi_set_icoeff_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &icoeff_p_f[0]
                            )

def omi_set_domainsize (int i_obs,
                        double[::1] domainsize
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_domainsize or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    domainsize : ndarray[tuple[ncoord], np.float64]
        setter value
    """
    cdef int ncoord
    ncoord = domainsize.shape[0]


    c__pdafomi_set_domainsize (&i_obs,
                               &ncoord,
                               &domainsize[0]
                              )

def omi_set_obs_err_type (int i_obs,
                          int obs_err_type
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_obs_err_type or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_err_type : int
        setter value
    """

    c__pdafomi_set_obs_err_type (&i_obs,
                                 &obs_err_type
                                )

def omi_set_use_global_obs (int i_obs,
                            int use_global_obs
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_use_global_obs or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    use_global_obs : int
        setter value
    """

    c__pdafomi_set_use_global_obs (&i_obs,
                                   &use_global_obs
                                  )

def omi_set_inno_omit (int i_obs,
                       double inno_omit
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_inno_omit or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit : float
        setter value
    """

    c__pdafomi_set_inno_omit (&i_obs,
                              &inno_omit
                             )

def omi_set_inno_omit_ivar (int i_obs,
                            double inno_omit_ivar
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_inno_omit_ivar or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit_ivar : float
        setter value
    """

    c__pdafomi_set_inno_omit_ivar (&i_obs,
                                   &inno_omit_ivar
                                  )

def omi_gather_obs (int i_obs,
                    double[::1] obs_p,
                    double[::1] ivar_obs_p,
                    double[:,:] ocoord_p,
                    double cradius
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/pdafomi_gather_obs or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
        pe-local observation vector
    ivar_obs_p : ndarray[tuple[dim_obs_p], np.float64]
        pe-local inverse observation error variance
    ocoord_p : ndarray[tuple[thisobs(i_obs)%ncoord, dim_obs_p], np.float64]
        pe-local observation coordinates
    cradius : float
        localization radius

    Returns
    -------
    dim_obs : int
        Full number of observations
    """
    cdef double[::1] ocoord_p_f = np.asfortranarray(ocoord_p).ravel(order="F")
    cdef int dim_obs_p
    _ = ocoord_p.shape[0]
    dim_obs_p = ocoord_p.shape[1]


    cdef int dim_obs

    c__pdafomi_gather_obs (&i_obs,
                           &dim_obs_p,
                           &obs_p[0],
                           &ivar_obs_p[0],
                           &ocoord_p_f[0],
                           &cradius,
                           &dim_obs
                          )

    return dim_obs

def omi_gather_obsstate (int i_obs,
                         double[::1] obsstate_p,
                         double[::1] obsstate_f
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_gather_obsstate or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    obsstate_p : ndarray[tuple[thisobs(i_obs)%dim_obs_p], np.float64]
        Vector of process-local observed state
    obsstate_f : ndarray[tuple[nobs_f_all], np.float64]
        Full observed vector for all types

    Returns
    -------
    obsstate_f : ndarray[tuple[nobs_f_all], np.float64]
         Full observed vector for all types
    """
    cdef int nobs_f_all
    nobs_f_all = obsstate_f.shape[0]


    c__pdafomi_gather_obsstate (&i_obs,
                                &obsstate_p[0],
                                &obsstate_f[0],
                                &nobs_f_all
                               )

    return np.asarray(obsstate_f).reshape((nobs_f_all), order='F')

def omi_set_domain_limits (double[:,:] lim_coords
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_domain_limits or PDAF source files 

    Parameters
    ----------
    lim_coords : ndarray[tuple[2, 2], np.float64]
        geographic coordinate array (1: longitude, 2: latitude)
    """
    cdef double[::1] lim_coords_f = np.asfortranarray(lim_coords).ravel(order="F")

    c__pdafomi_set_domain_limits (&lim_coords_f[0]
                                 )

def omi_set_debug_flag (int debugval
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_debug_flag or PDAF source files 

    Parameters
    ----------
    debugval : int
        Value for debugging flag
    """

    c__pdafomi_set_debug_flag (&debugval
                              )

def omi_deallocate_obs (int i_obs
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_deallocate_obs or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    """

    c__pdafomi_deallocate_obs (&i_obs
                              )

def omi_obs_op_gridpoint (int i_obs,
                          double[::1] state_p,
                          double[::1] obs_f_all
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_gridpoint or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
         Full observed state for all observation types (nobs_f_all)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_gridpoint (&i_obs,
                                 &state_p[0],
                                 &dim_p,
                                 &obs_f_all[0],
                                 &nobs_f_all
                                )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_gridavg (int i_obs,
                        int nrows,
                        double[::1] state_p,
                        double[::1] obs_f_all
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_gridavg or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
         Full observed state for all observation types (nobs_f_all)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_gridavg (&i_obs,
                               &nrows,
                               &state_p[0],
                               &dim_p,
                               &obs_f_all[0],
                               &nobs_f_all
                              )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_interp_lin (int i_obs,
                           int nrows,
                           double[::1] state_p,
                           double[::1] obs_f_all
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_interp_lin or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
         Full observed state for all observation types (nobs_f_all)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_interp_lin (&i_obs,
                                  &nrows,
                                  &state_p[0],
                                  &dim_p,
                                  &obs_f_all[0],
                                  &nobs_f_all
                                 )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_adj_gridavg (int i_obs,
                            int nrows,
                            double[::1] state_p,
                            double[::1] obs_f_all
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_gridavg or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[tuple[dim_p], np.float64]
         PE-local model state (dim_p)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_gridavg (&i_obs,
                                   &nrows,
                                   &state_p[0],
                                   &dim_p,
                                   &obs_f_all[0],
                                   &nobs_f_all
                                  )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_obs_op_adj_gridpoint (int i_obs,
                              double[::1] state_p,
                              double[::1] obs_f_all
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_gridpoint or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[tuple[dim_p], np.float64]
         PE-local model state (dim_p)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_gridpoint (&i_obs,
                                     &state_p[0],
                                     &dim_p,
                                     &obs_f_all[0],
                                     &nobs_f_all
                                    )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_obs_op_adj_interp_lin (int i_obs,
                               int nrows,
                               double[::1] state_p,
                               double[::1] obs_f_all
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_interp_lin or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p], np.float64]
        PE-local model state (dim_p)
    obs_f_all : ndarray[tuple[nobs_f_all], np.float64]
        Full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[tuple[dim_p], np.float64]
         PE-local model state (dim_p)
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_interp_lin (&i_obs,
                                      &nrows,
                                      &state_p[0],
                                      &dim_p,
                                      &obs_f_all[0],
                                      &nobs_f_all
                                     )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_get_interp_coeff_tri (double[:,:] gpc,
                              double[::1] oc,
                              double[::1] icoeff
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_tri or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[tuple[3, 2], np.float64]
        Coordinates of grid points; dim(3,2)
    oc : ndarray[tuple[2], np.float64]
        3 rows; each containing lon and lat coordinatesCoordinates of observation; dim(2)
    icoeff : ndarray[tuple[3], np.float64]
        Interpolation coefficients; dim(3)

    Returns
    -------
    icoeff : ndarray[tuple[3], np.float64]
         Interpolation coefficients; dim(3)
    """
    cdef double[::1] gpc_f = np.asfortranarray(gpc).ravel(order="F")

    c__pdafomi_get_interp_coeff_tri (&gpc_f[0],
                                     &oc[0],
                                     &icoeff[0]
                                    )

    return np.asarray(icoeff).reshape((3), order='F')

def omi_get_interp_coeff_lin1D (double[::1] gpc,
                                double oc,
                                double[::1] icoeff
                               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_lin1D or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[tuple[2], np.float64]
        Coordinates of grid points (dim=2)
    oc : float
        Coordinates of observation
    icoeff : ndarray[tuple[2], np.float64]
        Interpolation coefficients (dim=2)

    Returns
    -------
    icoeff : ndarray[tuple[2], np.float64]
         Interpolation coefficients (dim=2)
    """

    c__pdafomi_get_interp_coeff_lin1d (&gpc[0],
                                       &oc,
                                       &icoeff[0]
                                      )

    return np.asarray(icoeff).reshape((2), order='F')

def omi_get_interp_coeff_lin (double[:,:] gpc,
                              double[::1] oc,
                              double[::1] icoeff
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_lin or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[tuple[num_gp, n_dim], np.float64]
        Coordinates of grid points
    oc : ndarray[tuple[n_dim], np.float64]
        Coordinates of observation
    icoeff : ndarray[tuple[num_gp], np.float64]
        Interpolation coefficients (num_gp)

    Returns
    -------
    icoeff : ndarray[tuple[num_gp], np.float64]
         Interpolation coefficients (num_gp)
    """
    cdef double[::1] gpc_f = np.asfortranarray(gpc).ravel(order="F")
    cdef int num_gp, n_dim
    num_gp = gpc.shape[0]
    n_dim = gpc.shape[1]


    c__pdafomi_get_interp_coeff_lin (&num_gp,
                                     &n_dim,
                                     &gpc_f[0],
                                     &oc[0],
                                     &icoeff[0]
                                    )

    return np.asarray(icoeff).reshape((num_gp), order='F')

def omi_assimilate_3dvar (py__collect_state_pdaf,
                          py__distribute_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__cvt_pdaf,
                          py__cvt_adj_pdaf,
                          py__obs_op_lin_pdaf,
                          py__obs_op_adj_pdaf,
                          py__prepoststep_pdaf,
                          py__next_observation_pdaf,
                          int outflag
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_3dvar or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdafomi_assimilate_3dvar (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__distribute_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__cvt_pdaf,
                                 c__PDAFcython.c__cvt_adj_pdaf,
                                 c__PDAFcython.c__obs_op_lin_pdaf,
                                 c__PDAFcython.c__obs_op_adj_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 c__PDAFcython.c__next_observation_pdaf,
                                 &outflag
                                )

    return outflag

def omi_assimilate_en3dvar_estkf (py__collect_state_pdaf,
                                  py__distribute_state_pdaf,
                                  py__init_dim_obs_pdaf,
                                  py__obs_op_pdaf,
                                  py__cvt_ens_pdaf,
                                  py__cvt_adj_ens_pdaf,
                                  py__obs_op_lin_pdaf,
                                  py__obs_op_adj_pdaf,
                                  py__prepoststep_pdaf,
                                  py__next_observation_pdaf,
                                  int outflag
                                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_en3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdafomi_assimilate_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                         c__PDAFcython.c__distribute_state_pdaf,
                                         c__PDAFcython.c__init_dim_obs_pdaf,
                                         c__PDAFcython.c__obs_op_pdaf,
                                         c__PDAFcython.c__cvt_ens_pdaf,
                                         c__PDAFcython.c__cvt_adj_ens_pdaf,
                                         c__PDAFcython.c__obs_op_lin_pdaf,
                                         c__PDAFcython.c__obs_op_adj_pdaf,
                                         c__PDAFcython.c__prepoststep_pdaf,
                                         c__PDAFcython.c__next_observation_pdaf,
                                         &outflag
                                        )

    return outflag

def omi_assimilate_en3dvar_lestkf (py__collect_state_pdaf,
                                   py__distribute_state_pdaf,
                                   py__init_dim_obs_f_pdaf,
                                   py__obs_op_f_pdaf,
                                   py__cvt_ens_pdaf,
                                   py__cvt_adj_ens_pdaf,
                                   py__obs_op_lin_pdaf,
                                   py__obs_op_adj_pdaf,
                                   py__init_n_domains_p_pdaf,
                                   py__init_dim_l_pdaf,
                                   py__init_dim_obs_l_pdaf,
                                   py__g2l_state_pdaf,
                                   py__l2g_state_pdaf,
                                   py__prepoststep_pdaf,
                                   py__next_observation_pdaf,
                                   int outflag
                                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_en3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of full observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Full observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize local dimimension of obs. vector

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from local state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdafomi_assimilate_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                          c__PDAFcython.c__distribute_state_pdaf,
                                          c__PDAFcython.c__init_dim_obs_f_pdaf,
                                          c__PDAFcython.c__obs_op_f_pdaf,
                                          c__PDAFcython.c__cvt_ens_pdaf,
                                          c__PDAFcython.c__cvt_adj_ens_pdaf,
                                          c__PDAFcython.c__obs_op_lin_pdaf,
                                          c__PDAFcython.c__obs_op_adj_pdaf,
                                          c__PDAFcython.c__init_n_domains_p_pdaf,
                                          c__PDAFcython.c__init_dim_l_pdaf,
                                          c__PDAFcython.c__init_dim_obs_l_pdaf,
                                          c__PDAFcython.c__g2l_state_pdaf,
                                          c__PDAFcython.c__l2g_state_pdaf,
                                          c__PDAFcython.c__prepoststep_pdaf,
                                          c__PDAFcython.c__next_observation_pdaf,
                                          &outflag
                                         )

    return outflag

def omi_assimilate_global (py__collect_state_pdaf,
                           py__distribute_state_pdaf,
                           py__init_dim_obs_pdaf,
                           py__obs_op_pdaf,
                           py__prepoststep_pdaf,
                           py__next_observation_pdaf
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_global or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdafomi_assimilate_global (c__PDAFcython.c__collect_state_pdaf,
                                  c__PDAFcython.c__distribute_state_pdaf,
                                  c__PDAFcython.c__init_dim_obs_pdaf,
                                  c__PDAFcython.c__obs_op_pdaf,
                                  c__PDAFcython.c__prepoststep_pdaf,
                                  c__PDAFcython.c__next_observation_pdaf,
                                  &flag
                                 )

    return flag

def omi_assimilate_hyb3dvar_estkf (py__collect_state_pdaf,
                                   py__distribute_state_pdaf,
                                   py__init_dim_obs_pdaf,
                                   py__obs_op_pdaf,
                                   py__cvt_ens_pdaf,
                                   py__cvt_adj_ens_pdaf,
                                   py__cvt_pdaf,
                                   py__cvt_adj_pdaf,
                                   py__obs_op_lin_pdaf,
                                   py__obs_op_adj_pdaf,
                                   py__prepoststep_pdaf,
                                   py__next_observation_pdaf,
                                   int outflag
                                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_hyb3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply ensemble control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint ensemble control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdafomi_assimilate_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                          c__PDAFcython.c__distribute_state_pdaf,
                                          c__PDAFcython.c__init_dim_obs_pdaf,
                                          c__PDAFcython.c__obs_op_pdaf,
                                          c__PDAFcython.c__cvt_ens_pdaf,
                                          c__PDAFcython.c__cvt_adj_ens_pdaf,
                                          c__PDAFcython.c__cvt_pdaf,
                                          c__PDAFcython.c__cvt_adj_pdaf,
                                          c__PDAFcython.c__obs_op_lin_pdaf,
                                          c__PDAFcython.c__obs_op_adj_pdaf,
                                          c__PDAFcython.c__prepoststep_pdaf,
                                          c__PDAFcython.c__next_observation_pdaf,
                                          &outflag
                                         )

    return outflag

def omi_assimilate_hyb3dvar_lestkf (py__collect_state_pdaf,
                                    py__distribute_state_pdaf,
                                    py__init_dim_obs_f_pdaf,
                                    py__obs_op_f_pdaf,
                                    py__cvt_ens_pdaf,
                                    py__cvt_adj_ens_pdaf,
                                    py__cvt_pdaf,
                                    py__cvt_adj_pdaf,
                                    py__obs_op_lin_pdaf,
                                    py__obs_op_adj_pdaf,
                                    py__init_n_domains_p_pdaf,
                                    py__init_dim_l_pdaf,
                                    py__init_dim_obs_l_pdaf,
                                    py__g2l_state_pdaf,
                                    py__l2g_state_pdaf,
                                    py__prepoststep_pdaf,
                                    py__next_observation_pdaf,
                                    int outflag
                                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_hyb3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of full observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Full observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize local dimimension of obs. vector

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from local state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step, time and dimension of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    c__pdafomi_assimilate_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                           c__PDAFcython.c__distribute_state_pdaf,
                                           c__PDAFcython.c__init_dim_obs_f_pdaf,
                                           c__PDAFcython.c__obs_op_f_pdaf,
                                           c__PDAFcython.c__cvt_ens_pdaf,
                                           c__PDAFcython.c__cvt_adj_ens_pdaf,
                                           c__PDAFcython.c__cvt_pdaf,
                                           c__PDAFcython.c__cvt_adj_pdaf,
                                           c__PDAFcython.c__obs_op_lin_pdaf,
                                           c__PDAFcython.c__obs_op_adj_pdaf,
                                           c__PDAFcython.c__init_n_domains_p_pdaf,
                                           c__PDAFcython.c__init_dim_l_pdaf,
                                           c__PDAFcython.c__init_dim_obs_l_pdaf,
                                           c__PDAFcython.c__g2l_state_pdaf,
                                           c__PDAFcython.c__l2g_state_pdaf,
                                           c__PDAFcython.c__prepoststep_pdaf,
                                           c__PDAFcython.c__next_observation_pdaf,
                                           &outflag
                                          )

    return outflag

def omi_assimilate_lenkf (py__collect_state_pdaf,
                          py__distribute_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__prepoststep_pdaf,
                          py__localize_covar_pdaf,
                          py__next_observation_pdaf
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_lenkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__localize_covar_pdaf : Callable[dim_p:int, dim_obs:int, hp_p : ndarray[tuple[dim_obs, dim_p], np.float64], hph : ndarray[tuple[dim_obs, dim_obs], np.float64]]
        Apply localization to HP and HPH^T

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        dim_obs:int
            number of observations
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

        Returns
        -------
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdafomi_assimilate_lenkf (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__distribute_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 c__PDAFcython.c__localize_covar_pdaf,
                                 c__PDAFcython.c__next_observation_pdaf,
                                 &flag
                                )

    return flag

def omi_assimilate_local (py__collect_state_pdaf,
                          py__distribute_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__prepoststep_pdaf,
                          py__init_n_domains_p_pdaf,
                          py__init_dim_l_pdaf,
                          py__init_dim_obs_l_pdaf,
                          py__g2l_state_pdaf,
                          py__l2g_state_pdaf,
                          py__next_observation_pdaf
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_local or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdafomi_assimilate_local (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__distribute_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 c__PDAFcython.c__init_n_domains_p_pdaf,
                                 c__PDAFcython.c__init_dim_l_pdaf,
                                 c__PDAFcython.c__init_dim_obs_l_pdaf,
                                 c__PDAFcython.c__g2l_state_pdaf,
                                 c__PDAFcython.c__l2g_state_pdaf,
                                 c__PDAFcython.c__next_observation_pdaf,
                                 &flag
                                )

    return flag

def omi_generate_obs (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_f_pdaf,
                      py__obs_op_f_pdaf,
                      py__get_obs_f_pdaf,
                      py__prepoststep_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_generate_obs or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__distribute_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to distribute a state vector

        Parameters
        ----------
        dim_p:int
        state_p : ndarray[tuple[dim_p], np.float64]

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__get_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Provide observation vector to user

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__next_observation_pdaf : Callable[stepnow:int, nsteps:int, doexit:int, time:float]
        Provide time step and time of next observation

        Parameters
        ----------
        stepnow:int
            number of the current time step
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time

        Returns
        -------
        nsteps:int
            number of time steps until next obs
        doexit:int
            whether to exit forecasting (1 for exit)
        time:float
            current model (physical) time


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdafomi_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_f_pdaf,
                             c__PDAFcython.c__obs_op_f_pdaf,
                             c__PDAFcython.c__get_obs_f_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def omi_put_state_3dvar (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__cvt_pdaf,
                         py__cvt_adj_pdaf,
                         py__obs_op_lin_pdaf,
                         py__obs_op_adj_pdaf,
                         py__prepoststep_pdaf,
                         int outflag
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_3dvar or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    c__pdafomi_put_state_3dvar (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__init_dim_obs_pdaf,
                                c__PDAFcython.c__obs_op_pdaf,
                                c__PDAFcython.c__cvt_pdaf,
                                c__PDAFcython.c__cvt_adj_pdaf,
                                c__PDAFcython.c__obs_op_lin_pdaf,
                                c__PDAFcython.c__obs_op_adj_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                &outflag
                               )

    return outflag

def omi_put_state_en3dvar_estkf (py__collect_state_pdaf,
                                 py__init_dim_obs_pdaf,
                                 py__obs_op_pdaf,
                                 py__cvt_ens_pdaf,
                                 py__cvt_adj_ens_pdaf,
                                 py__obs_op_lin_pdaf,
                                 py__obs_op_adj_pdaf,
                                 py__prepoststep_pdaf,
                                 int outflag
                                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_en3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    c__pdafomi_put_state_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                        c__PDAFcython.c__init_dim_obs_pdaf,
                                        c__PDAFcython.c__obs_op_pdaf,
                                        c__PDAFcython.c__cvt_ens_pdaf,
                                        c__PDAFcython.c__cvt_adj_ens_pdaf,
                                        c__PDAFcython.c__obs_op_lin_pdaf,
                                        c__PDAFcython.c__obs_op_adj_pdaf,
                                        c__PDAFcython.c__prepoststep_pdaf,
                                        &outflag
                                       )

    return outflag

def omi_put_state_en3dvar_lestkf (py__collect_state_pdaf,
                                  py__init_dim_obs_f_pdaf,
                                  py__obs_op_f_pdaf,
                                  py__cvt_ens_pdaf,
                                  py__cvt_adj_ens_pdaf,
                                  py__obs_op_lin_pdaf,
                                  py__obs_op_adj_pdaf,
                                  py__init_n_domains_p_pdaf,
                                  py__init_dim_l_pdaf,
                                  py__init_dim_obs_l_pdaf,
                                  py__g2l_state_pdaf,
                                  py__l2g_state_pdaf,
                                  py__prepoststep_pdaf,
                                  int outflag
                                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_en3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of full observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Full observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize local dimimension of obs. vector

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from local state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    c__pdafomi_put_state_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                         c__PDAFcython.c__init_dim_obs_f_pdaf,
                                         c__PDAFcython.c__obs_op_f_pdaf,
                                         c__PDAFcython.c__cvt_ens_pdaf,
                                         c__PDAFcython.c__cvt_adj_ens_pdaf,
                                         c__PDAFcython.c__obs_op_lin_pdaf,
                                         c__PDAFcython.c__obs_op_adj_pdaf,
                                         c__PDAFcython.c__init_n_domains_p_pdaf,
                                         c__PDAFcython.c__init_dim_l_pdaf,
                                         c__PDAFcython.c__init_dim_obs_l_pdaf,
                                         c__PDAFcython.c__g2l_state_pdaf,
                                         c__PDAFcython.c__l2g_state_pdaf,
                                         c__PDAFcython.c__prepoststep_pdaf,
                                         &outflag
                                        )

    return outflag

def omi_put_state_generate_obs (py__collect_state_pdaf,
                                py__init_dim_obs_f_pdaf,
                                py__obs_op_f_pdaf,
                                py__get_obs_f_pdaf,
                                py__prepoststep_pdaf
                               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_generate_obs or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__get_obs_f_pdaf : Callable[step:int, dim_obs_f:int, observation_f : ndarray[tuple[dim_obs_f], np.float64]]
        Provide observation vector to user

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

        Returns
        -------
        observation_f : ndarray[tuple[dim_obs_f], np.float64]
            Full vector of synthetic observations (process-local)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__get_obs_f_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       &flag
                                      )

    return flag

def omi_put_state_global (py__collect_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__prepoststep_pdaf
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_global or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_global (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 &flag
                                )

    return flag

def omi_put_state_hyb3dvar_estkf (py__collect_state_pdaf,
                                  py__init_dim_obs_pdaf,
                                  py__obs_op_pdaf,
                                  py__cvt_ens_pdaf,
                                  py__cvt_adj_ens_pdaf,
                                  py__cvt_pdaf,
                                  py__cvt_adj_pdaf,
                                  py__obs_op_lin_pdaf,
                                  py__obs_op_adj_pdaf,
                                  py__prepoststep_pdaf,
                                  int outflag
                                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_hyb3dvar_estkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply ensemble control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint ensemble control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    c__pdafomi_put_state_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                         c__PDAFcython.c__init_dim_obs_pdaf,
                                         c__PDAFcython.c__obs_op_pdaf,
                                         c__PDAFcython.c__cvt_ens_pdaf,
                                         c__PDAFcython.c__cvt_adj_ens_pdaf,
                                         c__PDAFcython.c__cvt_pdaf,
                                         c__PDAFcython.c__cvt_adj_pdaf,
                                         c__PDAFcython.c__obs_op_lin_pdaf,
                                         c__PDAFcython.c__obs_op_adj_pdaf,
                                         c__PDAFcython.c__prepoststep_pdaf,
                                         &outflag
                                        )

    return outflag

def omi_put_state_hyb3dvar_lestkf (py__collect_state_pdaf,
                                   py__init_dim_obs_f_pdaf,
                                   py__obs_op_f_pdaf,
                                   py__cvt_ens_pdaf,
                                   py__cvt_adj_ens_pdaf,
                                   py__cvt_pdaf,
                                   py__cvt_adj_pdaf,
                                   py__obs_op_lin_pdaf,
                                   py__obs_op_adj_pdaf,
                                   py__init_n_domains_p_pdaf,
                                   py__init_dim_l_pdaf,
                                   py__init_dim_obs_l_pdaf,
                                   py__g2l_state_pdaf,
                                   py__l2g_state_pdaf,
                                   py__prepoststep_pdaf,
                                   int outflag
                                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_hyb3dvar_lestkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_f_pdaf : Callable[step:int, dim_obs_f:int]
        Initialize dimension of full observation vector

        Parameters
        ----------
        step:int
            Current time step
        dim_obs_f:int
            Size of the full observation vector

        Returns
        -------
        dim_obs_f:int
            Size of the full observation vector

    py__obs_op_f_pdaf : Callable[step:int, dim_p:int, dim_obs_f:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_f : ndarray[tuple[dim_obs_f], np.float64]]
        Full observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_f:int
            Size of full observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_f : ndarray[tuple[dim_obs_f], np.float64]
            Full observed state (i.e. the result after applying the observation operator to state_p)

    py__cvt_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cvec_ens:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], v_p : ndarray[tuple[dim_cvec_ens], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local dimension of state
        dim_ens:int
            Ensemble size
        dim_cvec_ens:int
            Dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        v_p : ndarray[tuple[dim_cvec_ens], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local state increment

    py__cvt_adj_ens_pdaf : Callable[iter:int, dim_p:int, dim_ens:int, dim_cv_ens_p:int, ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_ens:int
            Ensemble size
        dim_cv_ens_p:int
            PE-local dimension of control vector
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            PE-local ensemble
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local input vector
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cv_ens_p], np.float64]
            PE-local result vector

    py__cvt_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, cv_p : ndarray[tuple[dim_cvec], np.float64], Vv_p : ndarray[tuple[dim_p], np.float64]]
        Apply control vector transform matrix to control vector

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

        Returns
        -------
        Vv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)

    py__cvt_adj_pdaf : Callable[iter:int, dim_p:int, dim_cvec:int, Vcv_p : ndarray[tuple[dim_p], np.float64], cv_p : ndarray[tuple[dim_cvec], np.float64]]
        Apply adjoint control vector transform matrix

        Parameters
        ----------
        iter:int
            Iteration of optimization
        dim_p:int
            PE-local observation dimension
        dim_cvec:int
            Dimension of control vector
        Vcv_p : ndarray[tuple[dim_p], np.float64]
            PE-local result vector (state vector increment)
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

        Returns
        -------
        cv_p : ndarray[tuple[dim_cvec], np.float64]
            PE-local control vector

    py__obs_op_lin_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Linearized observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

    py__obs_op_adj_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Adjoint observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            PE-local dimension of state
        dim_obs_p:int
            Dimension of observed state
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            PE-local observed state

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            PE-local model state

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize local dimimension of obs. vector

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from local state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    c__pdafomi_put_state_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                          c__PDAFcython.c__init_dim_obs_f_pdaf,
                                          c__PDAFcython.c__obs_op_f_pdaf,
                                          c__PDAFcython.c__cvt_ens_pdaf,
                                          c__PDAFcython.c__cvt_adj_ens_pdaf,
                                          c__PDAFcython.c__cvt_pdaf,
                                          c__PDAFcython.c__cvt_adj_pdaf,
                                          c__PDAFcython.c__obs_op_lin_pdaf,
                                          c__PDAFcython.c__obs_op_adj_pdaf,
                                          c__PDAFcython.c__init_n_domains_p_pdaf,
                                          c__PDAFcython.c__init_dim_l_pdaf,
                                          c__PDAFcython.c__init_dim_obs_l_pdaf,
                                          c__PDAFcython.c__g2l_state_pdaf,
                                          c__PDAFcython.c__l2g_state_pdaf,
                                          c__PDAFcython.c__prepoststep_pdaf,
                                          &outflag
                                         )

    return outflag

def omi_put_state_lenkf (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__prepoststep_pdaf,
                         py__localize_covar_pdaf
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_lenkf or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__localize_covar_pdaf : Callable[dim_p:int, dim_obs:int, hp_p : ndarray[tuple[dim_obs, dim_p], np.float64], hph : ndarray[tuple[dim_obs, dim_obs], np.float64]]
        Apply localization to HP and HPH^T

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        dim_obs:int
            number of observations
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph

        Returns
        -------
        hp_p : ndarray[tuple[dim_obs, dim_p], np.float64]
            pe local part of matrix hp
        hph : ndarray[tuple[dim_obs, dim_obs], np.float64]
            matrix hph


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf

    cdef int flag

    c__pdafomi_put_state_lenkf (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__init_dim_obs_pdaf,
                                c__PDAFcython.c__obs_op_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__localize_covar_pdaf,
                                &flag
                               )

    return flag

def omi_put_state_local (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__prepoststep_pdaf,
                         py__init_n_domains_p_pdaf,
                         py__init_dim_l_pdaf,
                         py__init_dim_obs_l_pdaf,
                         py__g2l_state_pdaf,
                         py__l2g_state_pdaf
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_local or PDAF source files 

    Parameters
    ----------
    py__collect_state_pdaf : Callable[dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Routine to collect a state vector

        Parameters
        ----------
        dim_p:int
            pe-local state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            local state vector

    py__init_dim_obs_pdaf : Callable[step:int, dim_obs_p:int]
        Initialize dimension of observation vector

        Parameters
        ----------
        step:int
            current time step
        dim_obs_p:int
            dimension of observation vector

        Returns
        -------
        dim_obs_p:int
            dimension of observation vector

    py__obs_op_pdaf : Callable[step:int, dim_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], m_state_p : ndarray[tuple[dim_obs_p], np.float64]]
        Observation operator

        Parameters
        ----------
        step:int
            Current time step
        dim_p:int
            Size of state vector (local part in case of parallel decomposed state)
        dim_obs_p:int
            Size of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            Model state vector
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

        Returns
        -------
        m_state_p : ndarray[tuple[dim_obs_p], np.float64]
            Observed state vector (i.e. the result after applying the observation operator to state_p)

    py__prepoststep_pdaf : Callable[step:int, dim_p:int, dim_ens:int, dim_ens_p:int, dim_obs_p:int, state_p : ndarray[tuple[dim_p], np.float64], uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64], ens_p : ndarray[tuple[dim_p, dim_ens], np.float64], flag:int]
        User supplied pre/poststep routine

        Parameters
        ----------
        step:int
            current time step (negative for call after forecast)
        dim_p:int
            pe-local state dimension
        dim_ens:int
            size of state ensemble
        dim_ens_p:int
            pe-local size of ensemble
        dim_obs_p:int
            pe-local dimension of observation vector
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble
        flag:int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local forecast/analysis state(the array 'state_p' is not generally notinitialized in the case of seik.it can be used freely here.)
        uinv : ndarray[tuple[dim_ens-1, dim_ens-1], np.float64]
            inverse of matrix u
        ens_p : ndarray[tuple[dim_p, dim_ens], np.float64]
            pe-local state ensemble

    py__init_n_domains_p_pdaf : Callable[step:int, n_domains_p:int]
        Provide number of local analysis domains

        Parameters
        ----------
        step:int
            current time step
        n_domains_p:int
            pe-local number of analysis domains

        Returns
        -------
        n_domains_p:int
            pe-local number of analysis domains

    py__init_dim_l_pdaf : Callable[step:int, domain_p:int, dim_l:int]
        Init state dimension for local ana. domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension

        Returns
        -------
        dim_l:int
            local state dimension

    py__init_dim_obs_l_pdaf : Callable[domain_p:int, step:int, dim_obs_f:int, dim_obs_l:int]
        Initialize dim. of obs. vector for local ana. domain

        Parameters
        ----------
        domain_p:int
            index of current local analysis domain
        step:int
            current time step
        dim_obs_f:int
            full dimension of observation vector
        dim_obs_l:int
            local dimension of observation vector

        Returns
        -------
        dim_obs_l:int
            local dimension of observation vector

    py__g2l_state_pdaf : Callable[step:int, domain_p:int, dim_p:int, state_p : ndarray[tuple[dim_p], np.float64], dim_l:int, state_l : ndarray[tuple[dim_l], np.float64]]
        Get state on local ana. domain from full state

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

        Returns
        -------
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain

    py__l2g_state_pdaf : Callable[step:int, domain_p:int, dim_l:int, state_l : ndarray[tuple[dim_l], np.float64], dim_p:int, state_p : ndarray[tuple[dim_p], np.float64]]
        Init full state from state on local analysis domain

        Parameters
        ----------
        step:int
            current time step
        domain_p:int
            current local analysis domain
        dim_l:int
            local state dimension
        state_l : ndarray[tuple[dim_l], np.float64]
            state vector on local analysis domain
        dim_p:int
            pe-local full state dimension
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector

        Returns
        -------
        state_p : ndarray[tuple[dim_p], np.float64]
            pe-local full state vector


    Returns
    -------
    flag : int
        Status flag
    """
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf

    cdef int flag

    c__pdafomi_put_state_local (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__init_dim_obs_pdaf,
                                c__PDAFcython.c__obs_op_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__init_n_domains_p_pdaf,
                                c__PDAFcython.c__init_dim_l_pdaf,
                                c__PDAFcython.c__init_dim_obs_l_pdaf,
                                c__PDAFcython.c__g2l_state_pdaf,
                                c__PDAFcython.c__l2g_state_pdaf,
                                &flag
                               )

    return flag

def omi_init_obs_f_cb (int step,
                       int dim_obs_f,
                       double[::1] observation_f
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obs_f_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_f : int
        Dimension of full observation vector
    observation_f : ndarray[tuple[dim_obs_f], np.float64]
        Full observation vector

    Returns
    -------
    observation_f : ndarray[tuple[dim_obs_f], np.float64]
         Full observation vector
    """
    c__pdafomi_init_obs_f_cb (&step,
                              &dim_obs_f,
                              &observation_f[0]
                             )

    return np.asarray(observation_f).reshape((dim_obs_f), order='F')

def omi_init_obsvar_cb (int step,
                        int dim_obs_p,
                        double[::1] obs_p,
                        double meanvar
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obsvar_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
        PE-local observation vector
    meanvar : float
        Mean observation error variance

    Returns
    -------
    meanvar : float
        Mean observation error variance
    """
    c__pdafomi_init_obsvar_cb (&step,
                               &dim_obs_p,
                               &obs_p[0],
                               &meanvar
                              )

    return meanvar

def omi_g2l_obs_cb (int domain_p,
                    int step,
                    int dim_obs_f,
                    int dim_obs_l,
                    double[::1] ostate_f,
                    double[::1] ostate_l
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_g2l_obs_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain
    step : int
        Current time step
    dim_obs_f : int
        Dimension of full PE-local observation vector
    dim_obs_l : int
        Dimension of local observation vector
    ostate_f : ndarray[tuple[dim_obs_f], np.float64]
        Full PE-local obs.ervation vector
    ostate_l : ndarray[tuple[dim_obs_l], np.float64]
        Observation vector on local domain

    Returns
    -------
    ostate_l : ndarray[tuple[dim_obs_l], np.float64]
         Observation vector on local domain
    """
    c__pdafomi_g2l_obs_cb (&domain_p,
                           &step,
                           &dim_obs_f,
                           &dim_obs_l,
                           &ostate_f[0],
                           &ostate_l[0]
                          )

    return np.asarray(ostate_l).reshape((dim_obs_l), order='F')

def omi_init_obs_l_cb (int domain_p,
                       int step,
                       int dim_obs_l,
                       double[::1] observation_l
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obs_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain index
    step : int
        Current time step
    dim_obs_l : int
        Local dimension of observation vector
    observation_l : ndarray[tuple[dim_obs_l], np.float64]
        Local observation vector

    Returns
    -------
    observation_l : ndarray[tuple[dim_obs_l], np.float64]
         Local observation vector
    """
    c__pdafomi_init_obs_l_cb (&domain_p,
                              &step,
                              &dim_obs_l,
                              &observation_l[0]
                             )

    return np.asarray(observation_l).reshape((dim_obs_l), order='F')

def omi_init_obsvar_l_cb (int domain_p,
                          int step,
                          int dim_obs_l,
                          double[::1] obs_l,
                          double meanvar_l
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obsvar_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Local dimension of observation vector
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        Local observation vector
    meanvar_l : float
        Mean local observation error variance

    Returns
    -------
    meanvar_l : float
        Mean local observation error variance
    """
    c__pdafomi_init_obsvar_l_cb (&domain_p,
                                 &step,
                                 &dim_obs_l,
                                 &obs_l[0],
                                 &meanvar_l
                                )

    return meanvar_l

def omi_prodRinvA_l_cb (int domain_p,
                        int step,
                        int dim_obs_l,
                        int rank,
                        double[::1] obs_l,
                        double[::1,:] A_l,
                        double[::1,:] C_l
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Dimension of local observation vector
    rank : int
        Rank of initial covariance matrix
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        Local vector of observations
    A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
        Input matrix
    C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
        Output matrix

    Returns
    -------
    A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
         Input matrix
    C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
         Output matrix
    """
    cdef double[::1] A_l_f = np.asfortranarray(A_l).ravel(order="F")
    cdef double[::1] C_l_f = np.asfortranarray(C_l).ravel(order="F")
    c__pdafomi_prodrinva_l_cb (&domain_p,
                               &step,
                               &dim_obs_l,
                               &rank,
                               &obs_l[0],
                               &A_l_f[0],
                               &C_l_f[0]
                              )

    return np.asarray(A_l).reshape((dim_obs_l, rank), order='F'), np.asarray(C_l).reshape((dim_obs_l, rank), order='F')

def omi_likelihood_l_cb (int domain_p,
                         int step,
                         int dim_obs_l,
                         double[::1] obs_l,
                         double[::1] resid_l,
                         double lhood_l
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        PE-local dimension of obs. vector
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        PE-local vector of observations
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
        Input vector of residuum
    lhood_l : float
        Output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
         Input vector of residuum
    lhood_l : float
        Output vector - log likelihood
    """
    c__pdafomi_likelihood_l_cb (&domain_p,
                                &step,
                                &dim_obs_l,
                                &obs_l[0],
                                &resid_l[0],
                                &lhood_l
                               )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), lhood_l

def omi_prodRinvA_cb (int step,
                      int dim_obs_p,
                      int ncol,
                      double[::1] obs_p,
                      double[::1,:] A_p,
                      double[::1,:] C_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of PE-local observation vector
    ncol : int
        Number of columns in A_p and C_p
    obs_p : ndarray[tuple[dim_obs_p], np.float64]
        PE-local vector of observations
    A_p : ndarray[tuple[dim_obs_p, ncol], np.float64]
        Input matrix
    C_p : ndarray[tuple[dim_obs_p, ncol], np.float64]
        Output matrix

    Returns
    -------
    C_p : ndarray[tuple[dim_obs_p, ncol], np.float64]
         Output matrix
    """
    cdef double[::1] A_p_f = np.asfortranarray(A_p).ravel(order="F")
    cdef double[::1] C_p_f = np.asfortranarray(C_p).ravel(order="F")
    c__pdafomi_prodrinva_cb (&step,
                             &dim_obs_p,
                             &ncol,
                             &obs_p[0],
                             &A_p_f[0],
                             &C_p_f[0]
                            )

    return np.asarray(C_p).reshape((dim_obs_p, ncol), order='F')

def omi_likelihood_cb (int step,
                       int dim_obs,
                       double[::1] obs,
                       double[::1] resid,
                       double lhood
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs : int
        PE-local dimension of obs. vector
    obs : ndarray[tuple[dim_obs], np.float64]
        PE-local vector of observations
    resid : ndarray[tuple[dim_obs], np.float64]
        Input vector of residuum
    lhood : float
        Output vector - log likelihood

    Returns
    -------
    lhood : float
        Output vector - log likelihood
    """
    c__pdafomi_likelihood_cb (&step,
                              &dim_obs,
                              &obs[0],
                              &resid[0],
                              &lhood
                             )

    return lhood

def omi_add_obs_error_cb (int step,
                          int dim_obs_p,
                          double[::1,:] C_p
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_add_obs_error_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of PE-local observation vector
    C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
        Matrix to which R is added

    Returns
    -------
    C_p : ndarray[tuple[dim_obs_p, dim_obs_p], np.float64]
         Matrix to which R is added
    """
    cdef double[::1] C_p_f = np.asfortranarray(C_p).ravel(order="F")
    c__pdafomi_add_obs_error_cb (&step,
                                 &dim_obs_p,
                                 &C_p_f[0]
                                )

    return np.asarray(C_p).reshape((dim_obs_p, dim_obs_p), order='F')

def omi_init_obscovar_cb (int step,
                          int dim_obs,
                          int dim_obs_p,
                          double[::1,:] covar,
                          double[::1] m_state_p,
                          bint isdiag
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obscovar_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs : int
        Dimension of observation vector
    dim_obs_p : int
        PE-local dimension of obs. vector
    covar : ndarray[tuple[dim_obs, dim_obs], np.float64]
        Observation error covar. matrix
    m_state_p : ndarray[tuple[dim_obs_p], np.float64]
        Observation vector
    isdiag : bool
        Whether matrix R is diagonal

    Returns
    -------
    covar : ndarray[tuple[dim_obs, dim_obs], np.float64]
         Observation error covar. matrix
    isdiag : bool
        Whether matrix R is diagonal
    """
    cdef double[::1] covar_f = np.asfortranarray(covar).ravel(order="F")
    c__pdafomi_init_obscovar_cb (&step,
                                 &dim_obs,
                                 &dim_obs_p,
                                 &covar_f[0],
                                 &m_state_p[0],
                                 &isdiag
                                )

    return np.asarray(covar).reshape((dim_obs, dim_obs), order='F'), isdiag

def omi_init_obserr_f_cb (int step,
                          int dim_obs_f,
                          double[::1] obs_f,
                          double[::1] obserr_f
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obserr_f_cb or PDAF source files 

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_f : int
        Full dimension of observation vector
    obs_f : ndarray[tuple[dim_obs_f], np.float64]
        Full observation vector
    obserr_f : ndarray[tuple[dim_obs_f], np.float64]
        Full observation error stddev

    Returns
    -------
    obserr_f : ndarray[tuple[dim_obs_f], np.float64]
         Full observation error stddev
    """
    c__pdafomi_init_obserr_f_cb (&step,
                                 &dim_obs_f,
                                 &obs_f[0],
                                 &obserr_f[0]
                                )

    return np.asarray(obserr_f).reshape((dim_obs_f), order='F')

def omi_prodRinvA_hyb_l_cb (int domain_p,
                            int step,
                            int dim_obs_l,
                            int rank,
                            double[::1] obs_l,
                            double alpha,
                            double[::1,:] A_l,
                            double[::1,:] C_l
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_hyb_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Dimension of local observation vector
    rank : int
        Rank of initial covariance matrix
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        Local vector of observations
    alpha : float
        Hybrid weight
    A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
        Input matrix
    C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
        Output matrix

    Returns
    -------
    A_l : ndarray[tuple[dim_obs_l, rank], np.float64]
         Input matrix
    C_l : ndarray[tuple[dim_obs_l, rank], np.float64]
         Output matrix
    """
    cdef double[::1] A_l_f = np.asfortranarray(A_l).ravel(order="F")
    cdef double[::1] C_l_f = np.asfortranarray(C_l).ravel(order="F")
    c__pdafomi_prodrinva_hyb_l_cb (&domain_p,
                                   &step,
                                   &dim_obs_l,
                                   &rank,
                                   &obs_l[0],
                                   &alpha,
                                   &A_l_f[0],
                                   &C_l_f[0]
                                  )

    return np.asarray(A_l).reshape((dim_obs_l, rank), order='F'), np.asarray(C_l).reshape((dim_obs_l, rank), order='F')

def omi_likelihood_hyb_l_cb (int domain_p,
                             int step,
                             int dim_obs_l,
                             double[::1] obs_l,
                             double[::1] resid_l,
                             double alpha,
                             double lhood_l
                            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_hyb_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        PE-local dimension of obs. vector
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        PE-local vector of observations
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
        Input vector of residuum
    alpha : float
        Hybrid weight
    lhood_l : float
        Output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
         Input vector of residuum
    lhood_l : float
        Output vector - log likelihood
    """
    c__pdafomi_likelihood_hyb_l_cb (&domain_p,
                                    &step,
                                    &dim_obs_l,
                                    &obs_l[0],
                                    &resid_l[0],
                                    &alpha,
                                    &lhood_l
                                   )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), lhood_l

def omi_obsstats_l (int screen
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obsstats_l or PDAF source files 

    Parameters
    ----------
    screen : int
        Verbosity flag
    """

    c__pdafomi_obsstats_l (&screen
                          )

def omi_weights_l (int verbose,
                   int locweight,
                   double[::1] cradius,
                   double[::1] sradius,
                   double[:,:] matA,
                   double[::1] ivar_obs_l,
                   double[::1] dist_l
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_weights_l or PDAF source files 

    Parameters
    ----------
    verbose : int
        Verbosity flag
    locweight : int
        Localization weight type
    cradius : ndarray[tuple[nobs_l], np.float64]
        Localization cut-off radius
    sradius : ndarray[tuple[nobs_l], np.float64]
        support radius for weight functions
    matA : ndarray[tuple[nobs_l, ncols], np.float64]
    ivar_obs_l : ndarray[tuple[nobs_l], np.float64]
        Local vector of inverse obs. variances (nobs_l)
    dist_l : ndarray[tuple[nobs_l], np.float64]
        Local vector of obs. distances (nobs_l)

    Returns
    -------
    weight_l : ndarray[tuple[nobs_l], np.float64]
         Output: vector of weights
    """
    cdef double[::1] matA_f = np.asfortranarray(matA).ravel(order="F")
    cdef int nobs_l, ncols
    nobs_l = matA.shape[0]
    ncols = matA.shape[1]


    cdef double [::1] weight_l = np.zeros((nobs_l), dtype=np.float64).ravel()

    c__pdafomi_weights_l (&verbose,
                          &nobs_l,
                          &ncols,
                          &locweight,
                          &cradius[0],
                          &sradius[0],
                          &matA_f[0],
                          &ivar_obs_l[0],
                          &dist_l[0],
                          &weight_l[0]
                         )

    return np.asarray(weight_l).reshape((nobs_l), order='F')

def omi_weights_l_sgnl (int verbose,
                        int locweight,
                        double cradius,
                        double sradius,
                        double[:,:] matA,
                        double[::1] ivar_obs_l,
                        double[::1] dist_l
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_weights_l_sgnl or PDAF source files 

    Parameters
    ----------
    verbose : int
        Verbosity flag
    locweight : int
        Localization weight type
    cradius : float
        Localization cut-off radius
    sradius : float
        support radius for weight functions
    matA : ndarray[tuple[nobs_l, ncols], np.float64]
    ivar_obs_l : ndarray[tuple[nobs_l], np.float64]
        Local vector of inverse obs. variances (nobs_l)
    dist_l : ndarray[tuple[nobs_l], np.float64]
        Local vector of obs. distances (nobs_l)

    Returns
    -------
    weight_l : ndarray[tuple[nobs_l], np.float64]
         Output: vector of weights
    """
    cdef double[::1] matA_f = np.asfortranarray(matA).ravel(order="F")
    cdef int nobs_l, ncols
    nobs_l = matA.shape[0]
    ncols = matA.shape[1]


    cdef double [::1] weight_l = np.zeros((nobs_l), dtype=np.float64).ravel()

    c__pdafomi_weights_l_sgnl (&verbose,
                               &nobs_l,
                               &ncols,
                               &locweight,
                               &cradius,
                               &sradius,
                               &matA_f[0],
                               &ivar_obs_l[0],
                               &dist_l[0],
                               &weight_l[0]
                              )

    return np.asarray(weight_l).reshape((nobs_l), order='F')

def omi_check_error (int flag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_check_error or PDAF source files 

    Parameters
    ----------
    flag : int
        Error flag

    Returns
    -------
    flag : int
        Error flag
    """

    c__pdafomi_check_error (&flag
                           )

    return flag

def omi_gather_obsdims ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_gather_obsdims or PDAF source files 

    """
    c__pdafomi_gather_obsdims ()

def omi_obsstats (int screen
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obsstats or PDAF source files 

    Parameters
    ----------
    screen : int
        Verbosity flag
    """

    c__pdafomi_obsstats (&screen
                        )

def omi_init_dim_obs_l_iso (int i_obs,
                            double[::1] coords_l,
                            int locweight,
                            double cradius,
                            double sradius,
                            int cnt_obs_l
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_iso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord], np.float64]
        Coordinates of current analysis domain
    locweight : int
        Type of localization function
    cradius : float
        Localization cut-off radius (single or vector)
    sradius : float
        Support radius of localization function (single or vector)
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_iso (&i_obs,
                                   &ncoord,
                                   &coords_l[0],
                                   &locweight,
                                   &cradius,
                                   &sradius,
                                   &cnt_obs_l
                                  )

    return cnt_obs_l

def omi_init_dim_obs_l_noniso (int i_obs,
                               double[::1] coords_l,
                               int locweight,
                               double[::1] cradius,
                               double[::1] sradius,
                               int cnt_obs_l
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_noniso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord], np.float64]
        Coordinates of current analysis domain
    locweight : int
        Type of localization function
    cradius : ndarray[tuple[ncoord], np.float64]
        Vector of localization cut-off radii
    sradius : ndarray[tuple[ncoord], np.float64]
        Vector of support radii of localization function
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_noniso (&i_obs,
                                      &ncoord,
                                      &coords_l[0],
                                      &locweight,
                                      &cradius[0],
                                      &sradius[0],
                                      &cnt_obs_l
                                     )

    return cnt_obs_l

def omi_init_dim_obs_l_noniso_locweights (int i_obs,
                                          double[::1] coords_l,
                                          int[::1] locweights,
                                          double[::1] cradius,
                                          double[::1] sradius,
                                          int cnt_obs_l
                                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_noniso_locweights or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord], np.float64]
        Coordinates of current analysis domain
    locweights : ndarray[tuple[2], np.intc]
        Types of localization function
    cradius : ndarray[tuple[ncoord], np.float64]
        Vector of localization cut-off radii
    sradius : ndarray[tuple[ncoord], np.float64]
        Vector of support radii of localization function
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_noniso_locweights (&i_obs,
                                                 &ncoord,
                                                 &coords_l[0],
                                                 &locweights[0],
                                                 &cradius[0],
                                                 &sradius[0],
                                                 &cnt_obs_l
                                                )

    return cnt_obs_l

def omi_localize_covar_iso (int i_obs,
                            int locweight,
                            double cradius,
                            double sradius,
                            double[:,:] coords,
                            double[:,:] HP,
                            double[:,:] HPH
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_iso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observation type
    locweight : int
        Localization weight type
    cradius : float
        localization radius
    sradius : float
        support radius for weight functions
    coords : ndarray[tuple[ncoord, dim_p], np.float64]
        Coordinates of state vector elements
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
        Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
        Matrix HPH, dimension (nobs, nobs)

    Returns
    -------
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
         Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
         Matrix HPH, dimension (nobs, nobs)
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] HP_f = np.asfortranarray(HP).ravel(order="F")
    cdef double[::1] HPH_f = np.asfortranarray(HPH).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = HP.shape[0]
    _ = HP.shape[1]


    c__pdafomi_localize_covar_iso (&i_obs,
                                   &dim_p,
                                   &dim_obs,
                                   &ncoord,
                                   &locweight,
                                   &cradius,
                                   &sradius,
                                   &coords_f[0],
                                   &HP_f[0],
                                   &HPH_f[0]
                                  )

    return np.asarray(HP).reshape((dim_obs, dim_p), order='F'), np.asarray(HPH).reshape((dim_obs, dim_obs), order='F')

def omi_localize_covar_noniso (int i_obs,
                               int locweight,
                               double[::1] cradius,
                               double[::1] sradius,
                               double[:,:] coords,
                               double[:,:] HP,
                               double[:,:] HPH
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_noniso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        Data type with full observation
    locweight : int
        Localization weight type
    cradius : ndarray[tuple[ncoord], np.float64]
        Vector of localization cut-off radii
    sradius : ndarray[tuple[ncoord], np.float64]
        Vector of support radii of localization function
    coords : ndarray[tuple[ncoord, dim_p], np.float64]
        Coordinates of state vector elements
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
        Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
        Matrix HPH, dimension (nobs, nobs)

    Returns
    -------
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
         Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
         Matrix HPH, dimension (nobs, nobs)
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] HP_f = np.asfortranarray(HP).ravel(order="F")
    cdef double[::1] HPH_f = np.asfortranarray(HPH).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = HP.shape[0]
    _ = HP.shape[1]


    c__pdafomi_localize_covar_noniso (&i_obs,
                                      &dim_p,
                                      &dim_obs,
                                      &ncoord,
                                      &locweight,
                                      &cradius[0],
                                      &sradius[0],
                                      &coords_f[0],
                                      &HP_f[0],
                                      &HPH_f[0]
                                     )

    return np.asarray(HP).reshape((dim_obs, dim_p), order='F'), np.asarray(HPH).reshape((dim_obs, dim_obs), order='F')

def omi_localize_covar_noniso_locweights (int i_obs,
                                          int[::1] locweights,
                                          double[::1] cradius,
                                          double[::1] sradius,
                                          double[:,:] coords,
                                          double[:,:] HP,
                                          double[:,:] HPH
                                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_noniso_locweights or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observation type
    locweights : ndarray[tuple[2], np.intc]
        Types of localization function
    cradius : ndarray[tuple[ncoord], np.float64]
        Vector of localization cut-off radii
    sradius : ndarray[tuple[ncoord], np.float64]
        Vector of support radii of localization function
    coords : ndarray[tuple[ncoord, dim_p], np.float64]
        Coordinates of state vector elements
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
        Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
        Matrix HPH, dimension (nobs, nobs)

    Returns
    -------
    HP : ndarray[tuple[dim_obs, dim_p], np.float64]
         Matrix HP, dimension (nobs, dim)
    HPH : ndarray[tuple[dim_obs, dim_obs], np.float64]
         Matrix HPH, dimension (nobs, nobs)
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] HP_f = np.asfortranarray(HP).ravel(order="F")
    cdef double[::1] HPH_f = np.asfortranarray(HPH).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = HP.shape[0]
    _ = HP.shape[1]


    c__pdafomi_localize_covar_noniso_locweights (&i_obs,
                                                 &dim_p,
                                                 &dim_obs,
                                                 &ncoord,
                                                 &locweights[0],
                                                 &cradius[0],
                                                 &sradius[0],
                                                 &coords_f[0],
                                                 &HP_f[0],
                                                 &HPH_f[0]
                                                )

    return np.asarray(HP).reshape((dim_obs, dim_p), order='F'), np.asarray(HPH).reshape((dim_obs, dim_obs), order='F')

def omi_omit_by_inno_l_cb (int domain_p,
                           int dim_obs_l,
                           double[::1] resid_l,
                           double[::1] obs_l
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_omit_by_inno_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    dim_obs_l : int
        PE-local dimension of obs. vector
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
        Input vector of residuum
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
        Input vector of local observations

    Returns
    -------
    resid_l : ndarray[tuple[dim_obs_l], np.float64]
         Input vector of residuum
    obs_l : ndarray[tuple[dim_obs_l], np.float64]
         Input vector of local observations
    """
    c__pdafomi_omit_by_inno_l_cb (&domain_p,
                                  &dim_obs_l,
                                  &resid_l[0],
                                  &obs_l[0]
                                 )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), np.asarray(obs_l).reshape((dim_obs_l), order='F')

def omi_omit_by_inno_cb (int dim_obs_f,
                         double[::1] resid_f,
                         double[::1] obs_f
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_omit_by_inno_cb or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        Full dimension of obs. vector
    resid_f : ndarray[tuple[dim_obs_f], np.float64]
        Input vector of residuum
    obs_f : ndarray[tuple[dim_obs_f], np.float64]
        Input vector of full observations

    Returns
    -------
    resid_f : ndarray[tuple[dim_obs_f], np.float64]
         Input vector of residuum
    obs_f : ndarray[tuple[dim_obs_f], np.float64]
         Input vector of full observations
    """
    c__pdafomi_omit_by_inno_cb (&dim_obs_f,
                                &resid_f[0],
                                &obs_f[0]
                               )

    return np.asarray(resid_f).reshape((dim_obs_f), order='F'), np.asarray(obs_f).reshape((dim_obs_f), order='F')

