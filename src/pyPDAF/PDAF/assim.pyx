import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3


def get_state(int  steps, int  doexit, py__next_observation_pdaf,
    py__distribute_state_pdaf, py__prepoststep_pdaf, int  outflag):
    """get_state(steps:int, doexit:int, py__next_observation_pdaf:Callable, py__distribute_state_pdaf:Callable, py__prepoststep_pdaf:Callable, outflag:int) -> tuple[int, float, int, int]

    Distribute analysis state vector to an array.

    The primary purpose of this function is to distribute
    the analysis state vector to the model.
    This is attained by the user-supplied function
    :func:`py__distribute_state_pdaf`.
    One can also use this function to get the state vector
    for other purposes, e.g. to write the state vector to a file.

    In this function, the user-supplied function
    :func:`py__next_observation_pdaf` is executed
    to specify the number of forecast time steps
    until the next assimilation step.
    One can also use the user-supplied function to
    end the assimilation.

    In an online DA system, this function also execute
    the user-supplied function :func:`py__prepoststep_state_pdaf`,
    when this function is first called. The purpose of this design
    is to call this function right after :func:`pyPDAF.PDAF.init`
    to process the initial ensemble before using it to
    initialse model forecast. This user-supplied function
    will not be called afterwards.

    This function is also used in flexible parallel system
    where the number of ensemble members are greater than
    the parallel model tasks. In this case, this function
    is called multiple times to distribute the analysis ensemble.

    User-supplied function are executed in the following sequence:

        1. py__prepoststep_state_pdaf
           (only in online system when first called)
        2. py__distribute_state_pdaf
        3. py__next_observation_pdaf

    Parameters
    ----------
    steps : int
        Flag and number of time steps
    doexit : int
        Whether to exit from forecasts
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    outflag : int
        Status flag

    Returns
    -------
    steps : int
        Flag and number of time steps
    time : double
        current model time
    doexit : int
        Whether to exit from forecasts
    outflag : int
        Status flag
    """
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef double  time
    with nogil:
        c__pdaf_get_state(&steps, &time, &doexit,
                          pdaf_cb.c__next_observation_pdaf,
                          pdaf_cb.c__distribute_state_pdaf,
                          pdaf_cb.c__prepoststep_pdaf, &outflag)

    return steps, time, doexit, outflag


def assimilate_estkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__prodrinva_pdaf, py__init_obsvar_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`
    instead of this function.

    OMI functions need fewer user-supplied functions
    and improve DA efficiency.

    This function calls ESTKF
    (error space transform Kalman filter) [1]_.
    The ESTKF is a more efficient equivalent to the ETKF.

    The function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_estkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__init_obs_pdaf
        6. py__obs_op_pdaf (for each ensemble member)
        7. py__init_obsvar_pdaf (only relevant
           for adaptive forgetting factor schemes)
        8. py__prodRinvA_pdaf
        9. core DA algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_estkf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__prodrinva_pdaf,
                                 pdaf_cb.c__init_obsvar_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_estkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__prodrinva_pdaf,
    py__init_obsvar_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_estkf(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__prodrinva_pdaf,
                                    pdaf_cb.c__init_obsvar_pdaf, &outflag)

    return outflag


def assimilate_3dvar(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_3dvar`
    or :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DVar DA for a single step without OMI.
    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable
    transformation. This is a deterministic filtering
    scheme so no ensemble and
    parallelisation is needed.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_3dvar`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. core DA algorithm
        7. py__cvt_pdaf
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_3dvar`
       and :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_3dvar(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__prodrinva_pdaf,
                                 pdaf_cb.c__cvt_pdaf,
                                 pdaf_cb.c__cvt_adj_pdaf,
                                 pdaf_cb.c__obs_op_lin_pdaf,
                                 pdaf_cb.c__obs_op_adj_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_3dvar(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_3dvar(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__prodrinva_pdaf,
                                    pdaf_cb.c__cvt_pdaf,
                                    pdaf_cb.c__cvt_adj_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assimilate_en3dvar_lestkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf, py__init_obs_f_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf` or
    :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF.
    The background error covariance matrix is
    estimated by ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An LESTKF is used to generate ensemble perturbations.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_en3dvar_lestkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core DA algorithm
        7. py__cvt_ens_pdaf
        8. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__init_obs_pdaf
               (if global adaptive forgetting factor is used
               `type_forget=1` in :func:`pyPDAF.PDAF.init`)
            5. py__init_obsvar_pdaf
               (if global adaptive forgetting factor is used)
            6. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__g2l_state_pdaf
                4. py__g2l_obs_pdaf
                   (localise mean ensemble in observation space)
                5. py__init_obs_l_pdaf
                6. py__g2l_obs_pdaf
                   (localise each ensemble member
                   in observation space)
                7. py__init_obsvar_l_pdaf
                   (only called if local
                   adaptive forgetting factor
                   `type_forget=2` is used)
                8. py__prodRinvA_l_pdaf
                9. core DA algorithm
                10. py__l2g_state_pdaf
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`
       and
       :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
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
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_en3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__distribute_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__init_obs_pdaf,
                                          pdaf_cb.c__prodrinva_pdaf,
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
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__next_observation_pdaf,
                                          &outflag)

    return outflag


def assim_offline_en3dvar_lestkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf, py__init_obs_f_pdaf,
    py__init_obs_l_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
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
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_en3dvar_lestkf(pdaf_cb.c__init_dim_obs_pdaf,
                                             pdaf_cb.c__obs_op_pdaf,
                                             pdaf_cb.c__init_obs_pdaf,
                                             pdaf_cb.c__prodrinva_pdaf,
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
                                             pdaf_cb.c__prepoststep_pdaf,
                                             &outflag)

    return outflag


def assimilate_ensrf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obsvars_pdaf, py__localize_covar_serial_pdaf,
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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obsvars_pdaf = <void*>py__init_obsvars_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_ensrf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__init_obsvars_pdaf,
                                 pdaf_cb.c__localize_covar_serial_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_ensrf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obsvars_pdaf,
    py__localize_covar_serial_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obsvars_pdaf = <void*>py__init_obsvars_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_ensrf(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__init_obsvars_pdaf,
                                    pdaf_cb.c__localize_covar_serial_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assimilate_lknetf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prepoststep_pdaf, py__prodrinva_l_pdaf,
    py__prodrinva_hyb_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__g2l_obs_pdaf, py__init_obsvar_pdaf,
    py__init_obsvar_l_pdaf, py__likelihood_l_pdaf,
    py__likelihood_hyb_l_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate`
    or :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    A hybridised LETKF and LNETF [1]_ for a single DA step.
    The LNETF computes the distribution up to
    the second moment similar to Kalman filters but
    using a nonlinear weighting similar to
    particle filters. This leads to an equal weights
    assumption for the prior ensemble.
    The hybridisation with LETKF is expected to lead to
    improved performance for quasi-Gaussian problems.
    The function should be called at each model step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lknetf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf
           (for each ensemble member)
        6. py__init_obs_pdaf
           (if global adaptive forgetting factor `type_forget=1`
           is used in :func:`pyPDAF.PDAF.init`)
        7. py__init_obsvar_pdaf (if global adaptive
           forgetting factor is used)
        8. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__g2l_state_pdaf
            4. py__g2l_obs_pdaf
               (localise each ensemble member
               in observation space)
            5. py__init_obs_l_pdaf
            6. py__init_obsvar_l_pdaf
               (only called if local adaptive forgetting
               factor `type_forget=2` is used)
            7. py__prodRinvA_pdaf
            8. py__likelihood_l_pdaf
            9. core DA algorithm
            10. py__l2g_state_pdaf
        9. py__obs_op_pdaf
           (only called with `HKN` and `HNK` options
           called for each ensemble member)
        10. py__likelihood_hyb_l_pda
        11. py__init_obsvar_l_pdaf
            (only called if local adaptive forgetting
            factor `type_forget=2` is used)
        12. py__prodRinvA_hyb_l_pdaf
        13. py__prepoststep_state_pdaf
        14. py__distribute_state_pdaf
        15. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate`
       and :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`

    References
    ----------
    .. [1] Nerger, L.. (2022)
           Data assimilation for nonlinear systems with
           a hybrid nonlinear Kalman ensemble transform filter.
           Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221

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

    py__prodrinva_hyb_l_pdaf : Callable
        Provide product R^-1 A on local analysis domain with hybrid weight

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
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
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_lknetf(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__init_obs_pdaf,
                                  pdaf_cb.c__init_obs_l_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__prodrinva_l_pdaf,
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
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_lknetf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obs_l_pdaf, py__prepoststep_pdaf,
    py__prodrinva_l_pdaf, py__prodrinva_hyb_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__likelihood_l_pdaf, py__likelihood_hyb_l_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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

    py__prodrinva_hyb_l_pdaf : Callable
        Provide product R^-1 A on local analysis domain with hybrid weight

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
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
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_lknetf(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__init_obs_pdaf,
                                     pdaf_cb.c__init_obs_l_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf,
                                     pdaf_cb.c__prodrinva_l_pdaf,
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
                                     pdaf_cb.c__likelihood_hyb_l_pdaf, &outflag)

    return outflag


def assimilate_hyb3dvar_estkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_obsvar_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf`
    or :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    ESTKF in this implementation.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_hyb3dvar_estkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. the iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core 3DEnVar algorithm
        7. py__cvt_pdaf
        8. py__cvt_ens_pdaf
        9. Perform ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf
               (for ensemble mean)
            3. py__init_obs_pdaf
            4. py__obs_op_pdaf
               (for each ensemble member)
            5. py__init_obsvar_pdaf
               (only relevant for adaptive
               forgetting factor schemes)
            6. py__prodRinvA_pdaf
            7. core ESTKF algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf` and
       :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`.

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_hyb3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__distribute_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__init_obs_pdaf,
                                          pdaf_cb.c__prodrinva_pdaf,
                                          pdaf_cb.c__cvt_ens_pdaf,
                                          pdaf_cb.c__cvt_adj_ens_pdaf,
                                          pdaf_cb.c__cvt_pdaf,
                                          pdaf_cb.c__cvt_adj_pdaf,
                                          pdaf_cb.c__obs_op_lin_pdaf,
                                          pdaf_cb.c__obs_op_adj_pdaf,
                                          pdaf_cb.c__init_obsvar_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__next_observation_pdaf,
                                          &outflag)

    return outflag


def assim_offline_hyb3dvar_estkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_obsvar_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_hyb3dvar_estkf(pdaf_cb.c__init_dim_obs_pdaf,
                                             pdaf_cb.c__obs_op_pdaf,
                                             pdaf_cb.c__init_obs_pdaf,
                                             pdaf_cb.c__prodrinva_pdaf,
                                             pdaf_cb.c__cvt_pdaf,
                                             pdaf_cb.c__cvt_adj_pdaf,
                                             pdaf_cb.c__cvt_ens_pdaf,
                                             pdaf_cb.c__cvt_adj_ens_pdaf,
                                             pdaf_cb.c__obs_op_lin_pdaf,
                                             pdaf_cb.c__obs_op_adj_pdaf,
                                             pdaf_cb.c__init_obsvar_pdaf,
                                             pdaf_cb.c__prepoststep_pdaf,
                                             &outflag)

    return outflag


def assimilate_hyb3dvar_lestkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_dim_obs_f_pdaf,
    py__obs_op_f_pdaf, py__init_obs_f_pdaf, py__init_obs_l_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf, int  outflag):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf` or
    :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    LESTKF in this implementation.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_hyb3dvar_lestkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core DA algorithm
        7. py__cvt_pdaf
        8. py__cvt_ens_pdaf
        9. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__init_obs_pdaf
               (if global adaptive forgetting factor
               `type_forget=1` in :func:`pyPDAF.PDAF.init`)
            5. py__init_obsvar_pdaf
               (if global adaptive forgetting factor is used)
            6. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__g2l_state_pdaf
                4. py__g2l_obs_pdaf
                   (localise mean ensemble in observation space)
                5. py__init_obs_l_pdaf
                6. py__g2l_obs_pdaf
                   (localise each ensemble member
                   in observation space)
                7. py__init_obsvar_l_pdaf
                   (only called if local adaptive forgetting
                   factor `type_forget=2` is used)
                8. py__prodRinvA_l_pdaf
                9. core DA algorithm
                10. py__l2g_state_pdaf
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`
       and
       :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`

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
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
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
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf_assimilate_hyb3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__init_obs_pdaf,
                                           pdaf_cb.c__prodrinva_pdaf,
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
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assim_offline_hyb3dvar_lestkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_dim_obs_f_pdaf,
    py__obs_op_f_pdaf, py__init_obs_f_pdaf, py__init_obs_l_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
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
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_hyb3dvar_lestkf(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__init_obs_pdaf,
                                              pdaf_cb.c__prodrinva_pdaf,
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
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assimilate_lestkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prepoststep_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate`
    or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.

    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    Local ESTKF (error space transform Kalman filter) [1]_ for a single DA step without OMI.
    The LESTKF is a more efficient equivalent to the LETKF.

    This function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lestkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__init_obs_pdaf
           (if global adaptive forgetting factor
           `type_forget=1` is used
           in :func:`pyPDAF.PDAF.init`)
        7. py__init_obsvar_pdaf (if global adaptive
           forgetting factor is used)
        8. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__g2l_state_pdaf
            4. py__g2l_obs_pdaf (localise mean ensemble
               in observation space)
            5. py__init_obs_l_pdaf
            6. py__g2l_obs_pdaf
               (localise each ensemble member
               in observation space)
            7. py__init_obsvar_l_pdaf
               (only called if local adaptive
               forgetting factor `type_forget=2` is used)
            8. py__prodRinvA_l_pdaf
            9. core DA algorithm
            10. py__l2g_state_pdaf
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate`
       and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_lestkf(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__init_obs_pdaf,
                                  pdaf_cb.c__init_obs_l_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__prodrinva_l_pdaf,
                                  pdaf_cb.c__init_n_domains_p_pdaf,
                                  pdaf_cb.c__init_dim_l_pdaf,
                                  pdaf_cb.c__init_dim_obs_l_pdaf,
                                  pdaf_cb.c__g2l_state_pdaf,
                                  pdaf_cb.c__l2g_state_pdaf,
                                  pdaf_cb.c__g2l_obs_pdaf,
                                  pdaf_cb.c__init_obsvar_pdaf,
                                  pdaf_cb.c__init_obsvar_l_pdaf,
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_lestkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obs_l_pdaf, py__prepoststep_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_lestkf(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__init_obs_pdaf,
                                     pdaf_cb.c__init_obs_l_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf,
                                     pdaf_cb.c__prodrinva_l_pdaf,
                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                     pdaf_cb.c__init_dim_l_pdaf,
                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                     pdaf_cb.c__g2l_state_pdaf,
                                     pdaf_cb.c__l2g_state_pdaf,
                                     pdaf_cb.c__g2l_obs_pdaf,
                                     pdaf_cb.c__init_obsvar_pdaf,
                                     pdaf_cb.c__init_obsvar_l_pdaf, &outflag)

    return outflag


def assimilate_enkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__add_obs_err_pdaf, py__init_obs_covar_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Stochastic EnKF (ensemble Kalman filter) [1]_ for a single DA step without OMI.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_enkf` and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__add_obs_err_pdaf
        6. py__init_obs_pdaf
        7. py__init_obscovar_pdaf
        8. py__obs_op_pdaf (for each ensemble member)
        9. core DA algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`

    References
    ----------
    .. [1] Evensen, G. (1994),
           Sequential data assimilation with a nonlinear
           quasi-geostrophic model
           using Monte Carlo methods to forecast error statistics,
           J. Geophys. Res., 99(C5), 10143–10162,
           doi:10.1029/94JC00572.

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

    py__add_obs_err_pdaf : Callable
        Add obs error covariance R to HPH in EnKF

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
        Initialize obs. error cov. matrix R in EnKF

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_enkf(pdaf_cb.c__collect_state_pdaf,
                                pdaf_cb.c__distribute_state_pdaf,
                                pdaf_cb.c__init_dim_obs_pdaf,
                                pdaf_cb.c__obs_op_pdaf,
                                pdaf_cb.c__init_obs_pdaf,
                                pdaf_cb.c__prepoststep_pdaf,
                                pdaf_cb.c__add_obs_err_pdaf,
                                pdaf_cb.c__init_obs_covar_pdaf,
                                pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_enkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__add_obs_err_pdaf,
    py__init_obs_covar_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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

    py__add_obs_err_pdaf : Callable
        Add obs error covariance R to HPH in EnKF

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
        Initialize obs. error cov. matrix R in EnKF

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



    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_enkf(pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__init_obs_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__add_obs_err_pdaf,
                                   pdaf_cb.c__init_obs_covar_pdaf, &outflag)

    return outflag


def assimilate_letkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prepoststep_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate`
    or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.

    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    Local ensemble transform Kalman filter (LETKF) [1]_ for a single DA step without OMI.
    Implementation is based on [2]_.
    Note that the LESTKF is a more efficient equivalent
    to the LETKF.

    This function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_letkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__init_obs_pdaf
           (if global adaptive forgetting factor
           `type_forget=1` is used
           in :func:`pyPDAF.PDAF.init`)
        7. py__init_obsvar_pdaf (if global adaptive
           forgetting factor is used)
        8. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__g2l_state_pdaf
            4. py__g2l_obs_pdaf (localise mean ensemble
               in observation space)
            5. py__init_obs_l_pdaf
            6. py__g2l_obs_pdaf (localise each ensemble member
               in observation space)
            7. py__init_obsvar_l_pdaf
               (only called if local adaptive forgetting factor
               `type_forget=2` is used)
            8. py__prodRinvA_l_pdaf
            9. core DA algorithm
            10. py__l2g_state_pdaf
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate`
       and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`

    References
    ----------
    .. [1] Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007).
           Efficient data assimilation for spatiotemporal chaos:
           A local ensemble transform Kalman filter.
           Physica D: Nonlinear Phenomena, 230(1-2), 112-126.
    .. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_letkf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__init_obs_l_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__prodrinva_l_pdaf,
                                 pdaf_cb.c__init_n_domains_p_pdaf,
                                 pdaf_cb.c__init_dim_l_pdaf,
                                 pdaf_cb.c__init_dim_obs_l_pdaf,
                                 pdaf_cb.c__g2l_state_pdaf,
                                 pdaf_cb.c__l2g_state_pdaf,
                                 pdaf_cb.c__g2l_obs_pdaf,
                                 pdaf_cb.c__init_obsvar_pdaf,
                                 pdaf_cb.c__init_obsvar_l_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_letkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obs_l_pdaf, py__prepoststep_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_letkf(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__init_obs_l_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__prodrinva_l_pdaf,
                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                    pdaf_cb.c__init_dim_l_pdaf,
                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                    pdaf_cb.c__g2l_state_pdaf,
                                    pdaf_cb.c__l2g_state_pdaf,
                                    pdaf_cb.c__g2l_obs_pdaf,
                                    pdaf_cb.c__init_obsvar_pdaf,
                                    pdaf_cb.c__init_obsvar_l_pdaf, &outflag)

    return outflag


def assimilate_seik(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__prodrinva_pdaf, py__init_obsvar_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    This function will use singular evolutive
    interpolated Kalman filter [1]_ for a single DA step.
    The function should be called at each model step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_seik`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__init_obs_pdaf
        6. py__obs_op_pdaf (for each ensemble member)
        7. py__init_obsvar_pdaf (only relevant for
           adaptive forgetting factor schemes)
        8. py__prodRinvA_pdaf
        9. core DA algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`

    References
    ----------
    .. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).
           A singular evolutive extended Kalman
           filter for data assimilation
           in oceanography. Journal of Marine systems,
           16(3-4), 323-340.

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_seik(pdaf_cb.c__collect_state_pdaf,
                                pdaf_cb.c__distribute_state_pdaf,
                                pdaf_cb.c__init_dim_obs_pdaf,
                                pdaf_cb.c__obs_op_pdaf,
                                pdaf_cb.c__init_obs_pdaf,
                                pdaf_cb.c__prepoststep_pdaf,
                                pdaf_cb.c__prodrinva_pdaf,
                                pdaf_cb.c__init_obsvar_pdaf,
                                pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_seik(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__prodrinva_pdaf,
    py__init_obsvar_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_seik(pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__init_obs_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__prodrinva_pdaf,
                                   pdaf_cb.c__init_obsvar_pdaf, &outflag)

    return outflag


def assimilate_lnetf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prepoststep_pdaf, py__likelihood_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate`
    or :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Local Nonlinear Ensemble Transform Filter (LNETF) [1]_
    for a single DA step.
    The nonlinear filter computes the distribution up to
    the second moment similar to Kalman filters but
    it uses a nonlinear weighting similar to
    particle filters. This leads to an equal weights assumption
    for the prior ensemble at each step.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lnetf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__g2l_state_pdaf
            4. py__init_obs_l_pdaf
            5. py__g2l_obs_pdaf (localise each ensemble
               member in observation space)
            6. py__likelihood_l_pdaf
            7. core DA algorithm
            8. py__l2g_state_pdaf
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate`
       and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_lnetf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__init_obs_l_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__likelihood_l_pdaf,
                                 pdaf_cb.c__init_n_domains_p_pdaf,
                                 pdaf_cb.c__init_dim_l_pdaf,
                                 pdaf_cb.c__init_dim_obs_l_pdaf,
                                 pdaf_cb.c__g2l_state_pdaf,
                                 pdaf_cb.c__l2g_state_pdaf,
                                 pdaf_cb.c__g2l_obs_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_lnetf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obs_l_pdaf, py__prepoststep_pdaf,
    py__likelihood_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_lnetf(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__init_obs_l_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__likelihood_l_pdaf,
                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                    pdaf_cb.c__init_dim_l_pdaf,
                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                    pdaf_cb.c__g2l_state_pdaf,
                                    pdaf_cb.c__l2g_state_pdaf,
                                    pdaf_cb.c__g2l_obs_pdaf, &outflag)

    return outflag


def assimilate_prepost(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is used to preprocess and postprocess of the ensemble.

    No DA is performed in this function.
    Compared to :func:`pyPDAF.PDAF.prepost`,
    this function sets assimilation flag,
    which means that it is acted as an assimilation in PDAF.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_prepost`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf (preprocess, step < 0)
        3. py__prepoststep_state_pdaf (postprocess, step > 0)
        4. py__distribute_state_pdaf
        5. py__next_observation_pdaf

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
        c__pdaf_assimilate_prepost(pdaf_cb.c__collect_state_pdaf,
                                   pdaf_cb.c__distribute_state_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_en3dvar_estkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_obsvar_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`
    or :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step.
    The background error covariance matrix is estimated
    by an ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An ESTKF is used along with 3DEnVar
    to generate ensemble perturbations.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_en3dvar_estkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core 3DEnVar algorithm
        7. py__cvt_ens_pdaf
        8. ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf (for ensemble mean)
            3. py__init_obs_pdaf
            4. py__obs_op_pdaf (for each ensemble member)
            5. py__init_obsvar_pdaf
               (only relevant for adaptive
               forgetting factor schemes)
            6. py__prodRinvA_pdaf
            7. core ESTKF algorithm
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`
       and :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_en3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                         pdaf_cb.c__distribute_state_pdaf,
                                         pdaf_cb.c__init_dim_obs_pdaf,
                                         pdaf_cb.c__obs_op_pdaf,
                                         pdaf_cb.c__init_obs_pdaf,
                                         pdaf_cb.c__prodrinva_pdaf,
                                         pdaf_cb.c__cvt_ens_pdaf,
                                         pdaf_cb.c__cvt_adj_ens_pdaf,
                                         pdaf_cb.c__obs_op_lin_pdaf,
                                         pdaf_cb.c__obs_op_adj_pdaf,
                                         pdaf_cb.c__init_obsvar_pdaf,
                                         pdaf_cb.c__prepoststep_pdaf,
                                         pdaf_cb.c__next_observation_pdaf,
                                         &outflag)

    return outflag


def assim_offline_en3dvar_estkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_obsvar_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_en3dvar_estkf(pdaf_cb.c__init_dim_obs_pdaf,
                                            pdaf_cb.c__obs_op_pdaf,
                                            pdaf_cb.c__init_obs_pdaf,
                                            pdaf_cb.c__prodrinva_pdaf,
                                            pdaf_cb.c__cvt_ens_pdaf,
                                            pdaf_cb.c__cvt_adj_ens_pdaf,
                                            pdaf_cb.c__obs_op_lin_pdaf,
                                            pdaf_cb.c__obs_op_adj_pdaf,
                                            pdaf_cb.c__init_obsvar_pdaf,
                                            pdaf_cb.c__prepoststep_pdaf,
                                            &outflag)

    return outflag


def assimilate_netf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__likelihood_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    This function will use Nonlinear Ensemble
    Transform Filter (NETF) [1]_
    for a single DA step. The nonlinear filter
    computes the distribution up to
    the second moment similar to KF but using
    a nonlinear weighting similar to
    particle filter. This leads to an equal
    weights assumption for prior ensemble.
    The function should be called at each model step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_netf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__init_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__likelihood_pdaf
        7. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.

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
        Init. observation vector on local analysis domain

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_netf(pdaf_cb.c__collect_state_pdaf,
                                pdaf_cb.c__distribute_state_pdaf,
                                pdaf_cb.c__init_dim_obs_pdaf,
                                pdaf_cb.c__obs_op_pdaf,
                                pdaf_cb.c__init_obs_pdaf,
                                pdaf_cb.c__prepoststep_pdaf,
                                pdaf_cb.c__likelihood_pdaf,
                                pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_netf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__likelihood_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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
        Init. observation vector on local analysis domain

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_netf(pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__init_obs_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__likelihood_pdaf, &outflag)

    return outflag


def assimilate_pf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__likelihood_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    This function will use particle filter for a single DA step.
    This is a fully nonlinear filter, and may require
    a high number of ensemble members.
    A review of particle filter can be found at [1]_.
    The function should be called at each model step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_pf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__init_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__likelihood_pdaf
        7. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`

    References
    ----------
    .. [1] Van Leeuwen, P. J., Künsch, H. R.,
           Nerger, L., Potthast, R., & Reich, S. (2019).
           Particle filters for high‐dimensional geoscience
           applications:
           A review. Quarterly Journal of the Royal
           Meteorological Society, 145(723), 2335-2365.

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
        Init. observation vector on local analysis domain

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_pf(pdaf_cb.c__collect_state_pdaf,
                              pdaf_cb.c__distribute_state_pdaf,
                              pdaf_cb.c__init_dim_obs_pdaf,
                              pdaf_cb.c__obs_op_pdaf,
                              pdaf_cb.c__init_obs_pdaf,
                              pdaf_cb.c__prepoststep_pdaf,
                              pdaf_cb.c__likelihood_pdaf,
                              pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_pf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__likelihood_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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
        Init. observation vector on local analysis domain

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_pf(pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__likelihood_pdaf, &outflag)

    return outflag


def assimilate_lenkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__localize_covar_pdaf, py__add_obs_err_pdaf,
    py__init_obs_covar_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_lenkf`
    or :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Stochastic EnKF (ensemble Kalman filter)
    with covariance localisation [1]_
    for a single DA step without OMI.

    This is the only scheme for covariance localisation in PDAF.

    This function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lenkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__localize_pdaf
        6. py__add_obs_err_pdaf
        7. py__init_obs_pdaf
        8. py__init_obscovar_pdaf
        9. py__obs_op_pdaf (repeated to reduce storage)
        10. core DA algorith
        11. py__prepoststep_state_pdaf
        12. py__distribute_state_pdaf
        13. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_lenkf`
       and :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter
           Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

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

    py__add_obs_err_pdaf : Callable
        Add obs error covariance R to HPH in EnKF

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
        Initialize obs. error cov. matrix R in EnKF

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_lenkf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__localize_covar_pdaf,
                                 pdaf_cb.c__add_obs_err_pdaf,
                                 pdaf_cb.c__init_obs_covar_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_lenkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__localize_covar_pdaf,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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

    py__add_obs_err_pdaf : Callable
        Add obs error covariance R to HPH in EnKF

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
        Initialize obs. error cov. matrix R in EnKF

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



    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_lenkf(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__localize_covar_pdaf,
                                    pdaf_cb.c__add_obs_err_pdaf,
                                    pdaf_cb.c__init_obs_covar_pdaf, &outflag)

    return outflag


def assimilate_etkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__prepoststep_pdaf, py__prodrinva_pdaf, py__init_obsvar_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_global`
    or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`.

    PDAFlocal-OMI modules require fewer
    user-supplied functions and improved efficiency.

    Using ETKF (ensemble transform
    Kalman filter) [1]_ for a single DA step without OMI.
    The implementation is baed on [2]_.

    This function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_etkf` and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__init_obs_pdaf
        6. py__obs_op_pdaf (for each ensemble member)
        7. py__init_obsvar_pdaf (only relevant for
           adaptive forgetting factor schemes)
        8. py__prodRinvA_pdaf
        9. core DA algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_global`
       and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`

    References
    ----------
    .. [1] Bishop, C. H., B. J. Etherton, and S. J. Majumdar (2001)
           Adaptive Sampling with the Ensemble
           Transform Kalman Filter.
           Part I: Theoretical Aspects. Mon. Wea. Rev.,
           129, 420–436,
           doi: 10.1175/1520-0493(2001)129<0420:ASWTET>2.0.CO;2.
    .. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_etkf(pdaf_cb.c__collect_state_pdaf,
                                pdaf_cb.c__distribute_state_pdaf,
                                pdaf_cb.c__init_dim_obs_pdaf,
                                pdaf_cb.c__obs_op_pdaf,
                                pdaf_cb.c__init_obs_pdaf,
                                pdaf_cb.c__prepoststep_pdaf,
                                pdaf_cb.c__prodrinva_pdaf,
                                pdaf_cb.c__init_obsvar_pdaf,
                                pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_etkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__prepoststep_pdaf, py__prodrinva_pdaf,
    py__init_obsvar_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_etkf(pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__init_obs_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__prodrinva_pdaf,
                                   pdaf_cb.c__init_obsvar_pdaf, &outflag)

    return outflag


def assimilate_lseik(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_obs_pdaf,
    py__init_obs_l_pdaf, py__prepoststep_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_assimilate`
    or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.

    PDAF-OMI modules require fewer user-supplied
    functions and improved efficiency.

    Local singular evolutive interpolated Kalman filter [1]_
    for a single DA step.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lseik` and :func:`pyPDAF.PDAF.get_state`

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__init_obs_pdaf
           (if global adaptive forgetting factor `type_forget=1`
           is used in :func:`pyPDAF.PDAF.init`)
        7. py__init_obsvar_pdaf
           (if global adaptive forgetting factor is used)
        8. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__g2l_state_pdaf
            4. py__g2l_obs_pdaf (localise mean ensemble
               in observation space)
            5. py__init_obs_l_pdaf
            6. py__g2l_obs_pdaf
               (localise each ensemble member in observation space)
            7. py__init_obsvar_l_pdaf
               (only called if local adaptive forgetting
               factor `type_forget=2` is used)
            8. py__prodRinvA_l_pdaf
            9. core DA algorithm
            10. py__l2g_state_pdaf
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate`
       and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`

    References
    ----------
    .. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).
           A singular evolutive extended Kalman filter
           for data assimilation
           in oceanography. Journal of Marine systems,
           16(3-4), 323-340.

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
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assimilate_lseik(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__distribute_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_obs_pdaf,
                                 pdaf_cb.c__init_obs_l_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf,
                                 pdaf_cb.c__prodrinva_l_pdaf,
                                 pdaf_cb.c__init_n_domains_p_pdaf,
                                 pdaf_cb.c__init_dim_l_pdaf,
                                 pdaf_cb.c__init_dim_obs_l_pdaf,
                                 pdaf_cb.c__g2l_state_pdaf,
                                 pdaf_cb.c__l2g_state_pdaf,
                                 pdaf_cb.c__g2l_obs_pdaf,
                                 pdaf_cb.c__init_obsvar_pdaf,
                                 pdaf_cb.c__init_obsvar_l_pdaf,
                                 pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline_lseik(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_obs_pdaf, py__init_obs_l_pdaf, py__prepoststep_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__g2l_obs_pdaf, py__init_obsvar_pdaf, py__init_obsvar_l_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_obs_pdaf = <void*>py__init_obs_pdaf
    pdaf_cb.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    pdaf_cb.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    pdaf_cb.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_assim_offline_lseik(pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__init_obs_pdaf,
                                    pdaf_cb.c__init_obs_l_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__prodrinva_l_pdaf,
                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                    pdaf_cb.c__init_dim_l_pdaf,
                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                    pdaf_cb.c__g2l_state_pdaf,
                                    pdaf_cb.c__l2g_state_pdaf,
                                    pdaf_cb.c__g2l_obs_pdaf,
                                    pdaf_cb.c__init_obsvar_pdaf,
                                    pdaf_cb.c__init_obsvar_l_pdaf, &outflag)

    return outflag


def generate_obs(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf, py__init_obserr_f_pdaf,
    py__get_obs_f_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    """Generation of synthetic observations based on
    given error statistics and observation operator.

    When diagonal observation error covariance matrix is used,
    it is recommended to use
    :func:`pyPDAF.PDAF.omi_generate_obs` functionalities
    for fewer user-supplied functions and improved efficiency.

    The generated synthetic observations are based on
    each member of model forecast.
    Therefore, an ensemble of observations can be obtained.
    In a typical experiment,
    one may only need one ensemble member.
    The implementation strategy is similar to
    an assimilation step. This means that,
    one can reuse many user-supplied functions for
    assimilation and observation generation.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_generate_obs`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pda
        5. py__init_obserr_f_pdaf
        6. py__get_obs_f_pdaf
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

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
        Routine to provide time step, time and dimension of next observation

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
    pdaf_cb.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    pdaf_cb.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    pdaf_cb.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_generate_obs(pdaf_cb.c__collect_state_pdaf,
                             pdaf_cb.c__distribute_state_pdaf,
                             pdaf_cb.c__init_dim_obs_f_pdaf,
                             pdaf_cb.c__obs_op_f_pdaf,
                             pdaf_cb.c__init_obserr_f_pdaf,
                             pdaf_cb.c__get_obs_f_pdaf,
                             pdaf_cb.c__prepoststep_pdaf,
                             pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def generate_obs_offline(py__init_dim_obs_f_pdaf, py__obs_op_f_pdaf,
    py__init_obserr_f_pdaf, py__get_obs_f_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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

    py__get_obs_f_pdaf : Callable
        Provide observation vector

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    pdaf_cb.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    pdaf_cb.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf_generate_obs_offline(pdaf_cb.c__init_dim_obs_f_pdaf,
                                     pdaf_cb.c__obs_op_f_pdaf,
                                     pdaf_cb.c__init_obserr_f_pdaf,
                                     pdaf_cb.c__get_obs_f_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


