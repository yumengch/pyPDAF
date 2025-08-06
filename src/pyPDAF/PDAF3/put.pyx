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

def put_state_3dvar_nondiagr(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_3dvar_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__prodrinva_pdaf,
                                          pdaf_cb.c__cvt_pdaf,
                                          pdaf_cb.c__cvt_adj_pdaf,
                                          pdaf_cb.c__obs_op_lin_pdaf,
                                          pdaf_cb.c__obs_op_adj_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_en3dvar_estkf_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prodrinva_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_en3dvar_estkf_nondiagr(
                                                  pdaf_cb.c__collect_state_pdaf,
                                                  pdaf_cb.c__init_dim_obs_pdaf,
                                                  pdaf_cb.c__obs_op_pdaf,
                                                  pdaf_cb.c__prodrinva_pdaf,
                                                  pdaf_cb.c__cvt_ens_pdaf,
                                                  pdaf_cb.c__cvt_adj_ens_pdaf,
                                                  pdaf_cb.c__obs_op_lin_pdaf,
                                                  pdaf_cb.c__obs_op_adj_pdaf,
                                                  pdaf_cb.c__prepoststep_pdaf,
                                                  &outflag)

    return outflag


def put_state_en3dvar_lestkf_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prodrinva_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_en3dvar_lestkf_nondiagr(
                                                   pdaf_cb.c__collect_state_pdaf,
                                                   pdaf_cb.c__init_dim_obs_pdaf,
                                                   pdaf_cb.c__obs_op_pdaf,
                                                   pdaf_cb.c__prodrinva_pdaf,
                                                   pdaf_cb.c__cvt_ens_pdaf,
                                                   pdaf_cb.c__cvt_adj_ens_pdaf,
                                                   pdaf_cb.c__obs_op_lin_pdaf,
                                                   pdaf_cb.c__obs_op_adj_pdaf,
                                                   pdaf_cb.c__prodrinva_l_pdaf,
                                                   pdaf_cb.c__init_n_domains_p_pdaf,
                                                   pdaf_cb.c__init_dim_l_pdaf,
                                                   pdaf_cb.c__init_dim_obs_l_pdaf,
                                                   pdaf_cb.c__prepoststep_pdaf,
                                                   &outflag)

    return outflag


def put_state_hyb3dvar_estkf_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prodrinva_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_hyb3dvar_estkf_nondiagr(
                                                   pdaf_cb.c__collect_state_pdaf,
                                                   pdaf_cb.c__init_dim_obs_pdaf,
                                                   pdaf_cb.c__obs_op_pdaf,
                                                   pdaf_cb.c__prodrinva_pdaf,
                                                   pdaf_cb.c__cvt_ens_pdaf,
                                                   pdaf_cb.c__cvt_adj_ens_pdaf,
                                                   pdaf_cb.c__cvt_pdaf,
                                                   pdaf_cb.c__cvt_adj_pdaf,
                                                   pdaf_cb.c__obs_op_lin_pdaf,
                                                   pdaf_cb.c__obs_op_adj_pdaf,
                                                   pdaf_cb.c__prepoststep_pdaf,
                                                   &outflag)

    return outflag


def put_state_hyb3dvar_lestkf_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prodrinva_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_hyb3dvar_lestkf_nondiagr(
                                                    pdaf_cb.c__collect_state_pdaf,
                                                    pdaf_cb.c__init_dim_obs_pdaf,
                                                    pdaf_cb.c__obs_op_pdaf,
                                                    pdaf_cb.c__prodrinva_pdaf,
                                                    pdaf_cb.c__cvt_ens_pdaf,
                                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                                    pdaf_cb.c__cvt_pdaf,
                                                    pdaf_cb.c__cvt_adj_pdaf,
                                                    pdaf_cb.c__obs_op_lin_pdaf,
                                                    pdaf_cb.c__obs_op_adj_pdaf,
                                                    pdaf_cb.c__prodrinva_l_pdaf,
                                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                                    pdaf_cb.c__init_dim_l_pdaf,
                                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                                    pdaf_cb.c__prepoststep_pdaf,
                                                    &outflag)

    return outflag


def put_state_3dvar_all(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_3dvar_all(pdaf_cb.c__collect_state_pdaf,
                                     pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                     pdaf_cb.c__init_dim_l_pdaf,
                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_3dvar(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_put_state_3dvar`
    or :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`.

    PDAF-OMI modules require fewer user-supplied
    functions and improved efficiency.

    3DVar DA for a single DA step.

    Compared to :func:`pyPDAF.PDAF.assimilate_3dvar`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will not
    be assigned by user-supplied functions as well.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable transformation.
    This is a deterministic filtering scheme so no ensemble
    and parallelisation is needed.
    This function should be called at each model time step.

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

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_put_state_3dvar`
       and :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_3dvar(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__cvt_pdaf,
                                 pdaf_cb.c__cvt_adj_pdaf,
                                 pdaf_cb.c__obs_op_lin_pdaf,
                                 pdaf_cb.c__obs_op_adj_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_en3dvar(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_en3dvar(pdaf_cb.c__collect_state_pdaf,
                                   pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__cvt_ens_pdaf,
                                   pdaf_cb.c__cvt_adj_ens_pdaf,
                                   pdaf_cb.c__obs_op_lin_pdaf,
                                   pdaf_cb.c__obs_op_adj_pdaf,
                                   pdaf_cb.c__init_n_domains_p_pdaf,
                                   pdaf_cb.c__init_dim_l_pdaf,
                                   pdaf_cb.c__init_dim_obs_l_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_en3dvar_estkf(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`
    or :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step.

    Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_estkf`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will not be
    assigned by user-supplied functions as well.
    This function is typically used when there are
    not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    The background error covariance matrix is
    estimated by an ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An ESTKF is used along with 3DEnVar to
    generate ensemble perturbations.
    This function should be called at each model time step.

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
               (only relevant for adaptive forgetting factor schemes)
            6. py__prodRinvA_pdaf
            7. core ESTKF algorithm

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`
       and :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_en3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                         pdaf_cb.c__init_dim_obs_pdaf,
                                         pdaf_cb.c__obs_op_pdaf,
                                         pdaf_cb.c__cvt_ens_pdaf,
                                         pdaf_cb.c__cvt_adj_ens_pdaf,
                                         pdaf_cb.c__obs_op_lin_pdaf,
                                         pdaf_cb.c__obs_op_adj_pdaf,
                                         pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_en3dvar_lestkf(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`
    or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step without post-processing,
    distributing analysis, and setting next observation step,
    where the ensemble anomaly is generated by LESTKF.

    Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_lestkf`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will
    not be assigned by user-supplied functions as well.
    This function is typically used when there are
    not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    The background error covariance matrix is estimated by ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An LESTKF is used to generate ensemble perturbations.
    This function should be called at each model time step.

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
                   (localise each ensemble member in observation space)
                7. py__init_obsvar_l_pdaf
                   (only called if local adaptive forgetting factor
                   `type_forget=2` is used)
                8. py__prodRinvA_l_pdaf
                9. core DA algorithm
                10. py__l2g_state_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`
       and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_en3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__cvt_ens_pdaf,
                                          pdaf_cb.c__cvt_adj_ens_pdaf,
                                          pdaf_cb.c__obs_op_lin_pdaf,
                                          pdaf_cb.c__obs_op_adj_pdaf,
                                          pdaf_cb.c__init_n_domains_p_pdaf,
                                          pdaf_cb.c__init_dim_l_pdaf,
                                          pdaf_cb.c__init_dim_obs_l_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_hyb3dvar(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_hyb3dvar(pdaf_cb.c__collect_state_pdaf,
                                    pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__cvt_ens_pdaf,
                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                    pdaf_cb.c__cvt_pdaf,
                                    pdaf_cb.c__cvt_adj_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf,
                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                    pdaf_cb.c__init_dim_l_pdaf,
                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_hyb3dvar_estkf(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf,
    py__cvt_adj_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__prepoststep_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`
    or :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.

    Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_estkf`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will
    not be assigned by user-supplied functions as well.
    This function is typically used when there are
    not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    ESTKF in this implementation.
    This function should be called at each model time step.

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

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`
       and :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_hyb3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__cvt_ens_pdaf,
                                          pdaf_cb.c__cvt_adj_ens_pdaf,
                                          pdaf_cb.c__cvt_pdaf,
                                          pdaf_cb.c__cvt_adj_pdaf,
                                          pdaf_cb.c__obs_op_lin_pdaf,
                                          pdaf_cb.c__obs_op_adj_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_hyb3dvar_lestkf(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`
    or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step using
    non-diagnoal observation error covariance matrix
    without post-processing, distributing analysis,
    and setting next observation step, where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.

    Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_lestkf`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will
    not be assigned by user-supplied functions as well.
    This function is typically used when there are
    not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    LESTKF in this implementation.
    This function should be called at each model time step.

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

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`
       and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_hyb3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__cvt_ens_pdaf,
                                           pdaf_cb.c__cvt_adj_ens_pdaf,
                                           pdaf_cb.c__cvt_pdaf,
                                           pdaf_cb.c__cvt_adj_pdaf,
                                           pdaf_cb.c__obs_op_lin_pdaf,
                                           pdaf_cb.c__obs_op_adj_pdaf,
                                           pdaf_cb.c__init_n_domains_p_pdaf,
                                           pdaf_cb.c__init_dim_l_pdaf,
                                           pdaf_cb.c__init_dim_obs_l_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           &outflag)

    return outflag


def put_state_local_nondiagr(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prepoststep_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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

    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_local_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__init_n_domains_p_pdaf,
                                          pdaf_cb.c__init_dim_l_pdaf,
                                          pdaf_cb.c__init_dim_obs_l_pdaf,
                                          pdaf_cb.c__prodrinva_l_pdaf, &outflag)

    return outflag


def put_state_global_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prodrinva_pdaf,
    py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__prodrinva_pdaf : Callable
        Provide product of inverse of R with matrix A

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_global_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__prodrinva_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           &outflag)

    return outflag


def put_state_lnetf_nondiagr(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prepoststep_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__likelihood_l_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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

    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_lnetf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__init_n_domains_p_pdaf,
                                          pdaf_cb.c__init_dim_l_pdaf,
                                          pdaf_cb.c__init_dim_obs_l_pdaf,
                                          pdaf_cb.c__likelihood_l_pdaf,
                                          &outflag)

    return outflag


def put_state_lknetf_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prepoststep_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf,
    py__prodrinva_hyb_l_pdaf, py__likelihood_l_pdaf, py__likelihood_hyb_l_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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

    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A

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
        Product R^-1 A on local analysis domain with hybrid weight

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

    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization

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
        Compute likelihood and apply localization with tempering

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_lknetf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__init_n_domains_p_pdaf,
                                           pdaf_cb.c__init_dim_l_pdaf,
                                           pdaf_cb.c__init_dim_obs_l_pdaf,
                                           pdaf_cb.c__prodrinva_l_pdaf,
                                           pdaf_cb.c__prodrinva_hyb_l_pdaf,
                                           pdaf_cb.c__likelihood_l_pdaf,
                                           pdaf_cb.c__likelihood_hyb_l_pdaf,
                                           &outflag)

    return outflag


def put_state_enkf_nondiagr(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__add_obs_err_pdaf, py__init_obs_covar_pdaf,
    py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__init_obs_covar_pdaf : Callable
        Initialize mean observation error variance

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_enkf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                         pdaf_cb.c__init_dim_obs_pdaf,
                                         pdaf_cb.c__obs_op_pdaf,
                                         pdaf_cb.c__add_obs_err_pdaf,
                                         pdaf_cb.c__init_obs_covar_pdaf,
                                         pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_lenkf_nondiagr(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prepoststep_pdaf, py__localize_covar_pdaf,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Apply covariance localization

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
        Initialize mean observation error variance

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_lenkf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__localize_covar_pdaf,
                                          pdaf_cb.c__add_obs_err_pdaf,
                                          pdaf_cb.c__init_obs_covar_pdaf,
                                          &outflag)

    return outflag


def put_state_nonlin_nondiagr(py__collect_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__likelihood_pdaf,
    py__prepoststep_pdaf):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__likelihood_pdaf : Callable
        Compute likelihood

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


    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_nonlin_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__likelihood_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           &outflag)

    return outflag


def put_state(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf, int  outflag):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_put_state(pdaf_cb.c__collect_state_pdaf,
                           pdaf_cb.c__init_dim_obs_pdaf,
                           pdaf_cb.c__obs_op_pdaf,
                           pdaf_cb.c__init_n_domains_p_pdaf,
                           pdaf_cb.c__init_dim_l_pdaf,
                           pdaf_cb.c__init_dim_obs_l_pdaf,
                           pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_local(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__prepoststep_pdaf, int  outflag):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
        Initialize local dimimension of obs. vector

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
        Get local state from full state

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
        Init full state from local state

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_put_state_local(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__init_n_domains_p_pdaf,
                                 pdaf_cb.c__init_dim_l_pdaf,
                                 pdaf_cb.c__init_dim_obs_l_pdaf,
                                 pdaf_cb.c__g2l_state_pdaf,
                                 pdaf_cb.c__l2g_state_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_global(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prepoststep_pdaf, int  outflag):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_put_state_global(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_lenkf(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__localize_covar_pdaf, py__prepoststep_pdaf,
    int  outflag):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_put_state_lenkf`
    or :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Stochastic EnKF (ensemble Kalman filter)
    with covariance localisation [1]_
    for a single DA step without OMI.

    Compared to :func:`pyPDAF.PDAF.assimilate_lenkf`,
    this function has no :func:`get_state` call.
    This means that the analysis is not post-processed,
    and distributed to the model forecast
    by user-supplied functions. The next DA step will
    not be assigned by user-supplied functions as well.
    This function is typically used when there are
    not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    This is the only scheme for covariance localisation in PDAF.

    This function should be called at each model time step.

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

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_put_state_lenkf`
       and :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman
           Filter Technique.
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_put_state_lenkf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__localize_covar_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_ensrf(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__localize_covar_serial_pdaf, py__prepoststep_pdaf,
    int  outflag):
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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__localize_covar_serial_pdaf : Callable
        Apply localization to HP and BXY for single observation

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

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_put_state_ensrf(pdaf_cb.c__collect_state_pdaf,
                                 pdaf_cb.c__init_dim_obs_pdaf,
                                 pdaf_cb.c__obs_op_pdaf,
                                 pdaf_cb.c__localize_covar_serial_pdaf,
                                 pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def put_state_generate_obs(py__collect_state_pdaf, py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__get_obs_f_pdaf, py__prepoststep_pdaf):
    """Generation of synthetic observations
    based on given error statistics and observation operator
    without post-processing, distributing analysis,
    and setting next observation step.

    When diagonal observation error covariance matrix is used,
    it is recommended to use
    :func:`pyPDAF.PDAF.omi_generate_obs` functionalities
    for fewer user-supplied functions and improved efficiency.

    The generated synthetic observations are
    based on each member of model forecast.
    Therefore, an ensemble of observations can be obtained.
    In a typical experiment,
    one may only need one ensemble member.

    Compared to :func:`pyPDAF.PDAF.generate_obs`,
    this function has no :func:`get_state` call.
    This means that the next DA step will
    not be assigned by user-supplied functions.
    This function is typically used when there
    are not enough CPUs to run the ensemble in parallel,
    and some ensemble members have to be run serially.
    The :func:`pyPDAF.PDAF.get_state` function follows this
    function call to ensure the sequential DA.

    The implementation strategy is similar to
    an assimilation step. This means that,
    one can reuse many user-supplied functions for
    assimilation and observation generation.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pda
        5. py__init_obserr_f_pdaf
        6. py__get_obs_f_pdaf

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

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

        Callback Parameters
        -------------------
        step : int
                current time step

        Callback Returns
        ----------------
        dim_obs_p : int
                dimension of observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

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

    py__get_obs_f_pdaf : Callable
        Initialize observation vector

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
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_put_state_generate_obs(pdaf_cb.c__collect_state_pdaf,
                                        pdaf_cb.c__init_dim_obs_pdaf,
                                        pdaf_cb.c__obs_op_pdaf,
                                        pdaf_cb.c__get_obs_f_pdaf,
                                        pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


