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

def init_obs_f_cb(int  step, int  dim_obs_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_f : int
        Dimension of full observation vector

    Returns
    -------
    observation_f : ndarray[np.float64, ndim=1]
        Full observation vector
        Array shape: (dim_obs_f)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] observation_f_np = np.zeros((dim_obs_f), dtype=np.float64, order="F")
    cdef double [::1] observation_f = observation_f_np
    with nogil:
        c__pdafomi_init_obs_f_cb(&step, &dim_obs_f, &observation_f[0])

    return observation_f_np


def init_obsvar_cb(int  step, int  dim_obs_p, double [::1] obs_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        PE-local dimension of observation vector
    obs_p : ndarray[np.float64, ndim=1]
        PE-local observation vector
        Array shape: (dim_obs_p)

    Returns
    -------
    meanvar : double
        Mean observation error variance
    """
    cdef double  meanvar
    with nogil:
        c__pdafomi_init_obsvar_cb(&step, &dim_obs_p, &obs_p[0], &meanvar)

    return meanvar


def init_obsvars_f_cb(int  step, int  dim_obs_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_f : int
        Dimension of full observation vector

    Returns
    -------
    var_f : ndarray[np.float64, ndim=1]
        vector of observation error variances
        Array shape: (dim_obs_f)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] var_f_np = np.zeros((dim_obs_f), dtype=np.float64, order="F")
    cdef double [::1] var_f = var_f_np
    with nogil:
        c__pdafomi_init_obsvars_f_cb(&step, &dim_obs_f, &var_f[0])

    return var_f_np


def g2l_obs_cb(int  domain_p, int  step, int  dim_obs_f, int  dim_obs_l,
    double [::1] ostate_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    ostate_f : ndarray[np.float64, ndim=1]
        Full PE-local obs.ervation vector
        Array shape: (dim_obs_f)

    Returns
    -------
    ostate_l : ndarray[np.float64, ndim=1]
        Observation vector on local domain
        Array shape: (dim_obs_l)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] ostate_l_np = np.zeros((dim_obs_l), dtype=np.float64, order="F")
    cdef double [::1] ostate_l = ostate_l_np
    with nogil:
        c__pdafomi_g2l_obs_cb(&domain_p, &step, &dim_obs_f, &dim_obs_l,
                              &ostate_f[0], &ostate_l[0])

    return ostate_l_np


def init_obs_l_cb(int  domain_p, int  step, int  dim_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain index
    step : int
        Current time step
    dim_obs_l : int
        Local dimension of observation vector

    Returns
    -------
    observation_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] observation_l_np = np.zeros((dim_obs_l), dtype=np.float64, order="F")
    cdef double [::1] observation_l = observation_l_np
    with nogil:
        c__pdafomi_init_obs_l_cb(&domain_p, &step, &dim_obs_l,
                                 &observation_l[0])

    return observation_l_np


def init_obsvar_l_cb(int  domain_p, int  step, int  dim_obs_l,
    double [::1] obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Index of current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        Local dimension of observation vector
    obs_l : ndarray[np.float64, ndim=1]
        Local observation vector
        Array shape: (dim_obs_l)

    Returns
    -------
    meanvar_l : double
        Mean local observation error variance
    """
    cdef double  meanvar_l
    with nogil:
        c__pdafomi_init_obsvar_l_cb(&domain_p, &step, &dim_obs_l,
                                    &obs_l[0], &meanvar_l)

    return meanvar_l


def prodrinva_l_cb(int  domain_p, int  step, int  dim_obs_l, int  rank,
    double [::1] obs_l, double [::1,:] a_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    obs_l : ndarray[np.float64, ndim=1]
        Local vector of observations
        Array shape: (dim_obs_l)
    a_l : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (dim_obs_l, rank)

    Returns
    -------
    a_l : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (dim_obs_l, rank)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix
        Array shape: (dim_obs_l, rank)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_l_np = np.asarray(a_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_l_np = np.zeros((dim_obs_l, rank), dtype=np.float64, order="F")
    cdef double [::1,:] c_l = c_l_np
    with nogil:
        c__pdafomi_prodrinva_l_cb(&domain_p, &step, &dim_obs_l, &rank,
                                  &obs_l[0], &a_l[0,0], &c_l[0,0])

    return a_l_np, c_l_np


def likelihood_l_cb(int  domain_p, int  step, int  dim_obs_l,
    double [::1] obs_l, double [::1] resid_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        PE-local dimension of obs. vector
    obs_l : ndarray[np.float64, ndim=1]
        PE-local vector of observations
        Array shape: (dim_obs_l)
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)

    Returns
    -------
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)
    lhood_l : double
        Output vector - log likelihood
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_l_np = np.asarray(resid_l, dtype=np.float64, order="F")
    cdef double  lhood_l
    with nogil:
        c__pdafomi_likelihood_l_cb(&domain_p, &step, &dim_obs_l, &obs_l[0],
                                   &resid_l[0], &lhood_l)

    return resid_l_np, lhood_l


def prodrinva_hyb_l_cb(int  domain_p, int  step, int  dim_obs_l, int  rank,
    double [::1] obs_l, double  alpha, double [::1,:] a_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    obs_l : ndarray[np.float64, ndim=1]
        Local vector of observations
        Array shape: (dim_obs_l)
    alpha : double
        Hybrid weight
    a_l : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (dim_obs_l, rank)

    Returns
    -------
    a_l : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (dim_obs_l, rank)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix
        Array shape: (dim_obs_l, rank)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_l_np = np.asarray(a_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_l_np = np.zeros((dim_obs_l, rank), dtype=np.float64, order="F")
    cdef double [::1,:] c_l = c_l_np
    with nogil:
        c__pdafomi_prodrinva_hyb_l_cb(&domain_p, &step, &dim_obs_l, &rank,
                                      &obs_l[0], &alpha, &a_l[0,0], &c_l[0,0])

    return a_l_np, c_l_np


def likelihood_hyb_l_cb(int  domain_p, int  step, int  dim_obs_l,
    double [::1] obs_l, double [::1] resid_l, double  alpha):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    step : int
        Current time step
    dim_obs_l : int
        PE-local dimension of obs. vector
    obs_l : ndarray[np.float64, ndim=1]
        PE-local vector of observations
        Array shape: (dim_obs_l)
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)
    alpha : double
        Hybrid weight

    Returns
    -------
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)
    lhood_l : double
        Output vector - log likelihood
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_l_np = np.asarray(resid_l, dtype=np.float64, order="F")
    cdef double  lhood_l
    with nogil:
        c__pdafomi_likelihood_hyb_l_cb(&domain_p, &step, &dim_obs_l,
                                       &obs_l[0], &resid_l[0], &alpha, &lhood_l)

    return resid_l_np, lhood_l


def prodrinva_cb(int  step, int  dim_obs_p, int  ncol, double [::1] obs_p,
    double [::1,:] a_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of PE-local observation vector
    ncol : int
        Number of columns in A_p and C_p
    obs_p : ndarray[np.float64, ndim=1]
        PE-local vector of observations
        Array shape: (dim_obs_p)
    a_p : ndarray[np.float64, ndim=2]
        Input matrix
        Array shape: (dim_obs_p, ncol)

    Returns
    -------
    c_p : ndarray[np.float64, ndim=2]
        Output matrix
        Array shape: (dim_obs_p, ncol)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_p_np = np.zeros((dim_obs_p, ncol), dtype=np.float64, order="F")
    cdef double [::1,:] c_p = c_p_np
    with nogil:
        c__pdafomi_prodrinva_cb(&step, &dim_obs_p, &ncol, &obs_p[0],
                                &a_p[0,0], &c_p[0,0])

    return c_p_np


def likelihood_cb(int  step, int  dim_obs, double [::1] obs,
    double [::1] resid):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs : int
        PE-local dimension of obs. vector
    obs : ndarray[np.float64, ndim=1]
        PE-local vector of observations
        Array shape: (dim_obs)
    resid : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs)

    Returns
    -------
    lhood : double
        Output vector - log likelihood
    """
    cdef double  lhood
    with nogil:
        c__pdafomi_likelihood_cb(&step, &dim_obs, &obs[0], &resid[0], &lhood)

    return lhood


def add_obs_error_cb(int  step, int  dim_obs_p, double [::1,:] c_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Dimension of PE-local observation vector
    c_p : ndarray[np.float64, ndim=2]
        Matrix to which R is added
        Array shape: (dim_obs_p,dim_obs_p)

    Returns
    -------
    c_p : ndarray[np.float64, ndim=2]
        Matrix to which R is added
        Array shape: (dim_obs_p,dim_obs_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_p_np = np.asarray(c_p, dtype=np.float64, order="F")
    with nogil:
        c__pdafomi_add_obs_error_cb(&step, &dim_obs_p, &c_p[0,0])

    return c_p_np


def init_obscovar_cb(int  step, int  dim_obs, int  dim_obs_p,
    double [::1] m_state_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs : int
        Dimension of observation vector
    dim_obs_p : int
        PE-local dimension of obs. vector
    m_state_p : ndarray[np.float64, ndim=1]
        Observation vector
        Array shape: (dim_obs_p)

    Returns
    -------
    covar : ndarray[np.float64, ndim=2]
        Observation error covar. matrix
        Array shape: (dim_obs,dim_obs)
    isdiag : bint
        Whether matrix R is diagonal
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] covar_np = np.zeros((dim_obs,dim_obs), dtype=np.float64, order="F")
    cdef double [::1,:] covar = covar_np
    cdef bint  isdiag
    with nogil:
        c__pdafomi_init_obscovar_cb(&step, &dim_obs, &dim_obs_p,
                                    &covar[0,0], &m_state_p[0], &isdiag)

    return covar_np, isdiag


def init_obserr_f_cb(int  step, int  dim_obs_f, double [::1] obs_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_f : int
        Full dimension of observation vector
    obs_f : ndarray[np.float64, ndim=1]
        Full observation vector
        Array shape: (dim_obs_f)

    Returns
    -------
    obserr_f : ndarray[np.float64, ndim=1]
        Full observation error stddev
        Array shape: (dim_obs_f)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obserr_f_np = np.zeros((dim_obs_f), dtype=np.float64, order="F")
    cdef double [::1] obserr_f = obserr_f_np
    with nogil:
        c__pdafomi_init_obserr_f_cb(&step, &dim_obs_f, &obs_f[0], &obserr_f[0])

    return obserr_f_np


def localize_covar_cb(int  dim_p, int  dim_obs, double [::1,:] hp_p,
    double [::1,:] hph):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_p : int
        Process-local state dimension
    dim_obs : int
        Number of observations
    hp_p : ndarray[np.float64, ndim=2]
        Process-local part of matrix HP
        Array shape: (dim_obs, dim_p)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH
        Array shape: (dim_obs, dim_obs)

    Returns
    -------
    hp_p : ndarray[np.float64, ndim=2]
        Process-local part of matrix HP
        Array shape: (dim_obs, dim_p)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH
        Array shape: (dim_obs, dim_obs)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hp_p_np = np.asarray(hp_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hph_np = np.asarray(hph, dtype=np.float64, order="F")
    with nogil:
        c__pdafomi_localize_covar_cb(&dim_p, &dim_obs, &hp_p[0,0], &hph[0,0])

    return hp_p_np, hph_np


def localize_covar_serial_cb(int  iobs, int  dim_p, int  dim_obs,
    double [::1] hp_p, double [::1] hxy_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
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

    Returns
    -------
    hp_p : ndarray[np.float64, ndim=1]
        Process-local part of matrix HP for observation iobs
        Array shape: (dim_p)
    hxy_p : ndarray[np.float64, ndim=1]
        Process-local part of matrix HX(HX_all) for full observations
        Array shape: (dim_obs)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hp_p_np = np.asarray(hp_p, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxy_p_np = np.asarray(hxy_p, dtype=np.float64, order="F")
    with nogil:
        c__pdafomi_localize_covar_serial_cb(&iobs, &dim_p, &dim_obs,
                                            &hp_p[0], &hxy_p[0])

    return hp_p_np, hxy_p_np


def omit_by_inno_l_cb(int  domain_p, int  dim_obs_l, double [::1] resid_l,
    double [::1] obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    domain_p : int
        Current local analysis domain
    dim_obs_l : int
        PE-local dimension of obs. vector
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Input vector of local observations
        Array shape: (dim_obs_l)

    Returns
    -------
    resid_l : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_l)
    obs_l : ndarray[np.float64, ndim=1]
        Input vector of local observations
        Array shape: (dim_obs_l)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_l_np = np.asarray(resid_l, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_l_np = np.asarray(obs_l, dtype=np.float64, order="F")
    with nogil:
        c__pdafomi_omit_by_inno_l_cb(&domain_p, &dim_obs_l, &resid_l[0],
                                     &obs_l[0])

    return resid_l_np, obs_l_np


def omit_by_inno_cb(int  dim_obs_f, double [::1] resid_f, double [::1] obs_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs_f : int
        Full dimension of obs. vector
    resid_f : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_f)
    obs_f : ndarray[np.float64, ndim=1]
        Input vector of full observations
        Array shape: (dim_obs_f)

    Returns
    -------
    resid_f : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (dim_obs_f)
    obs_f : ndarray[np.float64, ndim=1]
        Input vector of full observations
        Array shape: (dim_obs_f)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_f_np = np.asarray(resid_f, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_np = np.asarray(obs_f, dtype=np.float64, order="F")
    with nogil:
        c__pdafomi_omit_by_inno_cb(&dim_obs_f, &resid_f[0], &obs_f[0])

    return resid_f_np, obs_f_np


def g2l_cb(int  step, int  domain_p, int  dim_p, double [::1] state_p,
    int  dim_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    domain_p : int
        Current local analysis domain
    dim_p : int
        PE-local full state dimension
    state_p : ndarray[np.float64, ndim=1]
        PE-local full state vector
        Array shape: (dim_p)
    dim_l : int
        Local state dimension

    Returns
    -------
    state_l : ndarray[np.float64, ndim=1]
        State vector on local analysis domain
        Array shape: (dim_l)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.zeros((dim_l), dtype=np.float64, order="F")
    cdef double [::1] state_l = state_l_np
    with nogil:
        c__pdaflocal_g2l_cb(&step, &domain_p, &dim_p, &state_p[0], &dim_l,
                            &state_l[0])

    return state_l_np


def l2g_cb(int  step, int  domain_p, int  dim_l, double [::1] state_l,
    int  dim_p, double [::1] state_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    step : int
        Current time step
    domain_p : int
        Current local analysis domain
    dim_l : int
        Local state dimension
    state_l : ndarray[np.float64, ndim=1]
        State vector on local analysis domain
        Array shape: (dim_l)
    dim_p : int
        PE-local full state dimension
    state_p : ndarray[np.float64, ndim=1]
        PE-local full state vector
        Array shape: (dim_p)

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local full state vector
        Array shape: (dim_p)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaflocal_l2g_cb(&step, &domain_p, &dim_l, &state_l[0], &dim_p,
                            &state_p[0])

    return state_p_np


