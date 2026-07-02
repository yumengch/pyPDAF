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

def set_globalobs(int  globalobs_in):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    globalobs_in : int
        Input value of globalobs

    Returns
    -------
    """
    c__pdafomi_set_globalobs(&globalobs_in)


def g2l_obs(int  i_obs, double [::1] obs_f_all, double [::1] obs_l_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obs_f_all : ndarray[np.float64, ndim=1]
        Full obs. vector of current obs. for all variables
        Array shape: (:)
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables
        Array shape: (:)

    Returns
    -------
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_l_all_cfi
    cdef size_t obs_l_all_nbytes = obs_l_all.nbytes
    cdef CFI_index_t obs_l_all_extent[1]
    obs_l_all_extent[0] = obs_l_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_l_all_np = np.asarray(obs_l_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    CFI_establish(<CFI_cdesc_t *> &obs_f_all_cfi, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_l_all_cfi, &obs_l_all[0], CFI_attribute_other,
                    CFI_type_double , obs_l_all_nbytes, 1, obs_l_all_extent)

    c__pdafomi_g2l_obs(&i_obs, <CFI_cdesc_t *> &obs_f_all_cfi, <CFI_cdesc_t *> &obs_l_all_cfi)

    return obs_l_all_np


def init_obs_l(int  i_obs, double [::1] obs_l_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables
        Array shape: (:)

    Returns
    -------
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_l_all_cfi
    cdef size_t obs_l_all_nbytes = obs_l_all.nbytes
    cdef CFI_index_t obs_l_all_extent[1]
    obs_l_all_extent[0] = obs_l_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_l_all_np = np.asarray(obs_l_all, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &obs_l_all_cfi, &obs_l_all[0], CFI_attribute_other,
                      CFI_type_double , obs_l_all_nbytes, 1, obs_l_all_extent)

    c__pdafomi_init_obs_l(&i_obs, <CFI_cdesc_t *> &obs_l_all_cfi)

    return obs_l_all_np


def init_obsvar_l(int  i_obs, double  meanvar_l, int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    meanvar_l : double
        Mean variance
    cnt_obs_l : int
        Observation counter

    Returns
    -------
    meanvar_l : double
        Mean variance
    cnt_obs_l : int
        Observation counter
    """
    c__pdafomi_init_obsvar_l(&i_obs, &meanvar_l, &cnt_obs_l)

    return meanvar_l, cnt_obs_l


def prodrinva_l(int  i_obs, int  nobs_all, int  ncols, double [::1,:] a_l,
    double [::1,:] c_l, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nobs_all : int
        Dimension of local obs. vector (all obs. types)
    ncols : int
        Rank of initial covariance matrix
    a_l : ndarray[np.float64, ndim=2]
        Input matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    verbose : int
        Verbosity flag

    Returns
    -------
    a_l : ndarray[np.float64, ndim=2]
        Input matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 a_l_cfi
    cdef size_t a_l_nbytes = a_l.nbytes
    cdef CFI_index_t a_l_extent[2]
    a_l_extent[0] = a_l.shape[0]
    a_l_extent[1] = a_l.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_l_np = np.asarray(a_l, dtype=np.float64, order="F")

    cdef CFI_cdesc_rank2 c_l_cfi
    cdef size_t c_l_nbytes = c_l.nbytes
    cdef CFI_index_t c_l_extent[2]
    c_l_extent[0] = c_l.shape[0]
    c_l_extent[1] = c_l.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_l_np = np.asarray(c_l, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &a_l_cfi, &a_l[0,0], CFI_attribute_other,
                      CFI_type_double , a_l_nbytes, 2, a_l_extent)

    CFI_establish(<CFI_cdesc_t *> &c_l_cfi, &c_l[0,0], CFI_attribute_other,
                    CFI_type_double , c_l_nbytes, 2, c_l_extent)

    c__pdafomi_prodrinva_l(&i_obs, &nobs_all, &ncols, <CFI_cdesc_t *> &a_l_cfi, <CFI_cdesc_t *> &c_l_cfi,
                            &verbose)

    return a_l_np, c_l_np


def prodrinva_hyb_l(int  i_obs, int  nobs_all, int  ncols, double  gamma,
    double [::1,:] a_l, double [::1,:] c_l, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nobs_all : int
        Dimension of local obs. vector (all obs. types)
    ncols : int
        Rank of initial covariance matrix
    gamma : double
        Hybrid weight
    a_l : ndarray[np.float64, ndim=2]
        Input matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    verbose : int
        Verbosity flag

    Returns
    -------
    a_l : ndarray[np.float64, ndim=2]
        Input matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    c_l : ndarray[np.float64, ndim=2]
        Output matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 a_l_cfi
    cdef size_t a_l_nbytes = a_l.nbytes
    cdef CFI_index_t a_l_extent[2]
    a_l_extent[0] = a_l.shape[0]
    a_l_extent[1] = a_l.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] a_l_np = np.asarray(a_l, dtype=np.float64, order="F")

    cdef CFI_cdesc_rank2 c_l_cfi
    cdef size_t c_l_nbytes = c_l.nbytes
    cdef CFI_index_t c_l_extent[2]
    c_l_extent[0] = c_l.shape[0]
    c_l_extent[1] = c_l.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_l_np = np.asarray(c_l, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &a_l_cfi, &a_l[0,0], CFI_attribute_other,
                      CFI_type_double , a_l_nbytes, 2, a_l_extent)

    CFI_establish(<CFI_cdesc_t *> &c_l_cfi, &c_l[0,0], CFI_attribute_other,
                    CFI_type_double , c_l_nbytes, 2, c_l_extent)

    c__pdafomi_prodrinva_hyb_l(&i_obs, &nobs_all, &ncols, &gamma,
                                <CFI_cdesc_t *> &a_l_cfi, <CFI_cdesc_t *> &c_l_cfi, &verbose)

    return a_l_np, c_l_np


def likelihood_l(int  i_obs, double [::1] resid_l_all, double  lhood_l,
    int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    resid_l_all : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (:)
    lhood_l : double
        Output vector - log likelihood
    verbose : int
        Verbosity flag

    Returns
    -------
    resid_l_all : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (:)
    lhood_l : double
        Output vector - log likelihood
    """
    cdef CFI_cdesc_rank1 resid_l_all_cfi
    cdef size_t resid_l_all_nbytes = resid_l_all.nbytes
    cdef CFI_index_t resid_l_all_extent[1]
    resid_l_all_extent[0] = resid_l_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_l_all_np = np.asarray(resid_l_all, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &resid_l_all_cfi, &resid_l_all[0], CFI_attribute_other,
                      CFI_type_double , resid_l_all_nbytes, 1, resid_l_all_extent)

    c__pdafomi_likelihood_l(&i_obs, <CFI_cdesc_t *> &resid_l_all_cfi, &lhood_l, &verbose)

    return resid_l_all_np, lhood_l


def likelihood_hyb_l(int  i_obs, double [::1] resid_l_all, double  gamma,
    double  lhood_l, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    resid_l_all : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (:)
    gamma : double
        Hybrid weight
    lhood_l : double
        Output vector - log likelihood
    verbose : int
        Verbosity flag

    Returns
    -------
    resid_l_all : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (:)
    lhood_l : double
        Output vector - log likelihood
    """
    cdef CFI_cdesc_rank1 resid_l_all_cfi
    cdef size_t resid_l_all_nbytes = resid_l_all.nbytes
    cdef CFI_index_t resid_l_all_extent[1]
    resid_l_all_extent[0] = resid_l_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] resid_l_all_np = np.asarray(resid_l_all, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &resid_l_all_cfi, &resid_l_all[0], CFI_attribute_other,
                      CFI_type_double , resid_l_all_nbytes, 1, resid_l_all_extent)

    c__pdafomi_likelihood_hyb_l(&i_obs, <CFI_cdesc_t *> &resid_l_all_cfi, &gamma,
                                &lhood_l, &verbose)

    return resid_l_all_np, lhood_l


def g2l_obs_internal(int  i_obs, double [::1] obs_f_one,
    int  offset_obs_l_all, double [::1] obs_l_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obs_f_one : ndarray[np.float64, ndim=1]
        Full obs. vector of current obs. type (nobs_f_one)
        Array shape: (:)
    offset_obs_l_all : int
        Offset of current observation in obs_l_all and ivar_l_all
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables (nobs_l_all)
        Array shape: (:)

    Returns
    -------
    obs_l_all : ndarray[np.float64, ndim=1]
        Local observation vector for all variables (nobs_l_all)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_l_all_cfi
    cdef size_t obs_l_all_nbytes = obs_l_all.nbytes
    cdef CFI_index_t obs_l_all_extent[1]
    obs_l_all_extent[0] = obs_l_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_l_all_np = np.asarray(obs_l_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_one_cfi
    cdef size_t obs_f_one_nbytes = obs_f_one.nbytes
    cdef CFI_index_t obs_f_one_extent[1]
    obs_f_one_extent[0] = obs_f_one.shape[0]
    CFI_establish(<CFI_cdesc_t *> &obs_f_one_cfi, &obs_f_one[0], CFI_attribute_other,
                      CFI_type_double , obs_f_one_nbytes, 1, obs_f_one_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_l_all_cfi, &obs_l_all[0], CFI_attribute_other,
                    CFI_type_double , obs_l_all_nbytes, 1, obs_l_all_extent)

    c__pdafomi_g2l_obs_internal(&i_obs, <CFI_cdesc_t *> &obs_f_one_cfi,
                                &offset_obs_l_all, <CFI_cdesc_t *> &obs_l_all_cfi)

    return obs_l_all_np


def comp_dist2(int  i_obs, double [::1] coordsa, double [::1] coordsb,
    int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain (ncoord)
        Array shape: (:)
    coordsb : ndarray[np.float64, ndim=1]
        Coordinates of observation (ncoord)
        Array shape: (:)
    verbose : int
        Control screen output

    Returns
    -------
    distance2 : double
        Squared distance
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    cdef CFI_cdesc_rank1 coordsb_cfi
    cdef size_t coordsb_nbytes = coordsb.nbytes
    cdef CFI_index_t coordsb_extent[1]
    coordsb_extent[0] = coordsb.shape[0]
    cdef double  distance2
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    CFI_establish(<CFI_cdesc_t *> &coordsb_cfi, &coordsb[0], CFI_attribute_other,
                    CFI_type_double , coordsb_nbytes, 1, coordsb_extent)

    c__pdafomi_comp_dist2(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, <CFI_cdesc_t *> &coordsb_cfi, &distance2,
                            &verbose)

    return distance2


def check_dist2(int  i_obs, double [::1] coordsa, double [::1] coordsb,
    int  verbose, int  cnt_obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain (ncoord)
        Array shape: (:)
    coordsb : ndarray[np.float64, ndim=1]
        Coordinates of observation (ncoord)
        Array shape: (:)
    verbose : int
        Control screen output
    cnt_obs : int
        Count number of local observations

    Returns
    -------
    distance2 : double
        Squared distance
    checkdist : bint
        Flag whether distance is within cut-off radius
    cnt_obs : int
        Count number of local observations
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    cdef CFI_cdesc_rank1 coordsb_cfi
    cdef size_t coordsb_nbytes = coordsb.nbytes
    cdef CFI_index_t coordsb_extent[1]
    coordsb_extent[0] = coordsb.shape[0]
    cdef double  distance2
    cdef bint  checkdist
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    CFI_establish(<CFI_cdesc_t *> &coordsb_cfi, &coordsb[0], CFI_attribute_other,
                    CFI_type_double , coordsb_nbytes, 1, coordsb_extent)

    c__pdafomi_check_dist2(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, <CFI_cdesc_t *> &coordsb_cfi,
                            &distance2, &checkdist, &verbose, &cnt_obs)

    return distance2, checkdist, cnt_obs


def check_dist2_noniso(int  i_obs, double [::1] coordsa,
    double [::1] coordsb, double [::1] dists, double  sradius,
    int  verbose, int  cnt_obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain (ncoord)
        Array shape: (:)
    coordsb : ndarray[np.float64, ndim=1]
        Coordinates of observation (ncoord)
        Array shape: (:)
    dists : ndarray[np.float64, ndim=1]
        Vector of distance in each coordinate direction
        Array shape: (:)
    sradius : double
        Directional support radius
    verbose : int
        Control screen output
    cnt_obs : int
        Count number of local observations

    Returns
    -------
    distance2 : double
        Squared distance
    dists : ndarray[np.float64, ndim=1]
        Vector of distance in each coordinate direction
        Array shape: (:)
    cradius : double
        Directional cut-off radius
    sradius : double
        Directional support radius
    checkdist : bint
        Flag whether distance is within cut-off radius
    cnt_obs : int
        Count number of local observations
    """
    cdef CFI_cdesc_rank1 dists_cfi
    cdef size_t dists_nbytes = dists.nbytes
    cdef CFI_index_t dists_extent[1]
    dists_extent[0] = dists.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] dists_np = np.asarray(dists, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    cdef CFI_cdesc_rank1 coordsb_cfi
    cdef size_t coordsb_nbytes = coordsb.nbytes
    cdef CFI_index_t coordsb_extent[1]
    coordsb_extent[0] = coordsb.shape[0]
    cdef double  distance2
    cdef double  cradius
    cdef bint  checkdist
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    CFI_establish(<CFI_cdesc_t *> &coordsb_cfi, &coordsb[0], CFI_attribute_other,
                    CFI_type_double , coordsb_nbytes, 1, coordsb_extent)

    CFI_establish(<CFI_cdesc_t *> &dists_cfi, &dists[0], CFI_attribute_other,
                    CFI_type_double , dists_nbytes, 1, dists_extent)

    c__pdafomi_check_dist2_noniso(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, <CFI_cdesc_t *> &coordsb_cfi,
                                    &distance2, <CFI_cdesc_t *> &dists_cfi, &cradius,
                                    &sradius, &checkdist, &verbose, &cnt_obs)

    return distance2, dists_np, cradius, sradius, checkdist, cnt_obs


def weights_l(int  verbose, int  nobs_l, int  ncols, int  locweight,
    double [::1] cradius, double [::1] sradius, double [::1,:] mata,
    double [::1] ivar_obs_l, double [::1] dist_l, double [::1] weight_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    verbose : int
        Verbosity flag
    nobs_l : int
        Number of local observations
    ncols : int

    locweight : int
        Localization weight type
    cradius : ndarray[np.float64, ndim=1]
        Localization cut-off radius
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        support radius for weight functions
        Array shape: (:)
    mata : ndarray[np.float64, ndim=2]

        Array shape: (:,:)
    ivar_obs_l : ndarray[np.float64, ndim=1]
        Local vector of inverse obs. variances (nobs_l)
        Array shape: (:)
    dist_l : ndarray[np.float64, ndim=1]
        Local vector of obs. distances (nobs_l)
        Array shape: (:)
    weight_l : ndarray[np.float64, ndim=1]
        Output: vector of weights
        Array shape: (:)

    Returns
    -------
    weight_l : ndarray[np.float64, ndim=1]
        Output: vector of weights
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 weight_l_cfi
    cdef size_t weight_l_nbytes = weight_l.nbytes
    cdef CFI_index_t weight_l_extent[1]
    weight_l_extent[0] = weight_l.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] weight_l_np = np.asarray(weight_l, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    cdef CFI_cdesc_rank2 mata_cfi
    cdef size_t mata_nbytes = mata.nbytes
    cdef CFI_index_t mata_extent[2]
    mata_extent[0] = mata.shape[0]
    mata_extent[1] = mata.shape[1]
    cdef CFI_cdesc_rank1 ivar_obs_l_cfi
    cdef size_t ivar_obs_l_nbytes = ivar_obs_l.nbytes
    cdef CFI_index_t ivar_obs_l_extent[1]
    ivar_obs_l_extent[0] = ivar_obs_l.shape[0]
    cdef CFI_cdesc_rank1 dist_l_cfi
    cdef size_t dist_l_nbytes = dist_l.nbytes
    cdef CFI_index_t dist_l_extent[1]
    dist_l_extent[0] = dist_l.shape[0]
    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    CFI_establish(<CFI_cdesc_t *> &mata_cfi, &mata[0,0], CFI_attribute_other,
                    CFI_type_double , mata_nbytes, 2, mata_extent)

    CFI_establish(<CFI_cdesc_t *> &ivar_obs_l_cfi, &ivar_obs_l[0], CFI_attribute_other,
                    CFI_type_double , ivar_obs_l_nbytes, 1, ivar_obs_l_extent)

    CFI_establish(<CFI_cdesc_t *> &dist_l_cfi, &dist_l[0], CFI_attribute_other,
                    CFI_type_double , dist_l_nbytes, 1, dist_l_extent)

    CFI_establish(<CFI_cdesc_t *> &weight_l_cfi, &weight_l[0], CFI_attribute_other,
                    CFI_type_double , weight_l_nbytes, 1, weight_l_extent)

    c__pdafomi_weights_l(&verbose, &nobs_l, &ncols, &locweight,
                            <CFI_cdesc_t *> &cradius_cfi, <CFI_cdesc_t *> &sradius_cfi,
                            <CFI_cdesc_t *> &mata_cfi,
                            <CFI_cdesc_t *> &ivar_obs_l_cfi,
                            <CFI_cdesc_t *> &dist_l_cfi,
                            <CFI_cdesc_t *> &weight_l_cfi)

    return weight_l_np


def weights_l_sgnl(int  verbose, int  nobs_l, int  ncols, int  locweight,
    double  cradius, double  sradius, double [::1,:] mata,
    double [::1] ivar_obs_l, double [::1] dist_l, double [::1] weight_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    verbose : int
        Verbosity flag
    nobs_l : int
        Number of local observations
    ncols : int

    locweight : int
        Localization weight type
    cradius : double
        Localization cut-off radius
    sradius : double
        support radius for weight functions
    mata : ndarray[np.float64, ndim=2]

        Array shape: (:,:)
    ivar_obs_l : ndarray[np.float64, ndim=1]
        Local vector of inverse obs. variances (nobs_l)
        Array shape: (:)
    dist_l : ndarray[np.float64, ndim=1]
        Local vector of obs. distances (nobs_l)
        Array shape: (:)
    weight_l : ndarray[np.float64, ndim=1]
        Output: vector of weights
        Array shape: (:)

    Returns
    -------
    weight_l : ndarray[np.float64, ndim=1]
        Output: vector of weights
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 weight_l_cfi
    cdef size_t weight_l_nbytes = weight_l.nbytes
    cdef CFI_index_t weight_l_extent[1]
    weight_l_extent[0] = weight_l.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] weight_l_np = np.asarray(weight_l, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 mata_cfi
    cdef size_t mata_nbytes = mata.nbytes
    cdef CFI_index_t mata_extent[2]
    mata_extent[0] = mata.shape[0]
    mata_extent[1] = mata.shape[1]
    cdef CFI_cdesc_rank1 ivar_obs_l_cfi
    cdef size_t ivar_obs_l_nbytes = ivar_obs_l.nbytes
    cdef CFI_index_t ivar_obs_l_extent[1]
    ivar_obs_l_extent[0] = ivar_obs_l.shape[0]
    cdef CFI_cdesc_rank1 dist_l_cfi
    cdef size_t dist_l_nbytes = dist_l.nbytes
    cdef CFI_index_t dist_l_extent[1]
    dist_l_extent[0] = dist_l.shape[0]
    CFI_establish(<CFI_cdesc_t *> &mata_cfi, &mata[0,0], CFI_attribute_other,
                      CFI_type_double , mata_nbytes, 2, mata_extent)

    CFI_establish(<CFI_cdesc_t *> &ivar_obs_l_cfi, &ivar_obs_l[0], CFI_attribute_other,
                    CFI_type_double , ivar_obs_l_nbytes, 1, ivar_obs_l_extent)

    CFI_establish(<CFI_cdesc_t *> &dist_l_cfi, &dist_l[0], CFI_attribute_other,
                    CFI_type_double , dist_l_nbytes, 1, dist_l_extent)

    CFI_establish(<CFI_cdesc_t *> &weight_l_cfi, &weight_l[0], CFI_attribute_other,
                    CFI_type_double , weight_l_nbytes, 1, weight_l_extent)

    c__pdafomi_weights_l_sgnl(&verbose, &nobs_l, &ncols, &locweight,
                                &cradius, &sradius, <CFI_cdesc_t *> &mata_cfi,
                                <CFI_cdesc_t *> &ivar_obs_l_cfi,
                                <CFI_cdesc_t *> &dist_l_cfi, <CFI_cdesc_t *> &weight_l_cfi)

    return weight_l_np


def omit_by_inno_l(int  i_obs, double [::1] inno_l,
    double [::1] obs_l_all, int  obsid, int  cnt_all, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    inno_l : ndarray[np.float64, ndim=1]
        Input vector of observation innovation
        Array shape: (:)
    obs_l_all : ndarray[np.float64, ndim=1]
        Input vector of local observations
        Array shape: (:)
    obsid : int
        ID of observation type
    cnt_all : int
        Count of omitted observation over all types
    verbose : int
        Verbosity flag

    Returns
    -------
    cnt_all : int
        Count of omitted observation over all types
    """
    cdef CFI_cdesc_rank1 inno_l_cfi
    cdef size_t inno_l_nbytes = inno_l.nbytes
    cdef CFI_index_t inno_l_extent[1]
    inno_l_extent[0] = inno_l.shape[0]
    cdef CFI_cdesc_rank1 obs_l_all_cfi
    cdef size_t obs_l_all_nbytes = obs_l_all.nbytes
    cdef CFI_index_t obs_l_all_extent[1]
    obs_l_all_extent[0] = obs_l_all.shape[0]
    CFI_establish(<CFI_cdesc_t *> &inno_l_cfi, &inno_l[0], CFI_attribute_other,
                      CFI_type_double , inno_l_nbytes, 1, inno_l_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_l_all_cfi, &obs_l_all[0], CFI_attribute_other,
                    CFI_type_double , obs_l_all_nbytes, 1, obs_l_all_extent)

    c__pdafomi_omit_by_inno_l(&i_obs, <CFI_cdesc_t *> &inno_l_cfi, <CFI_cdesc_t *> &obs_l_all_cfi,
                                &obsid, &cnt_all, &verbose)

    return cnt_all


def obsstats_l(int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    """
    c__pdafomi_obsstats_l(&screen)



def dealloc():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    c__pdafomi_dealloc()



def ocoord_all(int  ncoord, double [::1,:] oc_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    ncoord : int
        Number of coordinate directions
    oc_all : ndarray[np.float64, ndim=2]
        Array of observation coordinates size(ncoord, dim_obs)
        Array shape: (:,:)

    Returns
    -------
    oc_all : ndarray[np.float64, ndim=2]
        Array of observation coordinates size(ncoord, dim_obs)
        Array shape: (:,:)
    """
    cdef CFI_cdesc_rank2 oc_all_cfi
    cdef size_t oc_all_nbytes = oc_all.nbytes
    cdef CFI_index_t oc_all_extent[2]
    oc_all_extent[0] = oc_all.shape[0]
    oc_all_extent[1] = oc_all.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] oc_all_np = np.asarray(oc_all, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &oc_all_cfi, &oc_all[0,0], CFI_attribute_other,
                      CFI_type_double , oc_all_nbytes, 2, oc_all_extent)

    c__pdafomi_ocoord_all(&ncoord, <CFI_cdesc_t *> &oc_all_cfi)

    return oc_all_np


def local_weight(int  wtype, int  rtype, double  cradius, double  sradius,
    double  distance, int  nrows, int  ncols, double [::1,:] a,
    double  var_obs, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    c__pdafomi_local_weight(&wtype, &rtype, &cradius, &sradius,
                                &distance, &nrows, &ncols, &a[0,0],
                                &var_obs, &weight, &verbose)

    return weight


def check_dist2_loop(int  i_obs, double [::1] coordsa, int  cnt_obs,
    int  mode):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain (ncoord)
        Array shape: (:)
    cnt_obs : int
        Count number of local observations
    mode : int
        1: count local observations

    Returns
    -------
    cnt_obs : int
        Count number of local observations
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_loop(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_loop_opt(int  i_obs, double [::1] coordsa, int  cnt_obs,
    int  mode):
    """check_dist2_loop_opt(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Optimized variant of :func:`check_dist2_loop`.

    Parameters
    ----------
    i_obs : int
        Index into observation arrays.
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of the current analysis domain.
    cnt_obs : int
        Count of local observations.
    mode : int
        Local-observation processing mode.

    Returns
    -------
    cnt_obs : int
        Updated count of local observations.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_loop_opt(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_loop_sort(int  i_obs, double [::1] coordsa, int  cnt_obs,
    int  mode):
    """check_dist2_loop_sort(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Sorted-observation variant of :func:`check_dist2_loop`.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_loop_sort(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_loop_sort2(int  i_obs, double [::1] coordsa, int  cnt_obs,
    int  mode):
    """check_dist2_loop_sort2(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Upper-loop sorted-observation variant of :func:`check_dist2_loop`.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_loop_sort2(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def set_thisobs_l(int  i_obs, double [::1] coordsx, int  n_obs):
    """set_thisobs_l(i_obs: int, coordsx: np.ndarray, n_obs: int) -> int

    Store local observation metadata for one OMI observation type.

    Parameters
    ----------
    i_obs : int
        Index into observation arrays.
    coordsx : ndarray[np.float64, ndim=1]
        Coordinates of the current local analysis domain.
    n_obs : int
        Number of local observations.

    Returns
    -------
    n_obs : int
        Updated number of local observations.
    """
    cdef CFI_cdesc_rank1 coordsx_cfi
    cdef size_t coordsx_nbytes = coordsx.nbytes
    cdef CFI_index_t coordsx_extent[1]
    coordsx_extent[0] = coordsx.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsx_cfi, &coordsx[0], CFI_attribute_other,
                      CFI_type_double , coordsx_nbytes, 1, coordsx_extent)

    c__pdafomi_set_thisobs_l(&i_obs, <CFI_cdesc_t *> &coordsx_cfi, &n_obs)

    return n_obs


def check_dist2_noniso_loop(int  i_obs, double [::1] coordsa,
    int  cnt_obs, int  mode):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coordsa : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain (ncoord)
        Array shape: (:)
    cnt_obs : int
        Count number of local observations
    mode : int
        1: count local observations

    Returns
    -------
    cnt_obs : int
        Count number of local observations
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_noniso_loop(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_noniso_loop_opt(int  i_obs, double [::1] coordsa,
    int  cnt_obs, int  mode):
    """check_dist2_noniso_loop_opt(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Optimized non-isotropic variant of :func:`check_dist2_noniso_loop`.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_noniso_loop_opt(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_noniso_loop_sort(int  i_obs, double [::1] coordsa,
    int  cnt_obs, int  mode):
    """check_dist2_noniso_loop_sort(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Sorted-observation non-isotropic variant of
    :func:`check_dist2_noniso_loop`.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_noniso_loop_sort(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_dist2_noniso_loop_sort2(int  i_obs, double [::1] coordsa,
    int  cnt_obs, int  mode):
    """check_dist2_noniso_loop_sort2(i_obs: int, coordsa: np.ndarray, cnt_obs: int, mode: int) -> int

    Upper-loop sorted-observation non-isotropic variant of
    :func:`check_dist2_noniso_loop`.
    """
    cdef CFI_cdesc_rank1 coordsa_cfi
    cdef size_t coordsa_nbytes = coordsa.nbytes
    cdef CFI_index_t coordsa_extent[1]
    coordsa_extent[0] = coordsa.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coordsa_cfi, &coordsa[0], CFI_attribute_other,
                      CFI_type_double , coordsa_nbytes, 1, coordsa_extent)

    c__pdafomi_check_dist2_noniso_loop_sort2(&i_obs, <CFI_cdesc_t *> &coordsa_cfi, &cnt_obs, &mode)

    return cnt_obs


def check_noniso(int  i_obs, int  idx, int  cnt_obs,
    double [::1] dists, double  distance2, int  mode):
    """check_noniso(i_obs: int, idx: int, cnt_obs: int, dists: np.ndarray, distance2: float, mode: int) -> Tuple[int, np.ndarray]

    Check one observation against non-isotropic OMI localization.
    """
    cdef CFI_cdesc_rank1 dists_cfi
    cdef size_t dists_nbytes = dists.nbytes
    cdef CFI_index_t dists_extent[1]
    dists_extent[0] = dists.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] dists_np = np.asarray(dists, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &dists_cfi, &dists[0], CFI_attribute_other,
                      CFI_type_double , dists_nbytes, 1, dists_extent)

    c__pdafomi_check_noniso(&i_obs, &idx, &cnt_obs, <CFI_cdesc_t *> &dists_cfi,
                                &distance2, &mode)

    return cnt_obs, dists_np


def tree_idx_lower(int  row, double  tst, double [::1,:] points,
    int  npts, int  ilower, double  offset, int  level):
    """tree_idx_lower(row: int, tst: float, points: np.ndarray, npts: int, ilower: int, offset: float, level: int) -> Tuple[int, float, int]

    Return the lower index bound for sorted OMI coordinate searches.
    """
    cdef CFI_cdesc_rank2 points_cfi
    cdef size_t points_nbytes = points.nbytes
    cdef CFI_index_t points_extent[2]
    points_extent[0] = points.shape[0]
    points_extent[1] = points.shape[1]
    CFI_establish(<CFI_cdesc_t *> &points_cfi, &points[0,0], CFI_attribute_other,
                      CFI_type_double , points_nbytes, 2, points_extent)

    c__pdafomi_tree_idx_lower(&row, &tst, <CFI_cdesc_t *> &points_cfi,
                                  &npts, &ilower, &offset, &level)

    return ilower, offset, level


def tree_idx_upper(int  row, double  tst, double [::1,:] points,
    int  npts, int  iupper, double  offset, int  level):
    """tree_idx_upper(row: int, tst: float, points: np.ndarray, npts: int, iupper: int, offset: float, level: int) -> Tuple[int, float, int]

    Return the upper index bound for sorted OMI coordinate searches.
    """
    cdef CFI_cdesc_rank2 points_cfi
    cdef size_t points_nbytes = points.nbytes
    cdef CFI_index_t points_extent[2]
    points_extent[0] = points.shape[0]
    points_extent[1] = points.shape[1]
    CFI_establish(<CFI_cdesc_t *> &points_cfi, &points[0,0], CFI_attribute_other,
                      CFI_type_double , points_nbytes, 2, points_extent)

    c__pdafomi_tree_idx_upper(&row, &tst, <CFI_cdesc_t *> &points_cfi,
                                  &npts, &iupper, &offset, &level)

    return iupper, offset, level


def obs_op_gatheronly(int  i_obs, double [::1] state_p,
    double [::1] obs_f_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)

    Returns
    -------
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    CFI_establish(<CFI_cdesc_t *> &state_p_cfi, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_f_all_cfi, &obs_f_all[0], CFI_attribute_other,
                    CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

    c__pdafomi_obs_op_gatheronly(&i_obs, <CFI_cdesc_t *> &state_p_cfi,
                                    <CFI_cdesc_t *> &obs_f_all_cfi)

    return obs_f_all_np


def obs_op_adj_gatheronly(int  i_obs, double [::1] obs_f_all,
    double [::1] state_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)

    Returns
    -------
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    CFI_establish(<CFI_cdesc_t *> &obs_f_all_cfi, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

    CFI_establish(<CFI_cdesc_t *> &state_p_cfi, &state_p[0], CFI_attribute_other,
                    CFI_type_double , state_p_nbytes, 1, state_p_extent)

    c__pdafomi_obs_op_adj_gatheronly(&i_obs, <CFI_cdesc_t *> &obs_f_all_cfi,
                                        <CFI_cdesc_t *> &state_p_cfi)

    return obs_f_all_np


def init_obs_f(int  i_obs, int  dim_obs_f, double [::1] obsstate_f,
    int  offset):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim_obs_f : int
        Dimension of full observed state (all observed fields)
    obsstate_f : ndarray[np.float64, ndim=1]
        Full observation vector (dim_obs_f)
        Array shape: (:)
    offset : int
        input: offset of module-type observations in obsstate_f

    Returns
    -------
    obsstate_f : ndarray[np.float64, ndim=1]
        Full observation vector (dim_obs_f)
        Array shape: (:)
    offset : int
        input: offset of module-type observations in obsstate_f
    """
    cdef CFI_cdesc_rank1 obsstate_f_cfi
    cdef size_t obsstate_f_nbytes = obsstate_f.nbytes
    cdef CFI_index_t obsstate_f_extent[1]
    obsstate_f_extent[0] = obsstate_f.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obsstate_f_np = np.asarray(obsstate_f, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &obsstate_f_cfi, &obsstate_f[0], CFI_attribute_other,
                      CFI_type_double , obsstate_f_nbytes, 1, obsstate_f_extent)

    c__pdafomi_init_obs_f(&i_obs, &dim_obs_f, <CFI_cdesc_t *> &obsstate_f_cfi, &offset)

    return obsstate_f_np, offset


def init_obsvars_f(int  i_obs, int  dim_obs_f, double [::1] var_f,
    int  offset):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim_obs_f : int
        Dimension of full observed state (all observed fields)
    var_f : ndarray[np.float64, ndim=1]
        Full vector of observation variances (dim_obs_f)
        Array shape: (:)
    offset : int
        input: offset of module-type observations in obsstate_f

    Returns
    -------
    var_f : ndarray[np.float64, ndim=1]
        Full vector of observation variances (dim_obs_f)
        Array shape: (:)
    offset : int
        input: offset of module-type observations in obsstate_f
    """
    cdef CFI_cdesc_rank1 var_f_cfi
    cdef size_t var_f_nbytes = var_f.nbytes
    cdef CFI_index_t var_f_extent[1]
    var_f_extent[0] = var_f.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] var_f_np = np.asarray(var_f, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &var_f_cfi, &var_f[0], CFI_attribute_other,
                      CFI_type_double , var_f_nbytes, 1, var_f_extent)

    c__pdafomi_init_obsvars_f(&i_obs, &dim_obs_f, <CFI_cdesc_t *> &var_f_cfi, &offset)

    return var_f_np, offset


def init_obsvar_f(int  i_obs, double  meanvar, int  cnt_obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    meanvar : double
        Mean variance
    cnt_obs : int
        Observation counter

    Returns
    -------
    meanvar : double
        Mean variance
    cnt_obs : int
        Observation counter
    """
    c__pdafomi_init_obsvar_f(&i_obs, &meanvar, &cnt_obs)

    return meanvar, cnt_obs


def prodrinva(int  i_obs, int  ncols, double [::1,:] a_p, double [::1,:] c_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    ncols : int
        Number of columns in A_p and C_p
    a_p : ndarray[np.float64, ndim=2]
        Input matrix (nobs_f, ncols)
        Array shape: (:, :)
    c_p : ndarray[np.float64, ndim=2]
        Output matrix (nobs_f, ncols)
        Array shape: (:, :)

    Returns
    -------
    c_p : ndarray[np.float64, ndim=2]
        Output matrix (nobs_f, ncols)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 c_p_cfi
    cdef size_t c_p_nbytes = c_p.nbytes
    cdef CFI_index_t c_p_extent[2]
    c_p_extent[0] = c_p.shape[0]
    c_p_extent[1] = c_p.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] c_p_np = np.asarray(c_p, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 a_p_cfi
    cdef size_t a_p_nbytes = a_p.nbytes
    cdef CFI_index_t a_p_extent[2]
    a_p_extent[0] = a_p.shape[0]
    a_p_extent[1] = a_p.shape[1]
    CFI_establish(<CFI_cdesc_t *> &a_p_cfi, &a_p[0,0], CFI_attribute_other,
                      CFI_type_double , a_p_nbytes, 2, a_p_extent)

    CFI_establish(<CFI_cdesc_t *> &c_p_cfi, &c_p[0,0], CFI_attribute_other,
                    CFI_type_double , c_p_nbytes, 2, c_p_extent)

    c__pdafomi_prodrinva(&i_obs, &ncols, <CFI_cdesc_t *> &a_p_cfi, <CFI_cdesc_t *> &c_p_cfi)

    return c_p_np


def likelihood(int  i_obs, double [::1] resid, double  lhood):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    resid : ndarray[np.float64, ndim=1]
        Input vector of residuum
        Array shape: (:)
    lhood : double
        Output vector - log likelihood

    Returns
    -------
    lhood : double
        Output vector - log likelihood
    """
    cdef CFI_cdesc_rank1 resid_cfi
    cdef size_t resid_nbytes = resid.nbytes
    cdef CFI_index_t resid_extent[1]
    resid_extent[0] = resid.shape[0]
    CFI_establish(<CFI_cdesc_t *> &resid_cfi, &resid[0], CFI_attribute_other,
                      CFI_type_double , resid_nbytes, 1, resid_extent)

    c__pdafomi_likelihood(&i_obs, <CFI_cdesc_t *> &resid_cfi, &lhood)

    return lhood


def add_obs_error(int  i_obs, int  nobs_all, double [::1,:] matc):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nobs_all : int
        Number of observations
    matc : ndarray[np.float64, ndim=2]
        Input/Output matrix (nobs_f, rank)
        Array shape: (:, :)

    Returns
    -------
    matc : ndarray[np.float64, ndim=2]
        Input/Output matrix (nobs_f, rank)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 matc_cfi
    cdef size_t matc_nbytes = matc.nbytes
    cdef CFI_index_t matc_extent[2]
    matc_extent[0] = matc.shape[0]
    matc_extent[1] = matc.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] matc_np = np.asarray(matc, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &matc_cfi, &matc[0,0], CFI_attribute_other,
                      CFI_type_double , matc_nbytes, 2, matc_extent)

    c__pdafomi_add_obs_error(&i_obs, &nobs_all, <CFI_cdesc_t *> &matc_cfi)

    return matc_np


def init_obscovar(int  i_obs, int  nobs_all, double [::1,:] covar):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nobs_all : int
        Number of observations
    covar : ndarray[np.float64, ndim=2]
        Input/Output matrix (nobs_all, nobs_all)
        Array shape: (:, :)

    Returns
    -------
    covar : ndarray[np.float64, ndim=2]
        Input/Output matrix (nobs_all, nobs_all)
        Array shape: (:, :)
    isdiag : bint
        Whether matrix R is diagonal
    """
    cdef CFI_cdesc_rank2 covar_cfi
    cdef size_t covar_nbytes = covar.nbytes
    cdef CFI_index_t covar_extent[2]
    covar_extent[0] = covar.shape[0]
    covar_extent[1] = covar.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] covar_np = np.asarray(covar, dtype=np.float64, order="F")
    cdef bint  isdiag
    CFI_establish(<CFI_cdesc_t *> &covar_cfi, &covar[0,0], CFI_attribute_other,
                      CFI_type_double , covar_nbytes, 2, covar_extent)

    c__pdafomi_init_obscovar(&i_obs, &nobs_all, <CFI_cdesc_t *> &covar_cfi, &isdiag)

    return covar_np, isdiag


def init_obserr_f(int  i_obs, double [::1] obserr_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obserr_f : ndarray[np.float64, ndim=1]
        Full vector of observation errors
        Array shape: (:)

    Returns
    -------
    obserr_f : ndarray[np.float64, ndim=1]
        Full vector of observation errors
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obserr_f_cfi
    cdef size_t obserr_f_nbytes = obserr_f.nbytes
    cdef CFI_index_t obserr_f_extent[1]
    obserr_f_extent[0] = obserr_f.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obserr_f_np = np.asarray(obserr_f, dtype=np.float64, order="F")
    CFI_establish(<CFI_cdesc_t *> &obserr_f_cfi, &obserr_f[0], CFI_attribute_other,
                      CFI_type_double , obserr_f_nbytes, 1, obserr_f_extent)

    c__pdafomi_init_obserr_f(&i_obs, <CFI_cdesc_t *> &obserr_f_cfi)

    return obserr_f_np


def get_local_ids_obs_f(int  dim_obs_g, double  lradius,
    double [::1,:] oc_f, int [::1] id_lim, int  disttype,
    double [::1] domainsize):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs_g : int
        Global full number of observations
    lradius : double
        Localization radius (used is a constant one here)
    oc_f : ndarray[np.float64, ndim=2]
        observation coordinates (radians), row 1: lon, 2: lat
        Array shape: (:,:)
    id_lim : ndarray[np.intc, ndim=1]
        Indices of process-local full obs. in global full vector
        Array shape: (:)
    disttype : int
        type of distance computation
    domainsize : ndarray[np.float64, ndim=1]
        Global size of model domain
        Array shape: (:)

    Returns
    -------
    cnt_lim : int
        Number of full observation for local process domain
    id_lim : ndarray[np.intc, ndim=1]
        Indices of process-local full obs. in global full vector
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 id_lim_cfi
    cdef size_t id_lim_nbytes = id_lim.nbytes
    cdef CFI_index_t id_lim_extent[1]
    id_lim_extent[0] = id_lim.shape[0]
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] id_lim_np = np.asarray(id_lim, dtype=np.intc, order="F")
    cdef CFI_cdesc_rank2 oc_f_cfi
    cdef size_t oc_f_nbytes = oc_f.nbytes
    cdef CFI_index_t oc_f_extent[2]
    oc_f_extent[0] = oc_f.shape[0]
    oc_f_extent[1] = oc_f.shape[1]
    cdef CFI_cdesc_rank1 domainsize_cfi
    cdef size_t domainsize_nbytes = domainsize.nbytes
    cdef CFI_index_t domainsize_extent[1]
    domainsize_extent[0] = domainsize.shape[0]
    cdef int  cnt_lim
    CFI_establish(<CFI_cdesc_t *> &oc_f_cfi, &oc_f[0,0], CFI_attribute_other,
                      CFI_type_double , oc_f_nbytes, 2, oc_f_extent)

    CFI_establish(<CFI_cdesc_t *> &id_lim_cfi, &id_lim[0], CFI_attribute_other,
                    CFI_type_int , id_lim_nbytes, 1, id_lim_extent)

    CFI_establish(<CFI_cdesc_t *> &domainsize_cfi, &domainsize[0], CFI_attribute_other,
                    CFI_type_double , domainsize_nbytes, 1, domainsize_extent)

    c__pdafomi_get_local_ids_obs_f(&dim_obs_g, &lradius, <CFI_cdesc_t *> &oc_f_cfi,
                                    &cnt_lim, <CFI_cdesc_t *> &id_lim_cfi, &disttype,
                                    <CFI_cdesc_t *> &domainsize_cfi)

    return cnt_lim, id_lim_np


def limit_obs_f(int  i_obs, int  offset, double [::1] obs_f_one,
    double [::1] obs_f_lim):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    offset : int
        offset of this observation in obs_f_lim
    obs_f_one : ndarray[np.float64, ndim=1]
        Global full observation vector (nobs_f)
        Array shape: (:)
    obs_f_lim : ndarray[np.float64, ndim=1]
        full observation vector for process domains (nobs_lim)
        Array shape: (:)

    Returns
    -------
    obs_f_lim : ndarray[np.float64, ndim=1]
        full observation vector for process domains (nobs_lim)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_f_lim_cfi
    cdef size_t obs_f_lim_nbytes = obs_f_lim.nbytes
    cdef CFI_index_t obs_f_lim_extent[1]
    obs_f_lim_extent[0] = obs_f_lim.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_lim_np = np.asarray(obs_f_lim, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_one_cfi
    cdef size_t obs_f_one_nbytes = obs_f_one.nbytes
    cdef CFI_index_t obs_f_one_extent[1]
    obs_f_one_extent[0] = obs_f_one.shape[0]
    CFI_establish(<CFI_cdesc_t *> &obs_f_one_cfi, &obs_f_one[0], CFI_attribute_other,
                      CFI_type_double , obs_f_one_nbytes, 1, obs_f_one_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_f_lim_cfi, &obs_f_lim[0], CFI_attribute_other,
                    CFI_type_double , obs_f_lim_nbytes, 1, obs_f_lim_extent)

    c__pdafomi_limit_obs_f(&i_obs, &offset, <CFI_cdesc_t *> &obs_f_one_cfi, <CFI_cdesc_t *> &obs_f_lim_cfi)

    return obs_f_lim_np


def gather_dim_obs_f(int  dim_obs_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

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
    c__pdafomi_gather_dim_obs_f(&dim_obs_p, &dim_obs_f)

    return dim_obs_f


def gather_obs_f_flex(int  dim_obs_p, double [::1] obs_p, double [::1] obs_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension
    obs_p : ndarray[np.float64, ndim=1]
        PE-local vector
        Array shape: (:)
    obs_f : ndarray[np.float64, ndim=1]
        Full gathered vector
        Array shape: (:)

    Returns
    -------
    obs_f : ndarray[np.float64, ndim=1]
        Full gathered vector
        Array shape: (:)
    status : int
        Status flag: (0) no error
    """
    cdef CFI_cdesc_rank1 obs_f_cfi
    cdef size_t obs_f_nbytes = obs_f.nbytes
    cdef CFI_index_t obs_f_extent[1]
    obs_f_extent[0] = obs_f.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_np = np.asarray(obs_f, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_p_cfi
    cdef size_t obs_p_nbytes = obs_p.nbytes
    cdef CFI_index_t obs_p_extent[1]
    obs_p_extent[0] = obs_p.shape[0]
    cdef int  status
    CFI_establish(<CFI_cdesc_t *> &obs_p_cfi, &obs_p[0], CFI_attribute_other,
                      CFI_type_double , obs_p_nbytes, 1, obs_p_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_f_cfi, &obs_f[0], CFI_attribute_other,
                    CFI_type_double , obs_f_nbytes, 1, obs_f_extent)

    c__pdafomi_gather_obs_f_flex(&dim_obs_p, <CFI_cdesc_t *> &obs_p_cfi, <CFI_cdesc_t *> &obs_f_cfi, &status)

    return obs_f_np, status


def gather_obs_f2_flex(int  dim_obs_p, double [::1,:] coords_p,
    double [::1,:] coords_f, int  nrows):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    dim_obs_p : int
        PE-local observation dimension
    coords_p : ndarray[np.float64, ndim=2]
        PE-local array
        Array shape: (:,:)
    coords_f : ndarray[np.float64, ndim=2]
        Full gathered array
        Array shape: (:,:)
    nrows : int
        Number of rows in array

    Returns
    -------
    coords_f : ndarray[np.float64, ndim=2]
        Full gathered array
        Array shape: (:,:)
    status : int
        Status flag: (0) no error
    """
    cdef CFI_cdesc_rank2 coords_f_cfi
    cdef size_t coords_f_nbytes = coords_f.nbytes
    cdef CFI_index_t coords_f_extent[2]
    coords_f_extent[0] = coords_f.shape[0]
    coords_f_extent[1] = coords_f.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] coords_f_np = np.asarray(coords_f, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 coords_p_cfi
    cdef size_t coords_p_nbytes = coords_p.nbytes
    cdef CFI_index_t coords_p_extent[2]
    coords_p_extent[0] = coords_p.shape[0]
    coords_p_extent[1] = coords_p.shape[1]
    cdef int  status
    CFI_establish(<CFI_cdesc_t *> &coords_p_cfi, &coords_p[0,0], CFI_attribute_other,
                      CFI_type_double , coords_p_nbytes, 2, coords_p_extent)

    CFI_establish(<CFI_cdesc_t *> &coords_f_cfi, &coords_f[0,0], CFI_attribute_other,
                    CFI_type_double , coords_f_nbytes, 2, coords_f_extent)

    c__pdafomi_gather_obs_f2_flex(&dim_obs_p, <CFI_cdesc_t *> &coords_p_cfi,
                                    <CFI_cdesc_t *> &coords_f_cfi, &nrows, &status)
    return coords_f_np, status


def omit_by_inno(int  i_obs, double [::1] inno_f, double [::1] obs_f_all,
    int  obsid, int  cnt_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    inno_f : ndarray[np.float64, ndim=1]
        Input vector of observation innovation
        Array shape: (:)
    obs_f_all : ndarray[np.float64, ndim=1]
        Input vector of local observations
        Array shape: (:)
    obsid : int
        ID of observation type
    cnt_all : int
        Count of omitted observation over all types

    Returns
    -------
    cnt_all : int
        Count of omitted observation over all types
    """
    cdef CFI_cdesc_rank1 inno_f_cfi
    cdef size_t inno_f_nbytes = inno_f.nbytes
    cdef CFI_index_t inno_f_extent[1]
    inno_f_extent[0] = inno_f.shape[0]
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    CFI_establish(<CFI_cdesc_t *> &inno_f_cfi, &inno_f[0], CFI_attribute_other,
                      CFI_type_double , inno_f_nbytes, 1, inno_f_extent)

    CFI_establish(<CFI_cdesc_t *> &obs_f_all_cfi, &obs_f_all[0], CFI_attribute_other,
                    CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

    c__pdafomi_omit_by_inno(&i_obs, <CFI_cdesc_t *> &inno_f_cfi,
                            <CFI_cdesc_t *> &obs_f_all_cfi, &obsid,
                            &cnt_all)

    return cnt_all


def obsstats(int  screen):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    screen : int
        Verbosity flag

    Returns
    -------
    """
    c__pdafomi_obsstats(&screen)



def gather_obsdims():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    c__pdafomi_gather_obsdims()


def dealloc_local():
    r"""dealloc_local() -> None

    Deallocate PDAF-OMI local-observation storage.

    This releases the local ``obs_l`` array allocated by
    :func:`init_local`. It is useful when local-domain observation state should
    be reset between local analysis phases.

    Returns
    -------
    None
    """
    c__pdafomi_dealloc_local()
