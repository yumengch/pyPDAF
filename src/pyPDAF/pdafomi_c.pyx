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

def init(int  n_obs):
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
    n_obs : int
        number of observations

    Returns
    -------
    """
    with nogil:
        c__pdafomi_init(&n_obs)



def init_local():
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.
    """
    with nogil:
        c__pdafomi_init_local()



def check_error(int  flag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    flag : int
        Error flag

    Returns
    -------
    flag : int
        Error flag
    """
    with nogil:
        c__pdafomi_check_error(&flag)

    return flag


def gather_obs(int  i_obs, int  dim_obs_p, double [::1] obs_p,
    double [::1] ivar_obs_p, double [::1,:] ocoord_p, int  ncoord,
    double  lradius):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim_obs_p : int
        Number of process-local observation
    obs_p : ndarray[np.float64, ndim=1]
        Vector of process-local observations
        Array shape: (:)
    ivar_obs_p : ndarray[np.float64, ndim=1]
        Vector of process-local inverse observation error variance
        Array shape: (:)
    ocoord_p : ndarray[np.float64, ndim=2]
        Array of process-local observation coordinates
        Array shape: (:,:)
    ncoord : int
        Number of rows of coordinate array
    lradius : double
        Localization radius (the maximum radius used in this process domain)

    Returns
    -------
    dim_obs_f : int
        Full number of observations
    """
    cdef CFI_cdesc_rank1 obs_p_cfi
    cdef CFI_cdesc_t *obs_p_ptr = <CFI_cdesc_t *> &obs_p_cfi
    cdef size_t obs_p_nbytes = obs_p.nbytes
    cdef CFI_index_t obs_p_extent[1]
    obs_p_extent[0] = obs_p.shape[0]
    cdef CFI_cdesc_rank1 ivar_obs_p_cfi
    cdef CFI_cdesc_t *ivar_obs_p_ptr = <CFI_cdesc_t *> &ivar_obs_p_cfi
    cdef size_t ivar_obs_p_nbytes = ivar_obs_p.nbytes
    cdef CFI_index_t ivar_obs_p_extent[1]
    ivar_obs_p_extent[0] = ivar_obs_p.shape[0]
    cdef CFI_cdesc_rank2 ocoord_p_cfi
    cdef CFI_cdesc_t *ocoord_p_ptr = <CFI_cdesc_t *> &ocoord_p_cfi
    cdef size_t ocoord_p_nbytes = ocoord_p.nbytes
    cdef CFI_index_t ocoord_p_extent[2]
    ocoord_p_extent[0] = ocoord_p.shape[0]
    ocoord_p_extent[1] = ocoord_p.shape[1]
    cdef int  dim_obs_f
    with nogil:
        CFI_establish(obs_p_ptr, &obs_p[0], CFI_attribute_other,
                      CFI_type_double , obs_p_nbytes, 1, obs_p_extent)

        CFI_establish(ivar_obs_p_ptr, &ivar_obs_p[0], CFI_attribute_other,
                      CFI_type_double , ivar_obs_p_nbytes, 1, ivar_obs_p_extent)

        CFI_establish(ocoord_p_ptr, &ocoord_p[0,0], CFI_attribute_other,
                      CFI_type_double , ocoord_p_nbytes, 2, ocoord_p_extent)

        c__pdafomi_gather_obs(&i_obs, &dim_obs_p, obs_p_ptr,
                              ivar_obs_p_ptr, ocoord_p_ptr, &ncoord,
                              &lradius, &dim_obs_f)

    return dim_obs_f


def gather_obsstate(int  i_obs, double [::1] obsstate_p,
    double [::1] obsstate_f):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    obsstate_p : ndarray[np.float64, ndim=1]
        Vector of process-local observed state
        Array shape: (:)
    obsstate_f : ndarray[np.float64, ndim=1]
        Full observed vector for all types
        Array shape: (:)

    Returns
    -------
    obsstate_f : ndarray[np.float64, ndim=1]
        Full observed vector for all types
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obsstate_f_cfi
    cdef CFI_cdesc_t *obsstate_f_ptr = <CFI_cdesc_t *> &obsstate_f_cfi
    cdef size_t obsstate_f_nbytes = obsstate_f.nbytes
    cdef CFI_index_t obsstate_f_extent[1]
    obsstate_f_extent[0] = obsstate_f.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obsstate_f_np = np.asarray(obsstate_f, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obsstate_p_cfi
    cdef CFI_cdesc_t *obsstate_p_ptr = <CFI_cdesc_t *> &obsstate_p_cfi
    cdef size_t obsstate_p_nbytes = obsstate_p.nbytes
    cdef CFI_index_t obsstate_p_extent[1]
    obsstate_p_extent[0] = obsstate_p.shape[0]
    with nogil:
        CFI_establish(obsstate_p_ptr, &obsstate_p[0], CFI_attribute_other,
                      CFI_type_double , obsstate_p_nbytes, 1, obsstate_p_extent)

        CFI_establish(obsstate_f_ptr, &obsstate_f[0], CFI_attribute_other,
                      CFI_type_double , obsstate_f_nbytes, 1, obsstate_f_extent)

        c__pdafomi_gather_obsstate(&i_obs, obsstate_p_ptr, obsstate_f_ptr)

    return obsstate_f_np


def get_interp_coeff_tri(double [::1,:] gpc, double [::1] oc,
    double [::1] icoeff):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    gpc : ndarray[np.float64, ndim=2]
        Coordinates of grid points; dim(3,2)
        Array shape: (:,:)
    oc : ndarray[np.float64, ndim=1]
        Coordinates of observation; dim(2)
        Array shape: (:)
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients; dim(3)
        Array shape: (:)

    Returns
    -------
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients; dim(3)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 icoeff_cfi
    cdef CFI_cdesc_t *icoeff_ptr = <CFI_cdesc_t *> &icoeff_cfi
    cdef size_t icoeff_nbytes = icoeff.nbytes
    cdef CFI_index_t icoeff_extent[1]
    icoeff_extent[0] = icoeff.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] icoeff_np = np.asarray(icoeff, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 gpc_cfi
    cdef CFI_cdesc_t *gpc_ptr = <CFI_cdesc_t *> &gpc_cfi
    cdef size_t gpc_nbytes = gpc.nbytes
    cdef CFI_index_t gpc_extent[2]
    gpc_extent[0] = gpc.shape[0]
    gpc_extent[1] = gpc.shape[1]
    cdef CFI_cdesc_rank1 oc_cfi
    cdef CFI_cdesc_t *oc_ptr = <CFI_cdesc_t *> &oc_cfi
    cdef size_t oc_nbytes = oc.nbytes
    cdef CFI_index_t oc_extent[1]
    oc_extent[0] = oc.shape[0]
    with nogil:
        CFI_establish(gpc_ptr, &gpc[0,0], CFI_attribute_other,
                      CFI_type_double , gpc_nbytes, 2, gpc_extent)

        CFI_establish(oc_ptr, &oc[0], CFI_attribute_other,
                      CFI_type_double , oc_nbytes, 1, oc_extent)

        CFI_establish(icoeff_ptr, &icoeff[0], CFI_attribute_other,
                      CFI_type_double , icoeff_nbytes, 1, icoeff_extent)

        c__pdafomi_get_interp_coeff_tri(gpc_ptr, oc_ptr, icoeff_ptr)

    return icoeff_np


def get_interp_coeff_lin1d(double [::1] gpc, double  oc, double [::1] icoeff):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    gpc : ndarray[np.float64, ndim=1]
        Coordinates of grid points (dim=2)
        Array shape: (:)
    oc : double
        Coordinates of observation
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients (dim=2)
        Array shape: (:)

    Returns
    -------
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients (dim=2)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 icoeff_cfi
    cdef CFI_cdesc_t *icoeff_ptr = <CFI_cdesc_t *> &icoeff_cfi
    cdef size_t icoeff_nbytes = icoeff.nbytes
    cdef CFI_index_t icoeff_extent[1]
    icoeff_extent[0] = icoeff.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] icoeff_np = np.asarray(icoeff, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 gpc_cfi
    cdef CFI_cdesc_t *gpc_ptr = <CFI_cdesc_t *> &gpc_cfi
    cdef size_t gpc_nbytes = gpc.nbytes
    cdef CFI_index_t gpc_extent[1]
    gpc_extent[0] = gpc.shape[0]
    with nogil:
        CFI_establish(gpc_ptr, &gpc[0], CFI_attribute_other,
                      CFI_type_double , gpc_nbytes, 1, gpc_extent)

        CFI_establish(icoeff_ptr, &icoeff[0], CFI_attribute_other,
                      CFI_type_double , icoeff_nbytes, 1, icoeff_extent)

        c__pdafomi_get_interp_coeff_lin1d(gpc_ptr, &oc, icoeff_ptr)

    return icoeff_np


def get_interp_coeff_lin(int  num_gp, int  n_dim, double [::1,:] gpc,
    double [::1] oc, double [::1] icoeff):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    num_gp : int
        Length of icoeff
    n_dim : int
        Number of dimensions in interpolation
    gpc : ndarray[np.float64, ndim=2]
        Coordinates of grid points
        Array shape: (:,:)
    oc : ndarray[np.float64, ndim=1]
        Coordinates of observation
        Array shape: (:)
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients (num_gp)
        Array shape: (:)

    Returns
    -------
    icoeff : ndarray[np.float64, ndim=1]
        Interpolation coefficients (num_gp)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 icoeff_cfi
    cdef CFI_cdesc_t *icoeff_ptr = <CFI_cdesc_t *> &icoeff_cfi
    cdef size_t icoeff_nbytes = icoeff.nbytes
    cdef CFI_index_t icoeff_extent[1]
    icoeff_extent[0] = icoeff.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] icoeff_np = np.asarray(icoeff, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 gpc_cfi
    cdef CFI_cdesc_t *gpc_ptr = <CFI_cdesc_t *> &gpc_cfi
    cdef size_t gpc_nbytes = gpc.nbytes
    cdef CFI_index_t gpc_extent[2]
    gpc_extent[0] = gpc.shape[0]
    gpc_extent[1] = gpc.shape[1]
    cdef CFI_cdesc_rank1 oc_cfi
    cdef CFI_cdesc_t *oc_ptr = <CFI_cdesc_t *> &oc_cfi
    cdef size_t oc_nbytes = oc.nbytes
    cdef CFI_index_t oc_extent[1]
    oc_extent[0] = oc.shape[0]
    with nogil:
        CFI_establish(gpc_ptr, &gpc[0,0], CFI_attribute_other,
                      CFI_type_double , gpc_nbytes, 2, gpc_extent)

        CFI_establish(oc_ptr, &oc[0], CFI_attribute_other,
                      CFI_type_double , oc_nbytes, 1, oc_extent)

        CFI_establish(icoeff_ptr, &icoeff[0], CFI_attribute_other,
                      CFI_type_double , icoeff_nbytes, 1, icoeff_extent)

        c__pdafomi_get_interp_coeff_lin(&num_gp, &n_dim, gpc_ptr, oc_ptr,
                                        icoeff_ptr)

    return icoeff_np


def init_dim_obs_l_iso(int  i_obs, double [::1] coords_l, int  locweight,
    double  cradius, double  sradius, int  cnt_obs_l_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweight : int
        Type of localization function
    cradius : double
        Localization cut-off radius
    sradius : double
        Support radius of localization function
    cnt_obs_l_all : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l_all : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef CFI_cdesc_t *coords_l_ptr = <CFI_cdesc_t *> &coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    with nogil:
        CFI_establish(coords_l_ptr, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

        c__pdafomi_init_dim_obs_l_iso(&i_obs, coords_l_ptr, &locweight,
                                      &cradius, &sradius, &cnt_obs_l_all)

    return cnt_obs_l_all


def init_dim_obs_l_noniso(int  i_obs, double [::1] coords_l,
    int  locweight, double [::1] cradius, double [::1] sradius,
    int  cnt_obs_l_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweight : int
        Type of localization function
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    cnt_obs_l_all : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l_all : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef CFI_cdesc_t *coords_l_ptr = <CFI_cdesc_t *> &coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef CFI_cdesc_t *cradius_ptr = <CFI_cdesc_t *> &cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef CFI_cdesc_t *sradius_ptr = <CFI_cdesc_t *> &sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    with nogil:
        CFI_establish(coords_l_ptr, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

        CFI_establish(cradius_ptr, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

        CFI_establish(sradius_ptr, &sradius[0], CFI_attribute_other,
                      CFI_type_double , sradius_nbytes, 1, sradius_extent)

        c__pdafomi_init_dim_obs_l_noniso(&i_obs, coords_l_ptr, &locweight,
                                         cradius_ptr, sradius_ptr,
                                         &cnt_obs_l_all)

    return cnt_obs_l_all


def init_dim_obs_l_noniso_locweights(int  i_obs, double [::1] coords_l,
    int [::1] locweights, double [::1] cradius, double [::1] sradius,
    int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweights : ndarray[np.intc, ndim=1]
        Types of localization function
        Array shape: (:)
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef CFI_cdesc_t *coords_l_ptr = <CFI_cdesc_t *> &coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    cdef CFI_cdesc_rank1 locweights_cfi
    cdef CFI_cdesc_t *locweights_ptr = <CFI_cdesc_t *> &locweights_cfi
    cdef size_t locweights_nbytes = locweights.nbytes
    cdef CFI_index_t locweights_extent[1]
    locweights_extent[0] = locweights.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef CFI_cdesc_t *cradius_ptr = <CFI_cdesc_t *> &cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef CFI_cdesc_t *sradius_ptr = <CFI_cdesc_t *> &sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    with nogil:
        CFI_establish(coords_l_ptr, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

        CFI_establish(locweights_ptr, &locweights[0], CFI_attribute_other,
                      CFI_type_int , locweights_nbytes, 1, locweights_extent)

        CFI_establish(cradius_ptr, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

        CFI_establish(sradius_ptr, &sradius[0], CFI_attribute_other,
                      CFI_type_double , sradius_nbytes, 1, sradius_extent)

        c__pdafomi_init_dim_obs_l_noniso_locweights(&i_obs, coords_l_ptr,
                                                    locweights_ptr,
                                                    cradius_ptr,
                                                    sradius_ptr, &cnt_obs_l)

    return cnt_obs_l


def obs_op_gridpoint(int  i_obs, double [::1] state_p, double [::1] obs_f_all):
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
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    with nogil:
        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        c__pdafomi_obs_op_gridpoint(&i_obs, state_p_ptr, obs_f_all_ptr)

    return obs_f_all_np


def obs_op_gridavg(int  i_obs, int  nrows, double [::1] state_p,
    double [::1] obs_f_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nrows : int
        Number of values to be averaged
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
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    with nogil:
        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        c__pdafomi_obs_op_gridavg(&i_obs, &nrows, state_p_ptr, obs_f_all_ptr)

    return obs_f_all_np


def obs_op_interp_lin(int  i_obs, int  nrows, double [::1] state_p,
    double [::1] obs_f_all):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nrows : int
        Number of values to be averaged
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
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    with nogil:
        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        c__pdafomi_obs_op_interp_lin(&i_obs, &nrows, state_p_ptr, obs_f_all_ptr)

    return obs_f_all_np


def obs_op_adj_gridpoint(int  i_obs, double [::1] obs_f_all,
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
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    with nogil:
        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        c__pdafomi_obs_op_adj_gridpoint(&i_obs, obs_f_all_ptr, state_p_ptr)

    return state_p_np


def obs_op_adj_gridavg(int  i_obs, int  nrows, double [::1] obs_f_all,
    double [::1] state_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nrows : int
        Number of values to be averaged
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    with nogil:
        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        c__pdafomi_obs_op_adj_gridavg(&i_obs, &nrows, obs_f_all_ptr,
                                      state_p_ptr)

    return state_p_np


def obs_op_adj_interp_lin(int  i_obs, int  nrows, double [::1] obs_f_all,
    double [::1] state_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nrows : int
        Number of values to be averaged
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types (nobs_f_all)
        Array shape: (:)
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
        PE-local model state (dim_p)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 state_p_cfi
    cdef CFI_cdesc_t *state_p_ptr = <CFI_cdesc_t *> &state_p_cfi
    cdef size_t state_p_nbytes = state_p.nbytes
    cdef CFI_index_t state_p_extent[1]
    state_p_extent[0] = state_p.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    with nogil:
        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        CFI_establish(state_p_ptr, &state_p[0], CFI_attribute_other,
                      CFI_type_double , state_p_nbytes, 1, state_p_extent)

        c__pdafomi_obs_op_adj_interp_lin(&i_obs, &nrows, obs_f_all_ptr,
                                         state_p_ptr)

    return state_p_np


def observation_localization_weights(int  i_obs, int  ncols,
    double [::1,:] a_l, int dim_obs_l, int  verbose):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    ncols : int
        Rank of initial covariance matrix
    a_l : ndarray[np.float64, ndim=2]
        Input matrix (thisobs_l%dim_obs_l, ncols)
        Array shape: (:, :)
    dim_obs_l : int
        Dimension of local observation vector of the i_obs-th observation type
    verbose : int
        Verbosity flag

    Returns
    -------
    weight : ndarray[np.float64, ndim=1]
        > Localization weights
        Array shape: (thisobs_l(i_obs)%dim_obs_l)
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] weight_np = np.zeros((dim_obs_l), dtype=np.float64, order="F")
    cdef double [::1] weight = weight_np
    cdef CFI_cdesc_rank2 a_l_cfi
    cdef CFI_cdesc_t *a_l_ptr = <CFI_cdesc_t *> &a_l_cfi
    cdef size_t a_l_nbytes = a_l.nbytes
    cdef CFI_index_t a_l_extent[2]
    a_l_extent[0] = a_l.shape[0]
    a_l_extent[1] = a_l.shape[1]
    with nogil:
        CFI_establish(a_l_ptr, &a_l[0,0], CFI_attribute_other,
                      CFI_type_double , a_l_nbytes, 2, a_l_extent)

        c__pdafomi_observation_localization_weights(&i_obs, &ncols,
                                                    a_l_ptr, &weight[0],
                                                    &verbose)

    return weight_np


def set_debug_flag(int  debugval):
    """Activate the debug output of the PDAF.

    Starting from the use of this function,
    the debug infomation is sent to screen output.
    The screen output end when the debug flag is
    set to 0 by this function.

    For the sake of simplicity,
    we recommend using debugging output for
    a single local domain, e.g.
    `if domain_p == 1: pyPDAF.PDAF.set_debug_flag(1)`

    Parameters
    ----------
    debugval : int
        Value for debugging flag

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_debug_flag(&debugval)



def set_dim_obs_l(int  i_obs, int  cnt_obs_l_all, int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    cnt_obs_l_all : int
        Local dimension of observation vector over all obs. types
    cnt_obs_l : int
        Local dimension of single observation type vector

    Returns
    -------
    cnt_obs_l_all : int
        Local dimension of observation vector over all obs. types
    cnt_obs_l : int
        Local dimension of single observation type vector
    """
    with nogil:
        c__pdafomi_set_dim_obs_l(&i_obs, &cnt_obs_l_all, &cnt_obs_l)

    return cnt_obs_l_all, cnt_obs_l


def set_localization(int  i_obs, double  cradius, double  sradius,
    int  locweight):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    cradius : double
        Localization cut-off radius
    sradius : double
        Support radius of localization function
    locweight : int
        Type of localization function

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_localization(&i_obs, &cradius, &sradius, &locweight)



def set_localization_noniso(int  i_obs, int  nradii, double [::1] cradius,
    double [::1] sradius, int  locweight, int  locweight_v):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    nradii : int
        Number of radii to consider for localization
    cradius : ndarray[np.float64, ndim=1]
        Localization cut-off radius
        Array shape: (nradii)
    sradius : ndarray[np.float64, ndim=1]
        Support radius of localization function
        Array shape: (nradii)
    locweight : int
        Type of localization function
    locweight_v : int
        Type of localization function in vertical direction (only for nradii=3)

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_localization_noniso(&i_obs, &nradii, &cradius[0],
                                           &sradius[0], &locweight,
                                           &locweight_v)



def set_localize_covar_iso(int  i_obs, int  dim, int  ncoords,
    double [::1,:] coords, int  locweight, double  cradius, double  sradius):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    ncoords : int
        number of coordinate directions
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    locweight : int
        Localization weight type
    cradius : double
        localization radius
    sradius : double
        support radius for weight functions

    Returns
    -------
    """
    cdef CFI_cdesc_rank2 coords_cfi
    cdef CFI_cdesc_t *coords_ptr = <CFI_cdesc_t *> &coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    with nogil:
        CFI_establish(coords_ptr, &coords[0,0], CFI_attribute_other,
                      CFI_type_double , coords_nbytes, 2, coords_extent)

        c__pdafomi_set_localize_covar_iso(&i_obs, &dim, &ncoords,
                                          coords_ptr, &locweight, &cradius,
                                          &sradius)



def set_localize_covar_noniso(int  i_obs, int  dim, int  ncoords,
    double [::1,:] coords, int  locweight, double [::1] cradius,
    double [::1] sradius):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    ncoords : int
        number of coordinate directions
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    locweight : int
        Localization weight type
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)

    Returns
    -------
    """
    cdef CFI_cdesc_rank2 coords_cfi
    cdef CFI_cdesc_t *coords_ptr = <CFI_cdesc_t *> &coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef CFI_cdesc_t *cradius_ptr = <CFI_cdesc_t *> &cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef CFI_cdesc_t *sradius_ptr = <CFI_cdesc_t *> &sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    with nogil:
        CFI_establish(coords_ptr, &coords[0,0], CFI_attribute_other,
                      CFI_type_double , coords_nbytes, 2, coords_extent)

        CFI_establish(cradius_ptr, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

        CFI_establish(sradius_ptr, &sradius[0], CFI_attribute_other,
                      CFI_type_double , sradius_nbytes, 1, sradius_extent)

        c__pdafomi_set_localize_covar_noniso(&i_obs, &dim, &ncoords,
                                             coords_ptr, &locweight,
                                             cradius_ptr, sradius_ptr)



def set_localize_covar_noniso_locweights(int  i_obs, int  dim,
    int  ncoords, double [::1,:] coords, int [::1] locweights,
    double [::1] cradius, double [::1] sradius):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    ncoords : int
        number of coordinate directions
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    locweights : ndarray[np.intc, ndim=1]
        Types of localization function
        Array shape: (:)
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)

    Returns
    -------
    """
    cdef CFI_cdesc_rank2 coords_cfi
    cdef CFI_cdesc_t *coords_ptr = <CFI_cdesc_t *> &coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    cdef CFI_cdesc_rank1 locweights_cfi
    cdef CFI_cdesc_t *locweights_ptr = <CFI_cdesc_t *> &locweights_cfi
    cdef size_t locweights_nbytes = locweights.nbytes
    cdef CFI_index_t locweights_extent[1]
    locweights_extent[0] = locweights.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef CFI_cdesc_t *cradius_ptr = <CFI_cdesc_t *> &cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef CFI_cdesc_t *sradius_ptr = <CFI_cdesc_t *> &sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    with nogil:
        CFI_establish(coords_ptr, &coords[0,0], CFI_attribute_other,
                      CFI_type_double , coords_nbytes, 2, coords_extent)

        CFI_establish(locweights_ptr, &locweights[0], CFI_attribute_other,
                      CFI_type_int , locweights_nbytes, 1, locweights_extent)

        CFI_establish(cradius_ptr, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

        CFI_establish(sradius_ptr, &sradius[0], CFI_attribute_other,
                      CFI_type_double , sradius_nbytes, 1, sradius_extent)

        c__pdafomi_set_localize_covar_noniso_locweights(&i_obs, &dim,
                                                        &ncoords,
                                                        coords_ptr,
                                                        locweights_ptr,
                                                        cradius_ptr,
                                                        sradius_ptr)



def set_obs_diag(int  diag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    diag : int
        Value for observation diagnostics mode

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_obs_diag(&diag)



def set_domain_limits(double [::1,:] lim_coords):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    lim_coords : ndarray[np.float64, ndim=2]
        geographic coordinate array (1: longitude, 2: latitude)
        Array shape: (2,2)

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_domain_limits(&lim_coords[0,0])



def get_domain_limits_unstr(int  npoints_p, double [::1,:] coords_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    npoints_p : int
        number of process-local grid points
    coords_p : ndarray[np.float64, ndim=2]
        geographic coordinate array (row 1: longitude, 2: latitude)
        Array shape: (:,:)

    Returns
    -------
    """
    cdef CFI_cdesc_rank2 coords_p_cfi
    cdef CFI_cdesc_t *coords_p_ptr = <CFI_cdesc_t *> &coords_p_cfi
    cdef size_t coords_p_nbytes = coords_p.nbytes
    cdef CFI_index_t coords_p_extent[2]
    coords_p_extent[0] = coords_p.shape[0]
    coords_p_extent[1] = coords_p.shape[1]
    with nogil:
        CFI_establish(coords_p_ptr, &coords_p[0,0], CFI_attribute_other,
                      CFI_type_double , coords_p_nbytes, 2, coords_p_extent)

        c__pdafomi_get_domain_limits_unstr(&npoints_p, coords_p_ptr)



def store_obs_l_index(int  i_obs, int  idx, int  id_obs_l,
    double  distance, double  cradius_l, double  sradius_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    idx : int
        Element of local observation array to be filled
    id_obs_l : int
        Index of local observation in full observation array
    distance : double
        Distance between local analysis domain and observation
    cradius_l : double
        cut-off radius for this local observation
    sradius_l : double
        support radius for this local observation

    Returns
    -------
    """
    with nogil:
        c__pdafomi_store_obs_l_index(&i_obs, &idx, &id_obs_l, &distance,
                                     &cradius_l, &sradius_l)



def store_obs_l_index_vdist(int  i_obs, int  idx, int  id_obs_l,
    double  distance, double  cradius_l, double  sradius_l, double  vdist):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    idx : int
        Element of local observation array to be filled
    id_obs_l : int
        Index of local observation in full observation array
    distance : double
        Distance between local analysis domain and observation
    cradius_l : double
        cut-off radius for this local observation
    sradius_l : double
        support radius for this local observation
    vdist : double
        support radius in vertical direction for 2+1D factorized localization

    Returns
    -------
    """
    with nogil:
        c__pdafomi_store_obs_l_index_vdist(&i_obs, &idx, &id_obs_l,
                                           &distance, &cradius_l,
                                           &sradius_l, &vdist)



