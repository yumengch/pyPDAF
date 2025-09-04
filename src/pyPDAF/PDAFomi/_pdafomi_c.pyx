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

def init(int  n_obs):
    r"""Allocating an array of `obs_f` derived types instances.

    This function initialises the number of observation types,
    which should be called at the start of the DA system
    after :func:`pyPDAF.PDAF.init`.

    Parameters
    ----------
    n_obs : int
        number of observations
    """
    with nogil:
        c__pdafomi_init(&n_obs)



def init_local():
    r"""Allocating an array of `obs_l` derived types instances.

    This function initialises the number of observation types
    for each local analysis domain,
    which should be called at the start of the local analysis loop
    in :func:`py__init_dim_obs_l_pdaf`.

    """
    with nogil:
        c__pdafomi_init_local()



def check_error(int  flag):
    r"""This function returns the value of the PDAF-OMI internal error flag.

    Since PDAF-OMI executes internal routines in which errors could occur due to
    an inconsistent configuration of the observations. Directly returning
    an error flag as a subroutine argument is not always possible.
    For this reason there is this separate routine to check for the error flag.

    The errors that are checked by PDAF-OMI relate to the configuration of
    the observations, e.g. it is checked whether some dimensions are consistent.

    This can be useful if an error occurred PDAF-OMI also prints an error message,
    but it does not stop the program.


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


def gather_obs(int  i_obs, int  dim_obs_p, double[::1] obs_p,
    double[::1] ivar_obs_p, double[::1,:] ocoord_p, int  ncoord,
    double  lradius):
    r"""Gather the dimension of a given type of observation across
    multiple local domains/filter processors.

    This function can be used in the user-supplied function of
    :func:`py__init_dim_obs_f_pdaf`.

    This function does three things:
        1. Receiving observation dimension on each local process.
        2. Gather the total dimension of given observation type
           across local process and the displacement of PE-local
           observations relative to the total observation vector
        3. Set the observation vector, observation coordinates,
           the inverse of the observation variance, and localisation
           radius for this observation type.

    Parameters
    ----------
    i_obs : int
        index of observation type
    dim_obs_p: int
        PE-local dimension of observation vector
    obs_p : ndarray[tuple[dim_obs_p, ...], np.float64]
        PE-local observation vector
        The array dimension `dim_obs_p` is dimension of PE-local observation vector
    ivar_obs_p : ndarray[tuple[dim_obs_p, ...], np.float64]
        PE-local inverse of observation error variance
        The array dimension `dim_obs_p` is dimension of PE-local observation vector
    ocoord_p : ndarray[tuple[thisobs(i_obs)%ncoord, dim_obs_p, ...], np.float64]
        pe-local observation coordinates
        The 1st-th dimension dim_obs_p is dimension of PE-local observation vector
    ncoord: int
        Number of rows of coordinate array
    lradius : float
        localization radius

    Returns
    -------
    dim_obs : int
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
    r"""This function is used to implement custom observation operators.

    This function is used inside a custom observation operator.
    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Implementingyourownobservationoperator>`_

    Parameters
    ----------
    i_obs : int
        index of observations
    obsstate_p : ndarray[tuple[thisobs(i_obs)%dim_obs_p, ...], np.float64]
        Vector of process-local observed state
    obsstate_f : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed vector for all types
        The array dimension `nobs_f_all` is dimension of the observation

    Returns
    -------
    obsstate_f : ndarray[tuple[nobs_f_all, ...], np.float64]
         Full observed vector for all types

        The array dimension `nobs_f_all` is dimension of the observation
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
    r"""The coefficient for linear interpolation in 2D on unstructure triangular grid.

    The resulting coefficient is used in :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    This function is for triangular model grid interpolation coefficients determined as barycentric coordinates.

    Parameters
    ----------
    gpc : ndarray[tuple[3, 2, ...], np.float64]
        Coordinates of grid points with dimension of (3, 2).
        3 grid points surrounding the observation;
        each containing lon and lat coordinates.
        The order of the grid points in gcoords has to
        be consistent with the order of the indices specified in
        `id_obs_p` of `obs_f`.
    oc : ndarray[tuple[2, ...], np.float64]
        Coordinates of observation (targeted location); dim(2)
    icoeff : ndarray[tuple[3, ...], np.float64]
        Interpolation coefficients; dim(3)

    Returns
    -------
    icoeff : ndarray[tuple[3, ...], np.float64]
         Interpolation coefficients; dim(3)

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
    r"""The coefficient for linear interpolation in 1D.

    The resulting coefficient is used in :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    Parameters
    ----------
    gpc : ndarray[tuple[2, ...], np.float64]
        Coordinates of grid points surrounding the observations (dim=2)
    oc : float
        Coordinates of observation (targeted location)
    icoeff : ndarray[tuple[2, ...], np.float64]
        Interpolation coefficients (dim=2)

    Returns
    -------
    icoeff : ndarray[tuple[2, ...], np.float64]
         Interpolation coefficients (dim=2)

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
    r"""The coefficient for linear interpolation up to 3D.

    The resulting coefficient is used in :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    See introduction in `relevant PDAF-OMI wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin>`_

    Parameters
    ----------
    gpc : ndarray[tuple[num_gp, n_dim, ...], np.float64]
        Coordinates of grid points
        The order of the grid points in gcoords has to
        be consistent with the order of the indices specified in
        `id_obs_p` of `obs_f`.
        The 1st-th dimension num_gp is Length of icoeff
        The 2nd-th dimension n_dim is Number of dimensions in interpolation
    oc : ndarray[tuple[n_dim, ...], np.float64]
        Coordinates of observation
        The array dimension `n_dim` is Number of dimensions in interpolation
    icoeff : ndarray[tuple[num_gp, ...], np.float64]
        Interpolation coefficients (num_gp)
        The array dimension `num_gp` is Length of icoeff

    Returns
    -------
    icoeff : ndarray[tuple[num_gp, ...], np.float64]
         Interpolation coefficients (num_gp)

        The array dimension `num_gp` is Length of icoeff
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
    r"""Initialize the observation information corresponding to an isotropic local analysis domain.

    One can set localization parameters, like the localization radius, for each observation type.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord, ...], np.float64]
        Coordinates of current analysis domain
        The array dimension `ncoord` is number of coordinate dimension
    locweight : int
        Types of localization function
        0) unit weight; 1) exponential; 2) 5-th order polynomial;
        3) 5-th order polynomial with regulatioin using mean variance;
        4) 5-th order polynomial with regulatioin using variance of single observation point;
    cradius : float
        Vector of localization cut-off radii; observation weight=0 if distance > cradius
    sradius : float
        Vector of support radii of localization function.
        It has no impact if locweight=0; 	weight = exp(-d / sradius) if locweight=1;
        weight = 0 if d >= sradius else f(sradius, distance) if locweight in [2,3,4].
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
    r"""Initialize the observation information corresponding to a non-isotropic local analysis domain.

    One can set localization parameters, like the localization radius, for each observation type.

    Here, each dimension can use a different localisation radius.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_
    as well as `non-isotropic localisation page <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#Non-isotropiclocalization>`_

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord, ...], np.float64]
        Coordinates of current analysis domain
        The array dimension `ncoord` is number of coordinate dimension
    locweight : int
        Types of localization function
        0) unit weight; 1) exponential; 2) 5-th order polynomial;
        3) 5-th order polynomial with regulatioin using mean variance;
        4) 5-th order polynomial with regulatioin using variance of single observation point;
    cradius : ndarray[tuple[ncoord, ...], np.float64]
        Vector of localization cut-off radii; observation weight=0 if distance > cradius
        The array dimension `ncoord` is number of coordinate dimension
    sradius : ndarray[tuple[ncoord, ...], np.float64]
        Vector of support radii of localization function.
        It has no impact if locweight=0; 	weight = exp(-d / sradius) if locweight=1;
        weight = 0 if d >= sradius else f(sradius, distance) if locweight in [2,3,4].
        The array dimension `ncoord` is number of coordinate dimension
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
    r"""Initialize the observation information corresponding to a non-isotropic local analysis domain.

    One can set localization parameters, like the localization radius, for each observation type.

    Here, each dimension can use a different localisation radius and a different
    localisation weight.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_
    as well as `non-isotropic localisation page <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#Non-isotropiclocalization>`_

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord, ...], np.float64]
        Coordinates of current analysis domain
        The array dimension `ncoord` is number of coordinate dimension
    locweights : ndarray[tuple[2, ...], np.intc]
        Types of localization function
        0) unit weight; 1) exponential; 2) 5-th order polynomial;
        3) 5-th order polynomial with regulatioin using mean variance;
        4) 5-th order polynomial with regulatioin using variance of single observation point;
        The first dimension is horizontal weight function and the second is the vertical function
    cradius : ndarray[tuple[ncoord, ...], np.float64]
        Vector of localization cut-off radii for each dimension; observation weight=0 if distance > cradius
        The array dimension `ncoord` is number of coordinate dimension
    sradius : ndarray[tuple[ncoord, ...], np.float64]
        Vector of support radii of localization function for each dimension.
        It has no impact if locweight=0; 	weight = exp(-d / sradius) if locweight=1;
        weight = 0 if d >= sradius else f(sradius, distance) if locweight in [2,3,4].
        The array dimension `ncoord` is number of coordinate dimension
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
    r"""A (partial) identity observation operator

    This observation operator is used
    when observations and model use the same grid.

    The observations operator selects state vectors
    where observations are present based on properties given
    in `obs_f`, e.g., `id_obs_p`.

    The function is used in
    the user-supplied function :func:`py__obs_op_pdaf`.

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
         Full observed state for all observation types (nobs_f_all)

        The array dimension `nobs_f_all` is dimension of the observation
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
    r"""Observation operator that average values on given model grid points.

    The averaged model grid points are specified in `id_obs_p` property of `obs_f`,
    which can be set in :func:`pyPDAF.PDAF.omi_set_id_obs_p`.

    The function is used in the user-supplied function `py__obs_op_pdaf`.

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
         Full observed state for all observation types (nobs_f_all)

        The array dimension `nobs_f_all` is dimension of the observation
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

def obs_op_extern(int  i_obs, double [::1] ostate_p,
    double [::1] obs_f_all):
    """Observation operator for given observed model state.

    Application of observation operator for the case that
    a user performs the actual observation operator elsewhere,
    e.g., directly in the model during the forecast or offline.
    For this case, the user can provide the observed model state
    to this routine and it will just call :func:`pyPDAF.PDAFomi.gather_obsstate`,
    to obtain the full observed vector `obs_f_all`.

    This has to be called in all filter processes.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    ostate_p : ndarray[np.float64, ndim=1]
        PE-local observed model state
        Array shape: (:)
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed model state for all observation types
        Array shape: (:)

    Returns
    -------
    obs_f_all : ndarray[np.float64, ndim=1]
        Full observed state for all observation types
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 obs_f_all_cfi
    cdef CFI_cdesc_t *obs_f_all_ptr = <CFI_cdesc_t *> &obs_f_all_cfi
    cdef size_t obs_f_all_nbytes = obs_f_all.nbytes
    cdef CFI_index_t obs_f_all_extent[1]
    obs_f_all_extent[0] = obs_f_all.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_f_all_np = np.asarray(obs_f_all, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 ostate_p_cfi
    cdef CFI_cdesc_t *ostate_p_ptr = <CFI_cdesc_t *> &ostate_p_cfi
    cdef size_t ostate_p_nbytes = ostate_p.nbytes
    cdef CFI_index_t ostate_p_extent[1]
    ostate_p_extent[0] = ostate_p.shape[0]
    with nogil:
        CFI_establish(ostate_p_ptr, &ostate_p[0], CFI_attribute_other,
                      CFI_type_double , ostate_p_nbytes, 1, ostate_p_extent)

        CFI_establish(obs_f_all_ptr, &obs_f_all[0], CFI_attribute_other,
                      CFI_type_double , obs_f_all_nbytes, 1, obs_f_all_extent)

        c__pdafomi_obs_op_extern(&i_obs, ostate_p_ptr, obs_f_all_ptr)

    return obs_f_all_np


def obs_op_interp_lin(int  i_obs, int  nrows, double [::1] state_p,
    double [::1] obs_f_all):
    r"""Observation operator that linearly interpolates model grid values to observation location.

    The grid points used by linear interpolation is specified in `id_obs_p` of `obs_f`,
    which can be set by :func:`pyPDAF.PDAFomi.set_id_obs_p`.

    The function also requires `icoeff_p` attribute of `obs_f`,
    which can be set by :func:`pyPDAF.PDAFomi.set_icoeff_p`

    The interpolation coefficient can be obtained by :func:`pyPDAF.PDAFomi.get_interp_coeff_lin1D`,
    :func:`pyPDAF.PDAFomi.get_interp_coeff_lin`, and
    :func:`pyPDAF.PDAFomi.get_interp_coeff_tri`

    The details of interpolation setup can be found at
    `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_

    The function is used in the user-supplied function `py__obs_op_pdaf`.

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation

    Returns
    -------
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
         Full observed state for all observation types (nobs_f_all)

        The array dimension `nobs_f_all` is dimension of the observation
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
    r"""The adjoint observation operator of :func:`pyPDAF.PDAFomi.obs_op_gridpoint`.

    This function performs :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`,
    where :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{y}` is any state in observation space.

    Here :math:`\mathbf{H}` is :func:`pyPDAF.PDAFomi.obs_op_gridpoint`.

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
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
    r"""The adjoint observation operator of :func:`pyPDAF.PDAFomi.obs_op_gridavg`.

    This function performs :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`,
    where :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{y}` is any state in observation space.

    Here :math:`\mathbf{H}` is :func:`pyPDAF.PDAFomi.obs_op_gridavg`.

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
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
    r"""The adjoint observation operator of :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    This function performs :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`,
    where :math:`\mathbf{H}` is the observation operator and
    :math:`\mathbf{y}` is any state in observation space.

    Here :math:`\mathbf{H}` is :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        Number of values to be averaged
    obs_f_all : ndarray[tuple[nobs_f_all, ...], np.float64]
        Full observed state for all observation types (nobs_f_all)
        The array dimension `nobs_f_all` is dimension of the observation
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        :math:`\mathbf{H}^\mathrm{T}\mathbf{y}`
        PE-local model state (dim_p)
        The array dimension `dim_p` is dimension of model state
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
    r"""Returns a vector of observation localisation weights.

    The weights are based on specifications given by localisation setups and
    observation coordinates in OMI. This function is used in the case of
    non-diagonal observation error covariance matrix where one has to perform
    localisation in user-supplied functions, e.g.,
    :func:`pyPDAF.c__prodrinva_pdaf` or :func:`pyPDAF.c__likelihood_l_pdaf`.

    Here, `a_l` is typically the input array in :func:`pyPDAF.c__prodrinva_pdaf`.

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
    """Activate the debug output of the PDAFomi.

    Starting from the use of this function,
    the debug infomation is sent to screen output.
    The screen output end when the debug flag is
    set to 0 by this function.

    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_debugging>`_

    Parameters
    ----------
    debugval : int
        Value for debugging flag
    """
    with nogil:
        c__pdafomi_set_debug_flag(&debugval)



def set_dim_obs_l(int  i_obs, int  cnt_obs_l_all, int  cnt_obs_l):
    """Stores the local number of observations for OMI-internal initialisations.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_

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
    r"""Stores the isotropic localization parameters (cradius, sradius, locweight) in OMI.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_


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
    """
    with nogil:
        c__pdafomi_set_localization(&i_obs, &cradius, &sradius, &locweight)



def set_localization_noniso(int  i_obs, int  nradii, double [::1] cradius,
    double [::1] sradius, int  locweight, int  locweight_v):
    r"""Stores the non-isotropic localization parameters (cradius, sradius, locweight) in OMI.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_



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
    r"""Initialise local observation information for isotropic covariance localisation.

    This is only used in stochastic EnKF. This is called in user-supplied functions
    :func:`pyPDAF.c__init_dim_obs_pdaf`.

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
    r"""Initialise local observation information for non-isotropic covariance localisation.

    This is only used in stochastic EnKF. This is called in user-supplied functions
    :func:`pyPDAF.c__init_dim_obs_pdaf`.

    Here, localisation radii differ for each spatial dimension.

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
    r"""Initialise local observation information for non-isotropic covariance localisation.

    This is only used in stochastic EnKF. This is called in user-supplied functions
    :func:`pyPDAF.c__init_dim_obs_pdaf`.

    Here, both weighting function and localisation radii differ for each spatial dimension.

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
    """Activate or deactivate the observation diagnostics.

    By default, observation diagnostics are activated that stores
    additional information for diagnostics.
    However, as this functionality increases the required memory,
    it might be desirable to deactivate this functionality.

    This function is used deactivate the observation diagnostics.
    Once deactivated, one cannot use diagnostics in :mod:`pyPDAF.PDAFomi.diag`.
    It is also possible to re-activate the observation diagnostics at a later time.

    The function can be called by all processes, but it is sufficient to call it
    for those processes that handle observations, which usually are the filter processes.

    This function can be called after the initialization of PDAF in `pyPDAF.PDAF.init`.


    Parameters
    ----------
    diag : int
        Value for observation diagnostics mode
        - > 0: activates observation diagnostics
        - 0: deactivates observation diagnostics
    """
    with nogil:
        c__pdafomi_set_obs_diag(&diag)



def set_domain_limits(double [::1,:] lim_coords):
    r"""Set the domain limits for domain decomposed local domain.

    This is for the use of :func:`pyPDAF.PDAFomi.set_use_global_obs`.
    Currently, it only supports 2D limitations.

    See `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_set_domain_limits>`_


    Parameters
    ----------
    lim_coords : ndarray[tuple[2, 2, ...], np.float64]
        geographic coordinate array (1: longitude, 2: latitude)
    """
    with nogil:
        c__pdafomi_set_domain_limits(&lim_coords[0,0])



def get_domain_limits_unstr(int  npoints_p, double [::1,:] coords_p):
    r"""Set the domain limits for unstructured domain decomposed local domain.

    This is for the use of :func:`pyPDAF.PDAFomi.set_use_global_obs`.
    Currently, it only supports 2D limitations.

    See also `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_get_domain_limits_unstrc>`_

    Parameters
    ----------
    npoints_p : int
        number of process-local grid points
    coords_p : ndarray[np.float64, ndim=2]
        geographic coordinate array, dimension (2, npoints_p)
        (row 1: longitude, 2: latitude)
        ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)
        Array shape: (:,:)
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
    r"""Save local observation information in PDAF.

    One should check `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_
    before using this function.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    idx : int
        index of the valid local observations in current local analysis domain
    id_obs_l : int
        Index of local observation in full observation array
    distance : double
        Distance between local analysis domain and observation
    cradius_l : double
        cut-off radius for this local observation
    sradius_l : double
        support radius for this local observation
    """
    with nogil:
        c__pdafomi_store_obs_l_index(&i_obs, &idx, &id_obs_l, &distance,
                                     &cradius_l, &sradius_l)



def store_obs_l_index_vdist(int  i_obs, int  idx, int  id_obs_l,
    double  distance, double  cradius_l, double  sradius_l, double  vdist):
    r"""Save local observation information for 2+1D factorized localization in the vertical direction in PDAF.

    One should check `relevant PDAF wiki page <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_
    before using this function.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    idx : int
        index of the valid local observations in current local analysis domain
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
    """
    with nogil:
        c__pdafomi_store_obs_l_index_vdist(&i_obs, &idx, &id_obs_l,
                                           &distance, &cradius_l,
                                           &sradius_l, &vdist)



