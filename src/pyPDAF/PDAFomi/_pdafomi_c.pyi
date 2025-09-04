# pylint: disable=unused-argument
"""Stub file for PDAFomi module
"""
from typing import Tuple
import numpy as np


def init(n_obs:int) -> None:
    r"""Allocating an array of `obs_f` derived types instances.

    This function initialises the number of observation types,
    which should be called at the start of the DA system
    after :func:`pyPDAF.PDAF.init`.

    Parameters
    ----------
    n_obs : int
        number of observations
    """

def init_local() -> None:
    r"""Allocating an array of `obs_l` derived types instances.

    This function initialises the number of observation types
    for each local analysis domain,
    which should be called at the start of the local analysis loop
    in :func:`py__init_dim_obs_l_pdaf`.

    """

def check_error(flag: int) -> int:
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

def gather_obs(i_obs: int, dim_obs_p: int, obs_p: np.ndarray,
               ivar_obs_p: np.ndarray, ocoord_p: np.ndarray,
               ncoord: int, lradius: float) -> int:
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

def gather_obsstate(i_obs: int, obsstate_p: np.ndarray, obsstate_f: np.ndarray) -> np.ndarray:
    r"""This function is used to implement custom observation operators.

    This function is used inside a custom observation operator.
    See also `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Implementingyourownobservationoperator>`_

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

def get_interp_coeff_tri(gpc: np.ndarray, oc: np.ndarray, icoeff: np.ndarray) -> np.ndarray:
    r"""The coefficient for linear interpolation in 2D on unstructure triangular grid.

    The resulting coefficient is used in :func:`pyPDAF.PDAFomi.obs_op_interp_lin`.

    This function is for triangular model grid interpolation coefficients
    determined as barycentric coordinates.

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

def get_interp_coeff_lin1d(gpc: np.ndarray, oc: float, icoeff: np.ndarray) -> np.ndarray:
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

def get_interp_coeff_lin(num_gp: int, n_dim: int, gpc: np.ndarray,
                         oc: np.ndarray, icoeff: np.ndarray) -> np.ndarray:
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

def init_dim_obs_l_iso(i_obs: int, coords_l: np.ndarray, locweight: int,
                       cradius: float, sradius: float, cnt_obs_l_all: int) -> int:
    r"""Initialize the observation information corresponding to an isotropic local analysis domain.

    One can set localization parameters, like the localization radius, for each observation type.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_

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


def init_dim_obs_l_noniso(i_obs: int, coords_l: np.ndarray, locweight: int,
                          cradius: np.ndarray, sradius: np.ndarray,
                          cnt_obs_l_all: int) -> int:
    r"""Initialize the observation information corresponding to a non-isotropic
    local analysis domain.

    One can set localization parameters, like the localization radius, for
    each observation type.

    Here, each dimension can use a different localisation radius.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_
    as well as `non-isotropic localisation page
    <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#Non-isotropiclocalization>`_

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

def init_dim_obs_l_noniso_locweights(i_obs: int, coords_l: np.ndarray,
                                     locweights: np.ndarray, cradius: np.ndarray,
                                     sradius: np.ndarray, cnt_obs_l: int) -> int:
    r"""Initialize the observation information corresponding to a non-isotropic
    local analysis domain.

    One can set localization parameters, like the localization radius, for each
    observation type.

    Here, each dimension can use a different localisation radius and a different
    localisation weight.

    The function has to be called in user-supplied function of
    `init_dim_obs_l_OBTYPE` in each observation module if a
    domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used.

    It initialises the local observation information for PDAF-OMI for a
    single local analysis domain. This is used for isotropic localisation
    where the localisation radius is the same in all directions.

    See also `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_modules_PDAF3#init_dim_obs_l_OBSTYPE>`_
    as well as `non-isotropic localisation page
    <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#Non-isotropiclocalization>`_

    Parameters
    ----------
    i_obs : int
        index of observation type
    coords_l : ndarray[tuple[ncoord, ...], np.float64]
        Coordinates of current analysis domain
        The array dimension `ncoord` is number of coordinate dimension
    locweights : ndarray[tuple[2, ...], np.intc]
        Types of localization function
        - 0) unit weight; 1) exponential; 2) 5-th order polynomial;
        - 3) 5-th order polynomial with regulatioin using mean variance;
        - 4) 5-th order polynomial with regulatioin using variance of single observation point;
        The first dimension is horizontal weight function and the second is the vertical function
    cradius : ndarray[tuple[ncoord, ...], np.float64]
        Vector of localization cut-off radii for each dimension; observation
        weight=0 if distance > cradius
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

def obs_op_gridpoint(i_obs: int, state_p: np.ndarray,
                     obs_f_all: np.ndarray) -> np.ndarray:
    r"""A (partial) identity observation operator

    This observation operator is used
    when observations and model use the same grid.

    The observations operator selects state vectors
    where observations are present based on properties given
    in `obs_f`, e.g., `id_obs_p`.

    The function is used in
    the user-supplied function :func:`pyPDAF.c__obs_op_pdaf`.

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

def obs_op_gridavg(i_obs: int, nrows: int, state_p: np.ndarray,
                   obs_f_all: np.ndarray) -> np.ndarray:
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

def obs_op_extern(i_obs: int, ostate_p: np.ndarray, obs_f_all: np.ndarray) -> np.ndarray:
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

def obs_op_interp_lin(i_obs: int, nrows: int, state_p: np.ndarray,
                      obs_f_all: np.ndarray) -> np.ndarray:
    r"""Observation operator that linearly interpolates model grid values to observation location.

    The grid points used by linear interpolation is specified in `id_obs_p` of `obs_f`,
    which can be set by :func:`pyPDAF.PDAFomi.set_id_obs_p`.

    The function also requires `icoeff_p` attribute of `obs_f`,
    which can be set by :func:`pyPDAF.PDAFomi.set_icoeff_p`

    The interpolation coefficient can be obtained by :func:`pyPDAF.PDAFomi.get_interp_coeff_lin1D`,
    :func:`pyPDAF.PDAFomi.get_interp_coeff_lin`, and
    :func:`pyPDAF.PDAFomi.get_interp_coeff_tri`

    The details of interpolation setup can be found at
    `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_

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

def obs_op_adj_gridpoint(i_obs: int, obs_f_all: np.ndarray,
                         state_p: np.ndarray) -> np.ndarray:
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

def obs_op_adj_gridavg(i_obs: int, nrows: int, obs_f_all: np.ndarray,
                       state_p: np.ndarray) -> np.ndarray:
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

def obs_op_adj_interp_lin(i_obs: int, nrows: int, obs_f_all: np.ndarray,
                          state_p: np.ndarray) -> np.ndarray:
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

def observation_localization_weights(i_obs: int, ncols: int, a_l: np.ndarray,
                                     dim_obs_l: int, verbose: int) -> np.ndarray:
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

def set_debug_flag(debugval: int) -> None:
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

def set_dim_obs_l(i_obs: int, cnt_obs_l_all: int, cnt_obs_l: int) -> Tuple[int, int]:
    """Stores the local number of observations for OMI-internal initialisations.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_

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

def set_localization(i_obs: int, cradius: float, sradius: float, locweight: int) -> None:
    r"""Stores the isotropic localization parameters (cradius, sradius, locweight) in OMI.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_


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

def set_localization_noniso(i_obs: int, nradii: int, cradius: np.ndarray,
                            sradius: np.ndarray, locweight: int, locweight_v: int) -> None:
    r"""Stores the non-isotropic localization parameters (cradius, sradius, locweight) in OMI.

    This is used for alternative to :func:`pyPDF.PDAFomi.init_dim_obs_l`.

    See more details in `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_



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

def set_localize_covar_iso(i_obs: int, dim: int, ncoords: int, coords: np.ndarray,
                           locweight: int, cradius: float, sradius: float) -> None:
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

def set_localize_covar_noniso(i_obs: int, dim: int, ncoords: int,
                              coords: np.ndarray, locweight: int,
                              cradius: np.ndarray, sradius: np.ndarray) -> None:
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

def set_localize_covar_noniso_locweights(i_obs: int, dim: int, ncoords: int,
                                         coords: np.ndarray, locweights: np.ndarray,
                                         cradius: np.ndarray, sradius: np.ndarray) -> None:
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

def set_obs_diag(diag: int) -> None:
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

def set_domain_limits(lim_coords: np.ndarray) -> None:
    r"""Set the domain limits for domain decomposed local domain.

    This is for the use of :func:`pyPDAF.PDAFomi.set_use_global_obs`.
    Currently, it only supports 2D limitations.

    See `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_set_domain_limits>`_


    Parameters
    ----------
    lim_coords : ndarray[tuple[2, 2, ...], np.float64]
        geographic coordinate array (1: longitude, 2: latitude)
    """

def get_domain_limits_unstr(npoints_p: int, coords_p: np.ndarray) -> None:
    r"""Set the domain limits for unstructured domain decomposed local domain.

    This is for the use of :func:`pyPDAF.PDAFomi.set_use_global_obs`.
    Currently, it only supports 2D limitations.

    See also `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_get_domain_limits_unstrc>`_

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

def store_obs_l_index(i_obs: int, idx: int, id_obs_l: int, distance: float,
                      cradius_l: float, sradius_l: float) -> None:
    r"""Save local observation information in PDAF.

    One should check `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_
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

def store_obs_l_index_vdist(i_obs: int, idx: int, id_obs_l: int, distance: float,
                            cradius_l: float, sradius_l: float, vdist: float) -> None:
    r"""Save local observation information for 2+1D factorized localization
    in the vertical direction in PDAF.

    One should check `relevant PDAF wiki page
    <https://pdaf.awi.de/trac/wiki/OMI_search_local_observations>`_
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
