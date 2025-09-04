# pylint: disable=unused-argument
"""Stub file for PDAFomi setter module
"""
import numpy as np

def set_doassim(i_obs: int, doassim: int) -> None:
    r"""Setting the `doassim` attribute of `obs_f`
    for `i`-th observation type. This property must be
    explicitly set for OMI functionality.

    Properties of `obs_f` are typically set in user-supplied function
    `py__init_dim_obs_pdaf`.

    This is by default set to 0, which means that
    the given type of observation is not assimilated in the DA system.

    Parameters
    ----------
    i_obs : int
        index of observation types
    doassim : int
        0) do not assimilate;
        1) assimilate the observation type
    """

def set_disttype(i_obs: int, disttype: int) -> None:
    r"""Setting the observation localisation distance
    calculation method
    for `i`-th observation type. This is a mandatory property
    for OMI functionality.

    Properties of `obs_f` are typically set in user-supplied function
    `py__init_dim_obs_pdaf`.

    `disttype` determines the way the distance
    between observation and model grid is calculated in OMI.
    To perform distance computation, the observation coordinates should be
    given by `ocoord_p` argument
    when :func:`pyPDAF.PDAF.omi_gather_obs` is called.

    See also `PDAF distance computation
    <https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype>`_.

    Parameters
    ----------
    i_obs : int
        index of observations
    disttype : int
        Type of distance used for localisation
            - 0) Cartesian (any units)
            - 1) Cartesian periodic (any units)
            - 2) Approximation to geographic distance in metres using
              latitude and longitude expressed in radians
            - 3) Using Haversine formula to compute distance in metres
              between two points on the surface of a sphere
            - 10) 3D Cartesian distance where horizontal and vertical
              distances are treated separately
            - 11) 3D Cartesian periodic distance where horizontal and
              vertical distances are treated separately
            - 12) Same as 2) for horizontal distance but vertical
              distance is in units chosen by users where the horizontal
              and vertical distances are treated separately
            - 13) Same as 3) for horizontal distance but vertical
              distance is in units chosen by users where the horizontal
              and vertical distances are treated separately
    """

def set_ncoord(i_obs: int, ncoord: int) -> None:
    r"""Setting the number of spatial dimensions of observations
    for `i`-th observation type. This is a mandatory property
    for OMI functionality.

    Properties of `obs_f` are typically set in user-supplied function
    `py__init_dim_obs_pdaf`.

    `ncoord` gives the coordinate dimension of the observation.
    This information is used by observation distance computation
    for localisation.
    For example, `ncoord=2` for 2D observation coordinates.

    Parameters
    ----------
    i_obs : int
        index of observations
    ncoord : int
        Dimension of the observation coordinate
    """

def set_obs_err_type(i_obs: int, obs_err_type: int) -> None:
    r"""Setting the type of observation error distribution
    for `i`-th observation type. This property is optional
    unless a laplacian observation error distribution is used.

    The function is typically used in user-supplied
    function `py__init_dim_obs_pdaf`.

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_err_type : int
        type of observation error distribution
            0) Gaussian (default)
            1) double exponential (Laplacian)
    """

def set_use_global_obs(i_obs: int, use_global_obs: int) -> None:
    r"""Switch for only use process-local observations
    for `i`-th observation type.

    The function is typically used in user-supplied
    function `py__init_dim_obs_pdaf`.

    The filters can be performed in parallel
    based on the filtering communicator, `comm_filter`.
    This is typically the case for the domain-localised filters,
    e.g., LESTK, LETKF, LSEIK, LNETF.
    In this case, observation vectors are stored in
    process-local vectors, `obs_p`. Each local
    process (`obs_p`) only stores a section of the full
    observation vector. This typically corresponds to the
    local domain corresponding to the filtering process,
    based on model domain decomposition.

    By default, `use_global_obs=1`. This means that
    PDAF-OMI gathers the entire observation vector for all processes.
    One can choose to only use process-local observations
    for global filters, or within localisation radius for by setting `use_global_obs=0`.
    This can save computational cost used for
    observation distance calculations.

    However, it needs additional preparations to make
    PDAF-OMI aware of the limiting coordinates
    of a process sub-domain using
    :func:`pyPDAF.PDAFomi.set_domain_limits` or
    :func:`pyPDAF.PDAFomi.set_domain_limits_unstruc`.


    See Also
    --------
    https://pdaf.awi.de/trac/wiki/OMI_use_global_obs

    Parameters
    ----------
    i_obs : int
        index of observations
    use_global_obs : int
        Swith to use global observations or not
            0) Using process-local observations;
            1) using cross-process observations (default)
    """

def set_inno_omit(i_obs: int, inno_omit: float) -> None:
    r"""Setting innovation threshold for removing observation
    outliers. By default, no observations are omitted.

    This function is typically used in user-supplied
    function :func:`py__init_dim_obs_pdaf`.

    The observation omission is only activated when it is > 0.0.
    PDAF will omit observations where their squared
    the innovation of the ensemble mean is larger than
    the product of `inno_omit` and observation error variance.

    The observations are omitted by setting a very large
    observation error variance, i.e., a very small
    inverse of the observation error variance,  `inno_omit_ivar`.
    This can be set by :func:`pyPDAF.PDAF.omi_set_inno_omit_ivar`.

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit : float
        Threshold of innovation to be omitted
    """

def set_inno_omit_ivar(i_obs: int, inno_omit_ivar: float) -> None:
    r"""Setting the inverse of observation error variance for
    omitted observations.

    This should be set to a very small value relative to
    assimilated observations. By default, it is set to `1e-12`.

    This function is typically used in user-supplied function
    :func:`py__init_dim_obs_pdaf`.

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit_ivar : float
        Inverse of observation variance for omiited observations
    """

def set_id_obs_p(i_obs: int, nrows: int, dim_obs_p: int, id_obs_p: np.ndarray) -> None:
    r"""Setting the `id_obs_p` attribute of `obs_f`
    for `i`-th observation type. This is a mandatory property
    for OMI functionality.

    The function is typically used in user-supplied
    function `py__init_dim_obs_pdaf`.

    Here, `id_obs_p(nrows, dim_obs_p)` is a 2D array of integers.
    The value of `nrows` depends on the observation operator
    used for an observation.

    Examples:

    - `nrows=1`: observations are located on model grid point.
      In this case,
      `id_obs_p` stores the index of the state vector
      (starting from 1) corresponds to the observations,
      e.g. `id_obs_p[0, j] = i` means that the location
      and variable of the `i`-th element of the state vector
      is the same as the `j`-th observation.

    - `nrows=4`: each observation corresponds to
      4 indices of elements in the state vector.
      In this case,
      the location of these elements is used to perform bi-linear interpolation
      from model grid to observation location.
      For interpolation, this information is used in the
      :func:`pyPDAF.PDAF.omi_obs_op_interp_lin` functions.
      This information can also be used to
      perform a state vector averaging operator as
      observation operator in :func:`pyPDAF.PDAFomi.obs_op_gridavg`  When interpolation is needed,
      the weighting of the interpolation is done
      in the :func:`pyPDAF.PDAFomi.get_interp_coeff_lin`,
      :func:`pyPDAF.PDAFomi.get_interp_coeff_lin1D`,
      and :func:`pyPDAF.PDAFomi.get_interp_coeff_tri` functions.
      The details of interpolation setup can be found at
      `PDAF wiki page
      <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_.


    Parameters
    ----------
    i_obs : int
        index of observations
    id_obs_p : ndarray[tuple[nrows, dim_obs_p, ...], np.intc]
        indice corresponds to observations in the state vector
        The 1st-th dimension nrows is number of values to be averaged or used for interpolation
        The 2nd-th dimension dim_obs_p is dimension of PE local obs
    """

def set_icoeff_p(i_obs: int, nrows: int, dim_obs_p: int, icoeff_p: np.ndarray) -> None:
    r"""Setting the observation interpolation coefficient
    for `i`-th observation type. This property is optional
    unless interpolations needed in observation operators
    operator.

    The function is typically used in user-supplied
    function `py__init_dim_obs_pdaf`.

    `icoeff_p(nrows, dim_obs_p)` is a 2D array of real number
    used to interpolate state vector to point-wise observation grid.
    The `nrows` is the number of state vector used to interpolate
    to one observation location.

    A suite of functions are provided to obtain these coefficients,
    which depend on `obs_f` attribute of `id_obs_p` and
    observation coordinates.

    See also :func:`pyPDAF.PDAF.set_id_obs_p`:
        - :func:`pyPDAF.PDAFomi.get_interp_coeff_lin1D`
          1D interpolation coefficient
        - :func:`pyPDAF.PDAFomi.get_interp_coeff_lin`
          linear interpolation coefficient for 1, 2 and 3D rectangular grids
        - :func:`pyPDAF.PDAFomi.get_interp_coeff_tri`
          2D linear interpolation for triangular grids

    See also `PDAF documentation for OMI interpolations
    <https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_.

    Parameters
    ----------
    i_obs : int
        index of observations
    icoeff_p : ndarray[tuple[nrows, dim_obs_p, ...], np.float64]
        weighting coefficients for interpolations
        The 1st-th dimension nrows is number of state vector used to interpolate
        to one observation location
        The 2nd-th dimension dim_obs_p is dimension of PE local obs
    """

def set_domainsize(i_obs: int, ncoord: int, domainsize: np.ndarray) -> None:
    r"""Setting the domain periodicity
    attribute of `obs_f`
    for `i`-th observation type. This property is optional
    unless localisation is used.

    The function is typically used in user-supplied
    function `py__init_dim_obs_pdaf`.

    `domainsize(ncoord)` specifies the size of the domain
    in each spatial dimension.
    This information is used to compute the Cartesian disance
    with periodic boundary. That is `disttype = 1 or 11`
    Domain size must be positive.
    If the value of one dimension is `<=0`,
    no periodicity is assumed in that dimension.

    Parameters
    ----------
    i_obs : int
        index of observations
    ncoord: int
        Number of spatial dimensions
    domainsize : ndarray[tuple[ncoord, ...], np.float64]
        Size of the domain in each dimension
        The array dimension `ncoord` is state dimension
    """

def set_name(i_obs: int, obsname: str) -> None:
    """Set a name for given observation type

    Parameters
    ----------
    i_obs : int
        index of observation type
    obsname : str
        name of observation type
    """
