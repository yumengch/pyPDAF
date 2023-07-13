import pyPDAF.UserFunc as PDAFcython
cimport pyPDAF.UserFunc as c__PDAFcython

import numpy as np
import sys
from traceback import print_exception
import mpi4py.MPI as MPI
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


def init (int n_obs
         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    n_obs : int
        number of observations
    """

    c__init (&n_obs
            )

def set_doassim (int i_obs,
                 int doassim
                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

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

def set_disttype (int i_obs,
                  int disttype
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

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

def set_ncoord (int i_obs,
                int ncoord
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

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

def set_id_obs_p (int i_obs,
                  id_obs_p
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    id_obs_p : ndarray[int]
        setter value
    """
    cdef int[::1] id_obs_p_view = np.array(id_obs_p, dtype=np.intc).ravel(order='F')
    cdef int nrows, dim_obs_p
    nrows, dim_obs_p,  = id_obs_p.shape


    c__pdafomi_set_id_obs_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &id_obs_p_view[0]
                            )

def set_icoeff_p (int i_obs,
                  icoeff_p
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    icoeff_p : ndarray[float]
        setter value
    """
    cdef double[::1] icoeff_p_view = np.array(icoeff_p).ravel(order='F')
    cdef int nrows, dim_obs_p
    nrows, dim_obs_p,  = icoeff_p.shape


    c__pdafomi_set_icoeff_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &icoeff_p_view[0]
                            )

def set_domainsize (int i_obs,
                    domainsize
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    domainsize : ndarray[float]
        setter value
    """
    cdef double[::1] domainsize_view = np.array(domainsize).ravel(order='F')
    cdef int ncoord
    ncoord,  = domainsize.shape


    c__pdafomi_set_domainsize (&i_obs,
                               &ncoord,
                               &domainsize_view[0]
                              )

def set_obs_err_type (int i_obs,
                      int obs_err_type
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

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

def set_use_global_obs (int i_obs,
                        int use_global_obs
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

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

def gather_obs (int i_obs,
                obs_p,
                ivar_obs_p,
                ocoord_p,
                double local_range
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_p : ndarray[float]
        pe-local observation vector
    ivar_obs_p : ndarray[float]
        pe-local inverse observation error variance
    ocoord_p : ndarray[float]
        pe-local observation coordinates
    local_range : float
        localization radius

    Returns
    -------
    dim_obs : int
        full number of observations
    """
    cdef double[::1] obs_p_view = np.array(obs_p).ravel(order='F')
    cdef double[::1] ivar_obs_p_view = np.array(ivar_obs_p).ravel(order='F')
    cdef double[::1] ocoord_p_view = np.array(ocoord_p).ravel(order='F')
    cdef int dim_obs_p
    _, dim_obs_p,  = ocoord_p.shape


    cdef int dim_obs

    c__pdafomi_gather_obs (&i_obs,
                           &dim_obs_p,
                           &obs_p_view[0],
                           &ivar_obs_p_view[0],
                           &ocoord_p_view[0],
                           &local_range,
                           &dim_obs
                          )

    return dim_obs

def gather_obsstate (int i_obs,
                     obsstate_p,
                     obsstate_f
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    obsstate_p : ndarray[float]
        vector of process-local observed state
    obsstate_f : ndarray[float]
        full observed vector for all types

    Returns
    -------
    obsstate_f : ndarray[float]
        full observed vector for all types
    """
    cdef double[::1] obsstate_p_view = np.array(obsstate_p).ravel(order='F')
    cdef double[::1] obsstate_f_view = np.array(obsstate_f).ravel(order='F')
    cdef int nobs_f_all
    nobs_f_all,  = obsstate_f.shape


    c__pdafomi_gather_obsstate (&i_obs,
                                &obsstate_p_view[0],
                                &obsstate_f_view[0],
                                &nobs_f_all
                               )

    return np.asarray(obsstate_f_view).reshape((nobs_f_all), order='F')

def localize_covar (int i_obs,
                    int locweight,
                    double local_range,
                    double srange,
                    coords_p,
                    hp_p,
                    hph
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    locweight : int
        localization weight type
    local_range : float
        localization radius
    srange : float
        support radius for weight functions
    coords_p : ndarray[float]
        coordinates of state vector elements
    hp_p : ndarray[float]
        matrix hp, dimension (nobs, dim)
    hph : ndarray[float]
        matrix hph, dimension (nobs, nobs)

    Returns
    -------
    hp_p : ndarray[float]
        matrix hp, dimension (nobs, dim)
    hph : ndarray[float]
        matrix hph, dimension (nobs, nobs)
    """
    cdef double[::1] coords_p_view = np.array(coords_p).ravel(order='F')
    cdef double[::1] hp_p_view = np.array(hp_p).ravel(order='F')
    cdef double[::1] hph_view = np.array(hph).ravel(order='F')
    cdef int dim_obs, dim_p, dim_coords
    dim_coords, dim_p,  = coords_p.shape
    dim_obs, _,  = hp_p.shape


    c__pdafomi_localize_covar (&i_obs,
                               &dim_p,
                               &dim_obs,
                               &dim_coords,
                               &locweight,
                               &local_range,
                               &srange,
                               &coords_p_view[0],
                               &hp_p_view[0],
                               &hph_view[0]
                              )

    return np.asarray(hp_p_view).reshape((dim_obs,dim_p), order='F'), np.asarray(hph_view).reshape((dim_obs,dim_obs), order='F')

def set_domain_limits (lim_coords
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    lim_coords : ndarray[float]
        geographic coordinate array (1: longitude, 2: latitude)
    """
    cdef double[::1] lim_coords_view = np.array(lim_coords).ravel(order='F')

    c__pdafomi_set_domain_limits (&lim_coords_view[0]
                                 )

def set_debug_flag (int debugval
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    debugval : int
        value for debugging flag
    """

    c__pdafomi_set_debug_flag (&debugval
                              )

def deallocate_obs (int i_obs
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    """

    c__pdafomi_deallocate_obs (&i_obs
                              )

def init_dim_obs_l (int i_obs,
                    coords_l,
                    int locweight,
                    double local_range,
                    double srange,
                    int dim_obs_l
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    coords_l : ndarray[float]
        coordinates of current local analysis domain
    locweight : int
        type of localization function
    local_range : float
        localization radius
    srange : float
        support radius of localization function
    dim_obs_l : int
        local dimension of current observation vector

    Returns
    -------
    dim_obs_l : int
        local dimension of current observation vector
    """
    cdef double[::1] coords_l_view = np.array(coords_l).ravel(order='F')

    c__pdafomi_init_dim_obs_l (&i_obs,
                               &coords_l_view[0],
                               &locweight,
                               &local_range,
                               &srange,
                               &dim_obs_l
                              )

    return dim_obs_l

def obs_op_gridpoint (int i_obs,
                      state_p,
                      obs_f_all
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_gridpoint (&i_obs,
                                 &state_p_view[0],
                                 &dim_p,
                                 &obs_f_all_view[0],
                                 &nobs_f_all
                                )

    return np.asarray(obs_f_all_view).reshape((nobs_f_all), order='F')

def obs_op_gridavg (int i_obs,
                    int nrows,
                    state_p,
                    obs_f_all
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_gridavg (&i_obs,
                               &nrows,
                               &state_p_view[0],
                               &dim_p,
                               &obs_f_all_view[0],
                               &nobs_f_all
                              )

    return np.asarray(obs_f_all_view).reshape((nobs_f_all), order='F')

def obs_op_interp_lin (int i_obs,
                       int nrows,
                       state_p,
                       obs_f_all
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_interp_lin (&i_obs,
                                  &nrows,
                                  &state_p_view[0],
                                  &dim_p,
                                  &obs_f_all_view[0],
                                  &nobs_f_all
                                 )

    return np.asarray(obs_f_all_view).reshape((nobs_f_all), order='F')

def obs_op_adj_gridavg (int i_obs,
                        int nrows,
                        state_p,
                        obs_f_all
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_adj_gridavg (&i_obs,
                                   &nrows,
                                   &state_p_view[0],
                                   &dim_p,
                                   &obs_f_all_view[0],
                                   &nobs_f_all
                                  )

    return np.asarray(state_p_view).reshape((dim_p), order='F')

def obs_op_adj_gridpoint (int i_obs,
                          state_p,
                          obs_f_all
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_adj_gridpoint (&i_obs,
                                     &state_p_view[0],
                                     &dim_p,
                                     &obs_f_all_view[0],
                                     &nobs_f_all
                                    )

    return np.asarray(state_p_view).reshape((dim_p), order='F')

def obs_op_adj_interp_lin (int i_obs,
                           int nrows,
                           state_p,
                           obs_f_all
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p)
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all)

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p)
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef double[::1] obs_f_all_view = np.array(obs_f_all).ravel(order='F')
    cdef int nobs_f_all, dim_p
    dim_p,  = state_p.shape
    nobs_f_all,  = obs_f_all.shape


    c__pdafomi_obs_op_adj_interp_lin (&i_obs,
                                      &nrows,
                                      &state_p_view[0],
                                      &dim_p,
                                      &obs_f_all_view[0],
                                      &nobs_f_all
                                     )

    return np.asarray(state_p_view).reshape((dim_p), order='F')

def get_interp_coeff_tri (gpc,
                          oc,
                          icoeff
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points; dim(3,2)
    oc : ndarray[float]
        3 rows; each containing lon and lat coordinates
         coordinates of observation; dim(2)
    icoeff : ndarray[float]
        interpolation coefficients; dim(3)

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients; dim(3)
    """
    cdef double[::1] gpc_view = np.array(gpc).ravel(order='F')
    cdef double[::1] oc_view = np.array(oc).ravel(order='F')
    cdef double[::1] icoeff_view = np.array(icoeff).ravel(order='F')

    c__pdafomi_get_interp_coeff_tri (&gpc_view[0],
                                     &oc_view[0],
                                     &icoeff_view[0]
                                    )

    return np.asarray(icoeff_view).reshape((3), order='F')

def get_interp_coeff_lin1d (gpc,
                            double oc,
                            icoeff
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points (dim=2)
    oc : float
        coordinates of observation
    icoeff : ndarray[float]
        interpolation coefficients (dim=2)

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients (dim=2)
    """
    cdef double[::1] gpc_view = np.array(gpc).ravel(order='F')
    cdef double[::1] icoeff_view = np.array(icoeff).ravel(order='F')

    c__pdafomi_get_interp_coeff_lin1d (&gpc_view[0],
                                       &oc,
                                       &icoeff_view[0]
                                      )

    return np.asarray(icoeff_view).reshape((2), order='F')

def get_interp_coeff_lin (gpc,
                          oc,
                          icoeff
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points
    oc : ndarray[float]
        coordinates of observation
    icoeff : ndarray[float]
        interpolation coefficients (num_gp)

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients (num_gp)
    """
    cdef double[::1] gpc_view = np.array(gpc).ravel(order='F')
    cdef double[::1] oc_view = np.array(oc).ravel(order='F')
    cdef double[::1] icoeff_view = np.array(icoeff).ravel(order='F')
    cdef int num_gp, n_dim
    num_gp, n_dim,  = gpc.shape


    c__pdafomi_get_interp_coeff_lin (&num_gp,
                                     &n_dim,
                                     &gpc_view[0],
                                     &oc_view[0],
                                     &icoeff_view[0]
                                    )

    return np.asarray(icoeff_view).reshape((num_gp), order='F')

def assimilate_3dvar (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step, time and dimension of next observation
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_en3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step, time and dimension of next observation
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_en3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of full observation vector
    py__obs_op_f_pdaf : func
        full observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize local dimimension of obs. vector
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from local state
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step, time and dimension of next observation
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_global (py__collect_state_pdaf,
                       py__distribute_state_pdaf,
                       py__init_dim_obs_pdaf,
                       py__obs_op_pdaf,
                       py__prepoststep_pdaf,
                       py__next_observation_pdaf
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step and time of next observation

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_hyb3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_ens_pdaf : func
        apply ensemble control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint ensemble control vector transform matrix
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step, time and dimension of next observation
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_hyb3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of full observation vector
    py__obs_op_f_pdaf : func
        full observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize local dimimension of obs. vector
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from local state
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step, time and dimension of next observation
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_lenkf (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__prepoststep_pdaf,
                      py__localize_covar_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__localize_covar_pdaf : func
        apply localization to hp and hph^t
    py__next_observation_pdaf : func
        provide time step and time of next observation

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__localize_covar_pdaf = py__localize_covar_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def assimilate_local (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize dim. of obs. vector for local ana. domain
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from state on local analysis domain
    py__next_observation_pdaf : func
        provide time step and time of next observation

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def generate_obs (py__collect_state_pdaf,
                  py__distribute_state_pdaf,
                  py__init_dim_obs_f_pdaf,
                  py__obs_op_f_pdaf,
                  py__get_obs_f_pdaf,
                  py__prepoststep_pdaf,
                  py__next_observation_pdaf
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__get_obs_f_pdaf : func
        provide observation vector to user
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        provide time step and time of next observation

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__get_obs_f_pdaf = py__get_obs_f_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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

def put_state_3dvar (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__cvt_pdaf,
                     py__cvt_adj_pdaf,
                     py__obs_op_lin_pdaf,
                     py__obs_op_adj_pdaf,
                     py__prepoststep_pdaf,
                     int outflag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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

def put_state_en3dvar_estkf (py__collect_state_pdaf,
                             py__init_dim_obs_pdaf,
                             py__obs_op_pdaf,
                             py__cvt_ens_pdaf,
                             py__cvt_adj_ens_pdaf,
                             py__obs_op_lin_pdaf,
                             py__obs_op_adj_pdaf,
                             py__prepoststep_pdaf,
                             int outflag
                            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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

def put_state_en3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of full observation vector
    py__obs_op_f_pdaf : func
        full observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize local dimimension of obs. vector
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from local state
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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

def put_state_generate_obs (py__collect_state_pdaf,
                            py__init_dim_obs_f_pdaf,
                            py__obs_op_f_pdaf,
                            py__get_obs_f_pdaf,
                            py__prepoststep_pdaf
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__get_obs_f_pdaf : func
        provide observation vector to user
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__get_obs_f_pdaf = py__get_obs_f_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__get_obs_f_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       &flag
                                      )

    return flag

def put_state_global (py__collect_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__prepoststep_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_global (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 &flag
                                )

    return flag

def put_state_hyb3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__cvt_ens_pdaf : func
        apply ensemble control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint ensemble control vector transform matrix
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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

def put_state_hyb3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_f_pdaf : func
        initialize dimension of full observation vector
    py__obs_op_f_pdaf : func
        full observation operator
    py__cvt_ens_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize local dimimension of obs. vector
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from local state
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    outflag : int
        status flag

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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

def put_state_lenkf (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__prepoststep_pdaf,
                     py__localize_covar_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__localize_covar_pdaf : func
        apply localization to hp and hph^t

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__localize_covar_pdaf = py__localize_covar_pdaf

    cdef int flag

    c__pdafomi_put_state_lenkf (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__init_dim_obs_pdaf,
                                c__PDAFcython.c__obs_op_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__localize_covar_pdaf,
                                &flag
                               )

    return flag

def put_state_local (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__prepoststep_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__init_n_domains_p_pdaf : func
        provide number of local analysis domains
    py__init_dim_l_pdaf : func
        init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : func
        initialize dim. of obs. vector for local ana. domain
    py__g2l_state_pdaf : func
        get state on local ana. domain from full state
    py__l2g_state_pdaf : func
        init full state from state on local analysis domain

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf

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

def init_obs_f_cb (int step,
                   int dim_obs_f,
                   observation_f
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        dimension of full observation vector
    observation_f : ndarray[float]
        full observation vector

    Returns
    -------
    observation_f : ndarray[float]
        full observation vector
    """
    cdef double[::1] observation_f_view = np.array(observation_f).ravel(order='F')
    c__pdafomi_init_obs_f_cb (&step,
                              &dim_obs_f,
                              &observation_f_view[0]
                             )

    return np.asarray(observation_f_view).reshape((dim_obs_f), order='F')

def init_obsvar_cb (int step,
                    int dim_obs_p,
                    obs_p,
                    double meanvar
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        pe-local dimension of observation vector
    obs_p : ndarray[float]
        pe-local observation vector
    meanvar : float
        mean observation error variance

    Returns
    -------
    meanvar : float
        mean observation error variance
    """
    cdef double[::1] obs_p_view = np.array(obs_p).ravel(order='F')
    c__pdafomi_init_obsvar_cb (&step,
                               &dim_obs_p,
                               &obs_p_view[0],
                               &meanvar
                              )

    return meanvar

def g2l_obs_cb (int domain_p,
                int step,
                int dim_obs_f,
                int dim_obs_l,
                ostate_f,
                ostate_l
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_f : int
        dimension of full pe-local observation vector
    dim_obs_l : int
        dimension of local observation vector
    ostate_f : ndarray[float]
        full pe-local obs.ervation vector
    ostate_l : ndarray[float]
        observation vector on local domain

    Returns
    -------
    ostate_l : ndarray[float]
        observation vector on local domain
    """
    cdef double[::1] ostate_f_view = np.array(ostate_f).ravel(order='F')
    cdef double[::1] ostate_l_view = np.array(ostate_l).ravel(order='F')
    c__pdafomi_g2l_obs_cb (&domain_p,
                           &step,
                           &dim_obs_f,
                           &dim_obs_l,
                           &ostate_f_view[0],
                           &ostate_l_view[0]
                          )

    return np.asarray(ostate_l_view).reshape((dim_obs_l), order='F')

def init_obs_l_cb (int domain_p,
                   int step,
                   int dim_obs_l,
                   observation_l
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain index
    step : int
        current time step
    dim_obs_l : int
        local dimension of observation vector
    observation_l : ndarray[float]
        local observation vector

    Returns
    -------
    observation_l : ndarray[float]
        local observation vector
    """
    cdef double[::1] observation_l_view = np.array(observation_l).ravel(order='F')
    c__pdafomi_init_obs_l_cb (&domain_p,
                              &step,
                              &dim_obs_l,
                              &observation_l_view[0]
                             )

    return np.asarray(observation_l_view).reshape((dim_obs_l), order='F')

def init_obsvar_l_cb (int domain_p,
                      int step,
                      int dim_obs_l,
                      obs_l,
                      double meanvar_l
                     ):
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
    meanvar_l : float
        mean local observation error variance

    Returns
    -------
    meanvar_l : float
        mean local observation error variance
    """
    cdef double[::1] obs_l_view = np.array(obs_l).ravel(order='F')
    c__pdafomi_init_obsvar_l_cb (&domain_p,
                                 &step,
                                 &dim_obs_l,
                                 &obs_l_view[0],
                                 &meanvar_l
                                )

    return meanvar_l

def prodrinva_l_cb (int domain_p,
                    int step,
                    int dim_obs_l,
                    int rank,
                    obs_l,
                    a_l,
                    c_l
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        dimension of local observation vector
    rank : int
        rank of initial covariance matrix
    obs_l : ndarray[float]
        local vector of observations
    a_l : ndarray[float]
        input matrix
    c_l : ndarray[float]
        output matrix

    Returns
    -------
    a_l : ndarray[float]
        input matrix
    c_l : ndarray[float]
        output matrix
    """
    cdef double[::1] obs_l_view = np.array(obs_l).ravel(order='F')
    cdef double[::1] a_l_view = np.array(a_l).ravel(order='F')
    cdef double[::1] c_l_view = np.array(c_l).ravel(order='F')
    c__pdafomi_prodrinva_l_cb (&domain_p,
                               &step,
                               &dim_obs_l,
                               &rank,
                               &obs_l_view[0],
                               &a_l_view[0],
                               &c_l_view[0]
                              )

    return np.asarray(a_l_view).reshape((dim_obs_l,rank), order='F'), np.asarray(c_l_view).reshape((dim_obs_l,rank), order='F')

def likelihood_l_cb (int domain_p,
                     int step,
                     int dim_obs_l,
                     obs_l,
                     resid_l,
                     double lhood_l
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        pe-local dimension of obs. vector
    obs_l : ndarray[float]
        pe-local vector of observations
    resid_l : ndarray[float]
        input vector of residuum
    lhood_l : float
        output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[float]
        input vector of residuum
    lhood_l : float
        output vector - log likelihood
    """
    cdef double[::1] obs_l_view = np.array(obs_l).ravel(order='F')
    cdef double[::1] resid_l_view = np.array(resid_l).ravel(order='F')
    c__pdafomi_likelihood_l_cb (&domain_p,
                                &step,
                                &dim_obs_l,
                                &obs_l_view[0],
                                &resid_l_view[0],
                                &lhood_l
                               )

    return np.asarray(resid_l_view).reshape((dim_obs_l), order='F'), lhood_l

def prodrinva_cb (int step,
                  int dim_obs_p,
                  int ncol,
                  obs_p,
                  a_p,
                  c_p
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of pe-local observation vector
    ncol : int
        number of columns in a_p and c_p
    obs_p : ndarray[float]
        pe-local vector of observations
    a_p : ndarray[float]
        input matrix
    c_p : ndarray[float]
        output matrix

    Returns
    -------
    c_p : ndarray[float]
        output matrix
    """
    cdef double[::1] obs_p_view = np.array(obs_p).ravel(order='F')
    cdef double[::1] a_p_view = np.array(a_p).ravel(order='F')
    cdef double[::1] c_p_view = np.array(c_p).ravel(order='F')
    c__pdafomi_prodrinva_cb (&step,
                             &dim_obs_p,
                             &ncol,
                             &obs_p_view[0],
                             &a_p_view[0],
                             &c_p_view[0]
                            )

    return np.asarray(c_p_view).reshape((dim_obs_p,ncol), order='F')

def likelihood_cb (int step,
                   int dim_obs,
                   obs,
                   resid,
                   double lhood
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        pe-local dimension of obs. vector
    obs : ndarray[float]
        pe-local vector of observations
    resid : ndarray[float]
        input vector of residuum
    lhood : float
        output vector - log likelihood

    Returns
    -------
    lhood : float
        output vector - log likelihood
    """
    cdef double[::1] obs_view = np.array(obs).ravel(order='F')
    cdef double[::1] resid_view = np.array(resid).ravel(order='F')
    c__pdafomi_likelihood_cb (&step,
                              &dim_obs,
                              &obs_view[0],
                              &resid_view[0],
                              &lhood
                             )

    return lhood

def add_obs_error_cb (int step,
                      int dim_obs_p,
                      c_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of pe-local observation vector
    c_p : ndarray[float]
        matrix to which r is added

    Returns
    -------
    c_p : ndarray[float]
        matrix to which r is added
    """
    cdef double[::1] c_p_view = np.array(c_p).ravel(order='F')
    c__pdafomi_add_obs_error_cb (&step,
                                 &dim_obs_p,
                                 &c_p_view[0]
                                )

    return np.asarray(c_p_view).reshape((dim_obs_p,dim_obs_p), order='F')

def init_obscovar_cb (int step,
                      int dim_obs,
                      int dim_obs_p,
                      covar,
                      m_state_p,
                      bint isdiag
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    dim_obs_p : int
        pe-local dimension of obs. vector
    covar : ndarray[float]
        observation error covar. matrix
    m_state_p : ndarray[float]
        observation vector
    isdiag : bool
        whether matrix r is diagonal

    Returns
    -------
    covar : ndarray[float]
        observation error covar. matrix
    isdiag : bool
        whether matrix r is diagonal
    """
    cdef double[::1] covar_view = np.array(covar).ravel(order='F')
    cdef double[::1] m_state_p_view = np.array(m_state_p).ravel(order='F')
    c__pdafomi_init_obscovar_cb (&step,
                                 &dim_obs,
                                 &dim_obs_p,
                                 &covar_view[0],
                                 &m_state_p_view[0],
                                 &isdiag
                                )

    return np.asarray(covar_view).reshape((dim_obs,dim_obs), order='F'), isdiag

def init_obserr_f_cb (int step,
                      int dim_obs_f,
                      obs_f,
                      obserr_f
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        full dimension of observation vector
    obs_f : ndarray[float]
        full observation vector
    obserr_f : ndarray[float]
        full observation error stddev

    Returns
    -------
    obserr_f : ndarray[float]
        full observation error stddev
    """
    cdef double[::1] obs_f_view = np.array(obs_f).ravel(order='F')
    cdef double[::1] obserr_f_view = np.array(obserr_f).ravel(order='F')
    c__pdafomi_init_obserr_f_cb (&step,
                                 &dim_obs_f,
                                 &obs_f_view[0],
                                 &obserr_f_view[0]
                                )

    return np.asarray(obserr_f_view).reshape((dim_obs_f), order='F')

def prodrinva_hyb_l_cb (int domain_p,
                        int step,
                        int dim_obs_l,
                        int rank,
                        obs_l,
                        double alpha,
                        a_l,
                        c_l
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        dimension of local observation vector
    rank : int
        rank of initial covariance matrix
    obs_l : ndarray[float]
        local vector of observations
    alpha : float
        hybrid weight
    a_l : ndarray[float]
        input matrix
    c_l : ndarray[float]
        output matrix

    Returns
    -------
    a_l : ndarray[float]
        input matrix
    c_l : ndarray[float]
        output matrix
    """
    cdef double[::1] obs_l_view = np.array(obs_l).ravel(order='F')
    cdef double[::1] a_l_view = np.array(a_l).ravel(order='F')
    cdef double[::1] c_l_view = np.array(c_l).ravel(order='F')
    c__pdafomi_prodrinva_hyb_l_cb (&domain_p,
                                   &step,
                                   &dim_obs_l,
                                   &rank,
                                   &obs_l_view[0],
                                   &alpha,
                                   &a_l_view[0],
                                   &c_l_view[0]
                                  )

    return np.asarray(a_l_view).reshape((dim_obs_l,rank), order='F'), np.asarray(c_l_view).reshape((dim_obs_l,rank), order='F')

def likelihood_hyb_l_cb (int domain_p,
                         int step,
                         int dim_obs_l,
                         obs_l,
                         resid_l,
                         double alpha,
                         double lhood_l
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    domain_p : int
        current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        pe-local dimension of obs. vector
    obs_l : ndarray[float]
        pe-local vector of observations
    resid_l : ndarray[float]
        input vector of residuum
    alpha : float
        hybrid weight
    lhood_l : float
        output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[float]
        input vector of residuum
    lhood_l : float
        output vector - log likelihood
    """
    cdef double[::1] obs_l_view = np.array(obs_l).ravel(order='F')
    cdef double[::1] resid_l_view = np.array(resid_l).ravel(order='F')
    c__pdafomi_likelihood_hyb_l_cb (&domain_p,
                                    &step,
                                    &dim_obs_l,
                                    &obs_l_view[0],
                                    &resid_l_view[0],
                                    &alpha,
                                    &lhood_l
                                   )

    return np.asarray(resid_l_view).reshape((dim_obs_l), order='F'), lhood_l

