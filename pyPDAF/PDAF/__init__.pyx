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


def assimilate_3dvar (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_pdaf,
                      py__obs_op_pdaf,
                      py__init_obs_pdaf,
                      py__prodrinva_pdaf,
                      py__cvt_pdaf,
                      py__cvt_adj_pdaf,
                      py__obs_op_lin_pdaf,
                      py__obs_op_adj_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
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
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_3dvar (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__prodrinva_pdaf,
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
                              py__prodrinva_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_obsvar_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__distribute_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodrinva_pdaf,
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
                               py__prodrinva_pdaf,
                               py__cvt_ens_pdaf,
                               py__cvt_adj_ens_pdaf,
                               py__obs_op_lin_pdaf,
                               py__obs_op_adj_pdaf,
                               py__init_dim_obs_f_pdaf,
                               py__obs_op_f_pdaf,
                               py__init_obs_f_pdaf,
                               py__init_obs_l_pdaf,
                               py__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__init_obs_f_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__init_obs_f_pdaf = py__init_obs_f_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    c__pdaf_assimilate_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__distribute_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodrinva_pdaf,
                                       c__PDAFcython.c__cvt_ens_pdaf,
                                       c__PDAFcython.c__cvt_adj_ens_pdaf,
                                       c__PDAFcython.c__obs_op_lin_pdaf,
                                       c__PDAFcython.c__obs_op_adj_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__init_obs_f_pdaf,
                                       c__PDAFcython.c__init_obs_l_pdaf,
                                       c__PDAFcython.c__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__add_obs_err_pdaf : func
        add obs error covariance r to hph in enkf
    py__init_obs_covar_pdaf : func
        initialize obs. error cov. matrix r in enkf
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__add_obs_err_pdaf = py__add_obs_err_pdaf
    PDAFcython.py__init_obs_covar_pdaf = py__init_obs_covar_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
                      py__prodrinva_pdaf,
                      py__init_obsvar_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 hv
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_estkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodrinva_pdaf,
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
                     py__prodrinva_pdaf,
                     py__init_obsvar_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 hv
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_etkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_pdaf,
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
                               py__prodrinva_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__distribute_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodrinva_pdaf,
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
                                py__prodrinva_pdaf,
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
                                py__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__init_obs_f_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__init_obs_f_pdaf = py__init_obs_f_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int outflag

    c__pdaf_assimilate_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                        c__PDAFcython.c__distribute_state_pdaf,
                                        c__PDAFcython.c__init_dim_obs_pdaf,
                                        c__PDAFcython.c__obs_op_pdaf,
                                        c__PDAFcython.c__init_obs_pdaf,
                                        c__PDAFcython.c__prodrinva_pdaf,
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
                                        c__PDAFcython.c__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__localize_covar_pdaf : func
        apply localization to hp and hph^t
    py__add_obs_err_pdaf : func
        add obs error covariance r to hph in enkf
    py__init_obs_covar_pdaf : func
        initialize obs. error cov. matrix r in enkf
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__localize_covar_pdaf = py__localize_covar_pdaf
    PDAFcython.py__add_obs_err_pdaf = py__add_obs_err_pdaf
    PDAFcython.py__init_obs_covar_pdaf = py__init_obs_covar_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
                       py__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lestkf (c__PDAFcython.c__collect_state_pdaf,
                               c__PDAFcython.c__distribute_state_pdaf,
                               c__PDAFcython.c__init_dim_obs_pdaf,
                               c__PDAFcython.c__obs_op_pdaf,
                               c__PDAFcython.c__init_obs_pdaf,
                               c__PDAFcython.c__init_obs_l_pdaf,
                               c__PDAFcython.c__prepoststep_pdaf,
                               c__PDAFcython.c__prodrinva_l_pdaf,
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
                      py__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_letkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodrinva_l_pdaf,
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
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_l_pdaf : func
        compute observation likelihood for an ensemble member
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
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
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_l_pdaf = py__likelihood_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
                       py__prodrinva_l_pdaf,
                       py__prodrinva_hyb_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
    py__prodrinva_hyb_l_pdaf : func
        provide product r^-1 a on local analysis domain with hybrid weight
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__likelihood_l_pdaf : func
        compute likelihood
    py__likelihood_hyb_l_pdaf : func
        compute likelihood with hybrid weight
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__prodrinva_hyb_l_pdaf = py__prodrinva_hyb_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__likelihood_l_pdaf = py__likelihood_l_pdaf
    PDAFcython.py__likelihood_hyb_l_pdaf = py__likelihood_hyb_l_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lknetf (c__PDAFcython.c__collect_state_pdaf,
                               c__PDAFcython.c__distribute_state_pdaf,
                               c__PDAFcython.c__init_dim_obs_pdaf,
                               c__PDAFcython.c__obs_op_pdaf,
                               c__PDAFcython.c__init_obs_pdaf,
                               c__PDAFcython.c__init_obs_l_pdaf,
                               c__PDAFcython.c__prepoststep_pdaf,
                               c__PDAFcython.c__prodrinva_l_pdaf,
                               c__PDAFcython.c__prodrinva_hyb_l_pdaf,
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
                      py__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_lseik (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__distribute_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodrinva_l_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_pdaf : func
        compute observation likelihood for an ensemble member
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_pdaf = py__likelihood_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_pdaf : func
        compute observation likelihood for an ensemble member
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_pdaf = py__likelihood_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
                     py__prodrinva_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 hv
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_seek (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_pdaf,
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
                     py__prodrinva_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 hv
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
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_seik (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__distribute_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_pdaf,
                             c__PDAFcython.c__next_observation_pdaf,
                             &flag
                            )

    return flag

def assimilate_prepost (py__collect_state_pdaf,
                        py__distribute_state_pdaf,
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
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_prepost (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__distribute_state_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__next_observation_pdaf,
                                &flag
                               )

    return flag

def deallocate ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    """
    c__pdaf_deallocate ()

def diag_crps (int element,
               oens,
               obs
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    element : int
        id of element to be used. if element=0, mean values over all elements are computed
    oens : ndarray[float]
        state ensemble
    obs : ndarray[float]
        state ensemble

    Returns
    -------
    crps : float
        crps
    reli : float
        reliability
    resol : float
        resolution
    uncert : float
        uncertainty
    status : int
        status flag (0=success)
    """
    cdef double[::1] oens_view = np.array(oens).ravel(order='F')
    cdef double[::1] obs_view = np.array(obs).ravel(order='F')
    cdef int dim, dim_ens
    dim, dim_ens,  = oens.shape


    cdef double crps
    cdef double reli
    cdef double resol
    cdef double uncert
    cdef int status

    c__pdaf_diag_crps (&dim,
                       &dim_ens,
                       &element,
                       &oens_view[0],
                       &obs_view[0],
                       &crps,
                       &reli,
                       &resol,
                       &uncert,
                       &status
                      )

    return crps, reli, resol, uncert, status

def diag_effsample (weights
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    weights : ndarray[float]
        weights of the samples

    Returns
    -------
    effsample : float
        effecfive sample size
    """
    cdef double[::1] weights_view = np.array(weights).ravel(order='F')
    cdef int dim_sample
    dim_sample,  = weights.shape


    cdef double effsample

    c__pdaf_diag_effsample (&dim_sample,
                            &weights_view[0],
                            &effsample
                           )

    return effsample

def diag_ensstats (int element,
                   state,
                   ens
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    element : int
        id of element to be used. if element=0, mean values over all elements are computed
    state : ndarray[float]
        state vector
    ens : ndarray[float]
        state ensemble

    Returns
    -------
    skewness : float
        skewness of ensemble
    kurtosis : float
        kurtosis of ensemble
    status : int
        status flag (0=success)
    """
    cdef double[::1] state_view = np.array(state).ravel(order='F')
    cdef double[::1] ens_view = np.array(ens).ravel(order='F')
    cdef int dim, dim_ens
    dim, dim_ens,  = ens.shape


    cdef double skewness
    cdef double kurtosis
    cdef int status

    c__pdaf_diag_ensstats (&dim,
                           &dim_ens,
                           &element,
                           &state_view[0],
                           &ens_view[0],
                           &skewness,
                           &kurtosis,
                           &status
                          )

    return skewness, kurtosis, status

def diag_histogram (int ncall,
                    int element,
                    state,
                    ens,
                    hist
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    ncall : int
        number of calls to routine
    element : int
        element of vector used for histogram
    state : ndarray[float]
        if element=0, all elements are used
         state vector
    ens : ndarray[float]
        state ensemble
    hist : ndarray[int]
        histogram about the state

    Returns
    -------
    hist : ndarray[int]
        histogram about the state
    delta : float
        deviation measure from flat histogram
    status : int
        status flag (0=success)
    """
    cdef double[::1] state_view = np.array(state).ravel(order='F')
    cdef double[::1] ens_view = np.array(ens).ravel(order='F')
    cdef int[::1] hist_view = np.array(hist, dtype=np.intc).ravel(order='F')
    cdef int dim, dim_ens
    dim, dim_ens,  = ens.shape


    cdef double delta
    cdef int status

    c__pdaf_diag_histogram (&ncall,
                            &dim,
                            &dim_ens,
                            &element,
                            &state_view[0],
                            &ens_view[0],
                            &hist_view[0],
                            &delta,
                            &status
                           )

    return np.asarray(hist_view).reshape((dim_ens+1), order='F'), delta, status

def eofcovar (dim_fields,
              offsets,
              int remove_mstate,
              int do_mv,
              states,
              meanstate,
              int verbose
             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_fields : ndarray[int]
        size of each field
    offsets : ndarray[int]
        start position of each field
    remove_mstate : int
        1: subtract mean state from states
         before computing eofs; 0: don't remove
    do_mv : int
        1: do multivariate scaling; 0: no scaling
         nfields, dim_fields and offsets are only used if do_mv=1
    states : ndarray[float]
        state perturbations
    meanstate : ndarray[float]
        mean state (only changed if remove_mstate=1)
    verbose : int
        verbosity flag

    Returns
    -------
    states : ndarray[float]
        state perturbations
    stddev : ndarray[float]
        standard deviation of field variability
         without multivariate scaling (do_mv=0), it is stddev = 1.0
    svals : ndarray[float]
        singular values divided by sqrt(nstates-1)
    svec : ndarray[float]
        singular vectors
    meanstate : ndarray[float]
        mean state (only changed if remove_mstate=1)
    status : int
        status flag
    """
    cdef int[::1] dim_fields_view = np.array(dim_fields, dtype=np.intc).ravel(order='F')
    cdef int[::1] offsets_view = np.array(offsets, dtype=np.intc).ravel(order='F')
    cdef double[::1] states_view = np.array(states).ravel(order='F')
    cdef double[::1] meanstate_view = np.array(meanstate).ravel(order='F')
    cdef int dim_state, nstates, nfields
    dim_state, nstates,  = states.shape
    nfields,  = dim_fields.shape


    cdef double [::1] stddev_view = np.zeros((nfields)).ravel()
    cdef double [::1] svals_view = np.zeros((nstates)).ravel()
    cdef double [::1] svec_view = np.zeros((dim_state,nstates)).ravel()
    cdef int status

    c__pdaf_eofcovar (&dim_state,
                      &nstates,
                      &nfields,
                      &dim_fields_view[0],
                      &offsets_view[0],
                      &remove_mstate,
                      &do_mv,
                      &states_view[0],
                      &stddev_view[0],
                      &svals_view[0],
                      &svec_view[0],
                      &meanstate_view[0],
                      &verbose,
                      &status
                     )

    return np.asarray(states_view).reshape((dim_state,nstates), order='F'), np.asarray(stddev_view).reshape((nfields), order='F'), np.asarray(svals_view).reshape((nstates), order='F'), np.asarray(svec_view).reshape((dim_state,nstates), order='F'), np.asarray(meanstate_view).reshape((dim_state), order='F'), status

def force_analysis ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    """
    c__pdaf_force_analysis ()

def gather_dim_obs_f (int dim_obs_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_obs_p : int
        pe-local observation dimension

    Returns
    -------
    dim_obs_f : int
        full observation dimension
    """

    cdef int dim_obs_f

    c__pdaf_gather_dim_obs_f (&dim_obs_p,
                              &dim_obs_f
                             )

    return dim_obs_f

def gather_obs_f (obs_p,
                  int dimobs_f
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    obs_p : ndarray[float]
        pe-local vector
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    obs_f : ndarray[float]
        full gathered vector
    status : int
        status flag:
         (0) no error
         (1) when pdaf_gather_dim_obs_f not executed before
    """
    cdef double[::1] obs_p_view = np.array(obs_p).ravel(order='F')
    cdef int dimobs_p
    dimobs_p,  = obs_p.shape


    cdef double [::1] obs_f_view = np.zeros((dimobs_f)).ravel()
    cdef int status

    c__pdaf_gather_obs_f (&obs_p_view[0],
                          &dimobs_p,
                          &obs_f_view[0],
                          &dimobs_f,
                          &status
                         )

    return np.asarray(obs_f_view).reshape((dimobs_f), order='F'), status

def gather_obs_f2 (coords_p,
                   int dimobs_f
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    coords_p : ndarray[float]
        pe-local array
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    coords_f : ndarray[float]
        full gathered array
    status : int
        status flag:
         (0) no error
         (1) when pdaf_gather dim_obs_f not executed before
    """
    cdef double[::1] coords_p_view = np.array(coords_p).ravel(order='F')
    cdef int dimobs_p, nrows
    nrows, dimobs_p,  = coords_p.shape


    cdef double [::1] coords_f_view = np.zeros((nrows,dimobs_f)).ravel()
    cdef int status

    c__pdaf_gather_obs_f2 (&coords_p_view[0],
                           &dimobs_p,
                           &coords_f_view[0],
                           &dimobs_f,
                           &nrows,
                           &status
                          )

    return np.asarray(coords_f_view).reshape((nrows,dimobs_f), order='F'), status

def gather_obs_f2_flex (int dim_obs_f,
                        coords_p
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_obs_f : int
        full observation dimension
    coords_p : ndarray[float]
        pe-local array

    Returns
    -------
    coords_f : ndarray[float]
        full gathered array
    status : int
        status flag: (0) no error
    """
    cdef double[::1] coords_p_view = np.array(coords_p).ravel(order='F')
    cdef int dim_obs_p, nrows
    nrows, dim_obs_p,  = coords_p.shape


    cdef double [::1] coords_f_view = np.zeros((nrows,dim_obs_f)).ravel()
    cdef int status

    c__pdaf_gather_obs_f2_flex (&dim_obs_p,
                                &dim_obs_f,
                                &coords_p_view[0],
                                &coords_f_view[0],
                                &nrows,
                                &status
                               )

    return np.asarray(coords_f_view).reshape((nrows,dim_obs_f), order='F'), status

def gather_obs_f_flex (int dim_obs_f,
                       obs_p
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    dim_obs_f : int
        full observation dimension
    obs_p : ndarray[float]
        pe-local vector

    Returns
    -------
    obs_f : ndarray[float]
        full gathered vector
    status : int
        status flag: (0) no error
    """
    cdef double[::1] obs_p_view = np.array(obs_p).ravel(order='F')
    cdef int dim_obs_p
    dim_obs_p,  = obs_p.shape


    cdef double [::1] obs_f_view = np.zeros((dim_obs_f)).ravel()
    cdef int status

    c__pdaf_gather_obs_f_flex (&dim_obs_p,
                               &dim_obs_f,
                               &obs_p_view[0],
                               &obs_f_view[0],
                               &status
                              )

    return np.asarray(obs_f_view).reshape((dim_obs_f), order='F'), status

def generate_obs (py__collect_state_pdaf,
                  py__distribute_state_pdaf,
                  py__init_dim_obs_f_pdaf,
                  py__obs_op_f_pdaf,
                  py__get_obs_f_pdaf,
                  py__init_obserr_f_pdaf,
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
    py__init_obserr_f_pdaf : func
        initialize vector of observation error standard deviations
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
    PDAFcython.py__init_obserr_f_pdaf = py__init_obserr_f_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 


    Returns
    -------
    did_assim : int
        flag: (1) for assimilation; (0) else
    """

    cdef int did_assim

    c__pdaf_get_assim_flag (&did_assim
                           )

    return did_assim

def get_ensstats ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 


    Returns
    -------
    dims : ndarray[int]
        dimension of pointer
    c_skew_ptr : ndarray[float]
        pointer to skewness array
    c_kurt_ptr : ndarray[float]
        pointer to kurtosis array
    status : int
        status flag
    """

    cdef int [::1] dims_view = np.zeros((1), dtype=np.intc).ravel()
    cdef double* c_skew_ptr
    cdef double* c_kurt_ptr
    cdef int status

    c__pdaf_get_ensstats (&dims_view[0],
                          &c_skew_ptr,
                          &c_kurt_ptr,
                          &status
                         )

    dims = np.asarray(dims_view)
    return np.asarray(<double[:np.prod(dims)]> c_skew_ptr).reshape(dims, order='F'), \
           np.asarray(<double[:np.prod(dims)]> c_kurt_ptr).reshape(dims, order='F'), \
           status

def get_localfilter ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 


    Returns
    -------
    lfilter : int
        whether the filter is domain-localized
    """

    cdef int lfilter

    c__pdaf_get_localfilter (&lfilter
                            )

    return lfilter

def get_memberid (int memberid
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    memberid : int
        index in the local ensemble

    Returns
    -------
    memberid : int
        index in the local ensemble
    """

    c__pdaf_get_memberid (&memberid
                         )

    return memberid

def get_obsmemberid (int memberid
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    memberid : int
        index in the local observed ensemble

    Returns
    -------
    memberid : int
        index in the local observed ensemble
    """

    c__pdaf_get_obsmemberid (&memberid
                            )

    return memberid

def get_smootherens ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 


    Returns
    -------
    c_sens_point : ndarray[float]
        pointer to smoother array
    maxlag : int
        number of past timesteps processed in sens
    dims : ndarray[int]
        dimension of pointer
    status : int
        status flag
    """

    cdef double* c_sens_point
    cdef int maxlag
    cdef int [::1] dims_view = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_get_smootherens (&c_sens_point,
                             &maxlag,
                             &dims_view[0],
                             &status
                            )

    dims = np.asarray(dims_view)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           maxlag, status

def get_state (int steps,
               int doexit,
               py__next_observation_pdaf,
               py__distribute_state_pdaf,
               py__prepoststep_pdaf,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    steps : int
        flag and number of time steps
    doexit : int
        whether to exit from forecasts
    py__next_observation_pdaf : func
        provide time step and time of next observation
    py__distribute_state_pdaf : func
        routine to distribute a state vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    flag : int
        status flag

    Returns
    -------
    steps : int
        flag and number of time steps
    time : float
        current model time
    doexit : int
        whether to exit from forecasts
    flag : int
        status flag
    """
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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
          param_int,
          param_real,
          int comm_model,
          int comm_filter,
          int comm_couple,
          int task_id,
          int n_modeltasks,
          bint in_filterpe,
          py__init_ens_pdaf,
          int in_screen
         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    filtertype : int
        type of filter
    subtype : int
        sub-type of filter
    stepnull : int
        initial time step of assimilation
    param_int : ndarray[int]
        integer parameter array
    param_real : ndarray[float]
        real parameter array
    comm_model : int
        model communicator
    comm_filter : int
        filter communicator
    comm_couple : int
        coupling communicator
    task_id : int
        id of my ensemble task
    n_modeltasks : int
        number of parallel model tasks
    in_filterpe : bool
        is my pe a filter-pe?
    py__init_ens_pdaf : func
        user-supplied routine for ensemble initialization
    in_screen : int
        control screen output:

    Returns
    -------
    param_int : ndarray[int]
        integer parameter array
    param_real : ndarray[float]
        real parameter array
    flag : int
        status flag, 0: no error, error codes:
    """
    cdef int[::1] param_int_view = np.array(param_int, dtype=np.intc).ravel(order='F')
    cdef double[::1] param_real_view = np.array(param_real).ravel(order='F')
    cdef int dim_preal, dim_pint
    dim_pint,  = param_int.shape
    dim_preal,  = param_real.shape

    PDAFcython.py__init_ens_pdaf = py__init_ens_pdaf

    cdef int flag

    c__pdaf_init (&filtertype,
                  &subtype,
                  &stepnull,
                  &param_int_view[0],
                  &dim_pint,
                  &param_real_view[0],
                  &dim_preal,
                  &comm_model,
                  &comm_filter,
                  &comm_couple,
                  &task_id,
                  &n_modeltasks,
                  &in_filterpe,
                  c__PDAFcython.c__init_ens_pdaf,
                  &in_screen,
                  &flag
                 )

    return np.asarray(param_int_view).reshape((dim_pint), order='F'), np.asarray(param_real_view).reshape((dim_preal), order='F'), flag

def local_weight (int wtype,
                  int rtype,
                  double cradius,
                  double sradius,
                  double distance,
                  a,
                  double var_obs,
                  int verbose
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    wtype : int
        type of weight function
    rtype : int
        type of regulated weighting
    cradius : float
        cut-off radius
    sradius : float
        support radius
    distance : float
        distance to observation
    a : ndarray[float]
        input matrix
    var_obs : float
        observation variance
    verbose : int
        verbosity flag

    Returns
    -------
    weight : float
        weights
    """
    cdef double[::1] a_view = np.array(a).ravel(order='F')
    cdef int ncols, nrows
    nrows, ncols,  = a.shape


    cdef double weight

    c__pdaf_local_weight (&wtype,
                          &rtype,
                          &cradius,
                          &sradius,
                          &distance,
                          &nrows,
                          &ncols,
                          &a_view[0],
                          &var_obs,
                          &weight,
                          &verbose
                         )

    return weight

def local_weights (int wtype,
                   double cradius,
                   double sradius,
                   distance,
                   int verbose
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    wtype : int
        type of weight function
         (0): unit weight (=1 up to distance=cradius)
         (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
         (2): 5th order polynomial (gaspari&cohn 1999; 0 for distance>cradius)
    cradius : float
        parameter for cut-off
    sradius : float
        support radius
    distance : ndarray[float]
        array holding distances
    verbose : int
        verbosity flag

    Returns
    -------
    weight : ndarray[float]
        array for weights
    """
    cdef double[::1] distance_view = np.array(distance).ravel(order='F')
    cdef int dim
    dim,  = distance.shape


    cdef double [::1] weight_view = np.zeros((dim)).ravel()

    c__pdaf_local_weights (&wtype,
                           &cradius,
                           &sradius,
                           &dim,
                           &distance_view[0],
                           &weight_view[0],
                           &verbose
                          )

    return np.asarray(weight_view).reshape((dim), order='F')

def prepost (py__collect_state_pdaf,
             py__distribute_state_pdaf,
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
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__next_observation_pdaf : func
        routine to provide time step, time and dimension of next observation

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int outflag

    c__pdaf_prepost (c__PDAFcython.c__collect_state_pdaf,
                     c__PDAFcython.c__distribute_state_pdaf,
                     c__PDAFcython.c__prepoststep_pdaf,
                     c__PDAFcython.c__next_observation_pdaf,
                     &outflag
                    )

    return outflag

def print_info (int printtype
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    printtype : int
        type of screen output
    """

    c__pdaf_print_info (&printtype
                       )

def put_state_3dvar (py__collect_state_pdaf,
                     py__init_dim_obs_pdaf,
                     py__obs_op_pdaf,
                     py__init_obs_pdaf,
                     py__prodrinva_pdaf,
                     py__cvt_pdaf,
                     py__cvt_adj_pdaf,
                     py__obs_op_lin_pdaf,
                     py__obs_op_adj_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
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

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_3dvar (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prodrinva_pdaf,
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
                             py__prodrinva_pdaf,
                             py__cvt_ens_pdaf,
                             py__cvt_adj_ens_pdaf,
                             py__obs_op_lin_pdaf,
                             py__obs_op_adj_pdaf,
                             py__init_obsvar_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_en3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                     c__PDAFcython.c__init_dim_obs_pdaf,
                                     c__PDAFcython.c__obs_op_pdaf,
                                     c__PDAFcython.c__init_obs_pdaf,
                                     c__PDAFcython.c__prodrinva_pdaf,
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
                              py__prodrinva_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_dim_obs_f_pdaf,
                              py__obs_op_f_pdaf,
                              py__init_obs_f_pdaf,
                              py__init_obs_l_pdaf,
                              py__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__init_obs_f_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__init_obs_f_pdaf = py__init_obs_f_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_en3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodrinva_pdaf,
                                      c__PDAFcython.c__cvt_ens_pdaf,
                                      c__PDAFcython.c__cvt_adj_ens_pdaf,
                                      c__PDAFcython.c__obs_op_lin_pdaf,
                                      c__PDAFcython.c__obs_op_adj_pdaf,
                                      c__PDAFcython.c__init_dim_obs_f_pdaf,
                                      c__PDAFcython.c__obs_op_f_pdaf,
                                      c__PDAFcython.c__init_obs_f_pdaf,
                                      c__PDAFcython.c__init_obs_l_pdaf,
                                      c__PDAFcython.c__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__add_obs_err_pdaf : func
        add obs error covariance r to hph in enkf
    py__init_obs_covar_pdaf : func
        initialize obs. error cov. matrix r in enkf

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__add_obs_err_pdaf = py__add_obs_err_pdaf
    PDAFcython.py__init_obs_covar_pdaf = py__init_obs_covar_pdaf

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
                     py__prodrinva_pdaf,
                     py__init_obsvar_pdaf
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__init_obsvar_pdaf : func
        initialize mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_estkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_pdaf,
                             c__PDAFcython.c__init_obsvar_pdaf,
                             &flag
                            )

    return flag

def put_state_etkf (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__prodrinva_pdaf,
                    py__init_obsvar_pdaf
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__init_obsvar_pdaf : func
        initialize mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_etkf (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodrinva_pdaf,
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
    py__init_obserr_f_pdaf : func
        initialize vector of observation errors
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
    PDAFcython.py__init_obserr_f_pdaf = py__init_obserr_f_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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
                              py__prodrinva_pdaf,
                              py__cvt_pdaf,
                              py__cvt_adj_pdaf,
                              py__cvt_ens_pdaf,
                              py__cvt_adj_ens_pdaf,
                              py__obs_op_lin_pdaf,
                              py__obs_op_adj_pdaf,
                              py__init_obsvar_pdaf,
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_hyb3dvar_estkf (c__PDAFcython.c__collect_state_pdaf,
                                      c__PDAFcython.c__init_dim_obs_pdaf,
                                      c__PDAFcython.c__obs_op_pdaf,
                                      c__PDAFcython.c__init_obs_pdaf,
                                      c__PDAFcython.c__prodrinva_pdaf,
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
                               py__prodrinva_pdaf,
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
                               py__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__cvt_ens_pdaf : func
        apply control vector transform matrix (ensemble)
    py__cvt_adj_ens_pdaf : func
        apply adjoint control vector transform matrix (ensemble var)
    py__cvt_pdaf : func
        apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : func
        apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : func
        linearized observation operator
    py__obs_op_adj_pdaf : func
        adjoint observation operator
    py__init_dim_obs_f_pdaf : func
        initialize dimension of observation vector
    py__obs_op_f_pdaf : func
        observation operator
    py__init_obs_f_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__cvt_ens_pdaf = py__cvt_ens_pdaf
    PDAFcython.py__cvt_adj_ens_pdaf = py__cvt_adj_ens_pdaf
    PDAFcython.py__cvt_pdaf = py__cvt_pdaf
    PDAFcython.py__cvt_adj_pdaf = py__cvt_adj_pdaf
    PDAFcython.py__obs_op_lin_pdaf = py__obs_op_lin_pdaf
    PDAFcython.py__obs_op_adj_pdaf = py__obs_op_adj_pdaf
    PDAFcython.py__init_dim_obs_f_pdaf = py__init_dim_obs_f_pdaf
    PDAFcython.py__obs_op_f_pdaf = py__obs_op_f_pdaf
    PDAFcython.py__init_obs_f_pdaf = py__init_obs_f_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

    cdef int outflag

    c__pdaf_put_state_hyb3dvar_lestkf (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_pdaf,
                                       c__PDAFcython.c__obs_op_pdaf,
                                       c__PDAFcython.c__init_obs_pdaf,
                                       c__PDAFcython.c__prodrinva_pdaf,
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
                                       c__PDAFcython.c__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__localize_covar_pdaf : func
        apply localization to hp and hph^t
    py__add_obs_err_pdaf : func
        add obs error covariance r to hph in enkf
    py__init_obs_covar_pdaf : func
        initialize obs. error cov. matrix r in enkf

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__localize_covar_pdaf = py__localize_covar_pdaf
    PDAFcython.py__add_obs_err_pdaf = py__add_obs_err_pdaf
    PDAFcython.py__init_obs_covar_pdaf = py__init_obs_covar_pdaf

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
                      py__prodrinva_l_pdaf,
                      py__init_n_domains_p_pdaf,
                      py__init_dim_l_pdaf,
                      py__init_dim_obs_l_pdaf,
                      py__g2l_state_pdaf,
                      py__l2g_state_pdaf,
                      py__g2l_obs_pdaf,
                      py__init_obsvar_pdaf,
                      py__init_obsvar_l_pdaf
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_lestkf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodrinva_l_pdaf,
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
                     py__prodrinva_l_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf,
                     py__g2l_obs_pdaf,
                     py__init_obsvar_pdaf,
                     py__init_obsvar_l_pdaf
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_letkf (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__init_obs_l_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_l_pdaf : func
        compute observation likelihood for an ensemble member
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_l_pdaf = py__likelihood_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf

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
                      py__prodrinva_l_pdaf,
                      py__prodrinva_hyb_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
    py__prodrinva_hyb_l_pdaf : func
        provide product r^-1 a on local analysis domain with hybrid weight
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance
    py__likelihood_l_pdaf : func
        compute likelihood
    py__likelihood_hyb_l_pdaf : func
        compute likelihood with hybrid weight

    Returns
    -------
    outflag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__prodrinva_hyb_l_pdaf = py__prodrinva_hyb_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf
    PDAFcython.py__likelihood_l_pdaf = py__likelihood_l_pdaf
    PDAFcython.py__likelihood_hyb_l_pdaf = py__likelihood_hyb_l_pdaf

    cdef int outflag

    c__pdaf_put_state_lknetf (c__PDAFcython.c__collect_state_pdaf,
                              c__PDAFcython.c__init_dim_obs_pdaf,
                              c__PDAFcython.c__obs_op_pdaf,
                              c__PDAFcython.c__init_obs_pdaf,
                              c__PDAFcython.c__init_obs_l_pdaf,
                              c__PDAFcython.c__prepoststep_pdaf,
                              c__PDAFcython.c__prodrinva_l_pdaf,
                              c__PDAFcython.c__prodrinva_hyb_l_pdaf,
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
                     py__prodrinva_l_pdaf,
                     py__init_n_domains_p_pdaf,
                     py__init_dim_l_pdaf,
                     py__init_dim_obs_l_pdaf,
                     py__g2l_state_pdaf,
                     py__l2g_state_pdaf,
                     py__g2l_obs_pdaf,
                     py__init_obsvar_pdaf,
                     py__init_obsvar_l_pdaf
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
    py__init_obs_pdaf : func
        initialize pe-local observation vector
    py__init_obs_l_pdaf : func
        init. observation vector on local analysis domain
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_l_pdaf : func
        provide product r^-1 a on local analysis domain
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
    py__g2l_obs_pdaf : func
        restrict full obs. vector to local analysis domain
    py__init_obsvar_pdaf : func
        initialize mean observation error variance
    py__init_obsvar_l_pdaf : func
        initialize local mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__init_obs_l_pdaf = py__init_obs_l_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_l_pdaf = py__prodrinva_l_pdaf
    PDAFcython.py__init_n_domains_p_pdaf = py__init_n_domains_p_pdaf
    PDAFcython.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    PDAFcython.py__init_dim_obs_l_pdaf = py__init_dim_obs_l_pdaf
    PDAFcython.py__g2l_state_pdaf = py__g2l_state_pdaf
    PDAFcython.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__g2l_obs_pdaf = py__g2l_obs_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf
    PDAFcython.py__init_obsvar_l_pdaf = py__init_obsvar_l_pdaf

    cdef int flag

    c__pdaf_put_state_lseik (c__PDAFcython.c__collect_state_pdaf,
                             c__PDAFcython.c__init_dim_obs_pdaf,
                             c__PDAFcython.c__obs_op_pdaf,
                             c__PDAFcython.c__init_obs_pdaf,
                             c__PDAFcython.c__init_obs_l_pdaf,
                             c__PDAFcython.c__prepoststep_pdaf,
                             c__PDAFcython.c__prodrinva_l_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_pdaf : func
        compute observation likelihood for an ensemble member

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_pdaf = py__likelihood_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__init_dim_obs_pdaf : func
        initialize dimension of observation vector
    py__obs_op_pdaf : func
        observation operator
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__likelihood_pdaf : func
        compute observation likelihood for an ensemble member

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__likelihood_pdaf = py__likelihood_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    py__collect_state_pdaf : func
        routine to collect a state vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf

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
                    py__prodrinva_pdaf
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 hv

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf

    cdef int flag

    c__pdaf_put_state_seek (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodrinva_pdaf,
                            &flag
                           )

    return flag

def put_state_seik (py__collect_state_pdaf,
                    py__init_dim_obs_pdaf,
                    py__obs_op_pdaf,
                    py__init_obs_pdaf,
                    py__prepoststep_pdaf,
                    py__prodrinva_pdaf,
                    py__init_obsvar_pdaf
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
    py__init_obs_pdaf : func
        initialize observation vector
    py__prepoststep_pdaf : func
        user supplied pre/poststep routine
    py__prodrinva_pdaf : func
        provide product r^-1 a
    py__init_obsvar_pdaf : func
        initialize mean observation error variance

    Returns
    -------
    flag : int
        status flag
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__init_dim_obs_pdaf = py__init_dim_obs_pdaf
    PDAFcython.py__obs_op_pdaf = py__obs_op_pdaf
    PDAFcython.py__init_obs_pdaf = py__init_obs_pdaf
    PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
    PDAFcython.py__prodrinva_pdaf = py__prodrinva_pdaf
    PDAFcython.py__init_obsvar_pdaf = py__init_obsvar_pdaf

    cdef int flag

    c__pdaf_put_state_seik (c__PDAFcython.c__collect_state_pdaf,
                            c__PDAFcython.c__init_dim_obs_pdaf,
                            c__PDAFcython.c__obs_op_pdaf,
                            c__PDAFcython.c__init_obs_pdaf,
                            c__PDAFcython.c__prepoststep_pdaf,
                            c__PDAFcython.c__prodrinva_pdaf,
                            c__PDAFcython.c__init_obsvar_pdaf,
                            &flag
                           )

    return flag

def reset_forget (double forget_in
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    forget_in : float
        new value of forgetting factor
    """

    c__pdaf_reset_forget (&forget_in
                         )

def sampleens (modes,
               svals,
               state,
               int verbose,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    modes : ndarray[float]
        array of eof modes
    svals : ndarray[float]
        vector of singular values
    state : ndarray[float]
        pe-local model state
    verbose : int
        verbosity flag
    flag : int
        status flag

    Returns
    -------
    modes : ndarray[float]
        array of eof modes
    state : ndarray[float]
        pe-local model state
    ens : ndarray[float]
        state ensemble
    flag : int
        status flag
    """
    cdef double[::1] modes_view = np.array(modes).ravel(order='F')
    cdef double[::1] svals_view = np.array(svals).ravel(order='F')
    cdef double[::1] state_view = np.array(state).ravel(order='F')
    cdef int dim, dim_ens
    dim, dim_ens,  = modes.shape
    dim_ens = dim_ens + 1


    cdef double [::1] ens_view = np.zeros((dim,dim_ens)).ravel()

    c__pdaf_sampleens (&dim,
                       &dim_ens,
                       &modes_view[0],
                       &svals_view[0],
                       &state_view[0],
                       &ens_view[0],
                       &verbose,
                       &flag
                      )

    return np.asarray(modes_view).reshape((dim,dim_ens-1), order='F'), np.asarray(state_view).reshape((dim), order='F'), np.asarray(ens_view).reshape((dim,dim_ens), order='F'), flag

def set_debug_flag (int debugval
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    debugval : int
        value of debugging flag; print debug information for >0
    """

    c__pdaf_set_debug_flag (&debugval
                           )

def set_ens_pointer ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 


    Returns
    -------
    c_ens_point : ndarray[float]
        pointer to smoother array
    dims : ndarray[int]
        dimension of the pointer
    status : int
        status flag
    """

    cdef double* c_ens_point
    cdef int [::1] dims_view = np.zeros((2), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_ens_pointer (&c_ens_point,
                             &dims_view[0],
                             &status
                            )

    dims = np.asarray(dims_view)
    return np.asarray(<double[:np.prod(dims)]> c_ens_point).reshape(dims, order='F'), \
           status

def set_memberid (int memberid
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    memberid : int
        index in the local ensemble

    Returns
    -------
    memberid : int
        index in the local ensemble
    """

    c__pdaf_set_memberid (&memberid
                         )

    return memberid

def set_smootherens (int maxlag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    maxlag : int
        number of past timesteps processed in sens

    Returns
    -------
    c_sens_point : ndarray[float]
        pointer to smoother array
    dims : ndarray[int]
        dimension of the pointer
    status : int
        status flag
    """

    cdef double* c_sens_point
    cdef int [::1] dims_view = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_smootherens (&c_sens_point,
                             &maxlag,
                             &dims_view[0],
                             &status
                            )

    dims = np.asarray(dims_view)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           status

def seik_ttimesa (a
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    a : ndarray[float]
        input matrix

    Returns
    -------
    b : ndarray[float]
        output matrix (ta)
    """
    cdef double[::1] a_view = np.array(a).ravel(order='F')
    cdef int dim_col, rank
    rank, dim_col,  = a.shape


    cdef double [::1] b_view = np.zeros((rank+1,dim_col)).ravel()

    c__pdaf_seik_ttimesa (&rank,
                          &dim_col,
                          &a_view[0],
                          &b_view[0]
                         )

    return np.asarray(b_view).reshape((rank+1,dim_col), order='F')

def etkf_tleft (a
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    a : ndarray[float]
        input/output matrix

    Returns
    -------
    a : ndarray[float]
        input/output matrix
    """
    cdef double[::1] a_view = np.array(a).ravel(order='F')
    cdef int dim, dim_ens
    dim_ens, dim,  = a.shape


    c__pdaf_etkf_tleft (&dim_ens,
                        &dim,
                        &a_view[0]
                       )

    return np.asarray(a_view).reshape((dim_ens,dim), order='F')

def estkf_omegaa (a
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    a : ndarray[float]
        input matrix

    Returns
    -------
    b : ndarray[float]
        output matrix (ta)
    """
    cdef double[::1] a_view = np.array(a).ravel(order='F')
    cdef int dim_col, rank
    rank, dim_col,  = a.shape


    cdef double [::1] b_view = np.zeros((rank+1,dim_col)).ravel()

    c__pdaf_estkf_omegaa (&rank,
                          &dim_col,
                          &a_view[0],
                          &b_view[0]
                         )

    return np.asarray(b_view).reshape((rank+1,dim_col), order='F')

def enkf_omega (seed,
                omega,
                double norm,
                int otype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    seed : ndarray[int]
        seed for random number generation
    omega : ndarray[float]
        random matrix
    norm : float
        norm for ensemble transformation
    otype : int
        type of omega
    screen : int
        verbosity flag

    Returns
    -------
    omega : ndarray[float]
        random matrix
    norm : float
        norm for ensemble transformation
    """
    cdef int[::1] seed_view = np.array(seed, dtype=np.intc).ravel(order='F')
    cdef double[::1] omega_view = np.array(omega).ravel(order='F')
    cdef int dim_ens, r
    dim_ens, r,  = omega.shape


    c__pdaf_enkf_omega (&seed_view[0],
                        &r,
                        &dim_ens,
                        &omega_view[0],
                        &norm,
                        &otype,
                        &screen
                       )

    return np.asarray(omega_view).reshape((dim_ens,r), order='F'), norm

def seik_omega (omega,
                int omegatype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    omega : ndarray[float]
        matrix omega
    omegatype : int
        select type of omega
    screen : int
        verbosity flag

    Returns
    -------
    omega : ndarray[float]
        matrix omega
    """
    cdef double[::1] omega_view = np.array(omega).ravel(order='F')
    cdef int rank
    rank, _,  = omega.shape
    rank = rank - 1


    c__pdaf_seik_omega (&rank,
                        &omega_view[0],
                        &omegatype,
                        &screen
                       )

    return np.asarray(omega_view).reshape((rank+1,rank), order='F')

def incremental (int steps,
                 py__dist_stateinc_pdaf
                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    steps : int
        time steps over which increment is distributed
    py__dist_stateinc_pdaf : func
        add state increment during integration
    """
    PDAFcython.py__dist_stateinc_pdaf = py__dist_stateinc_pdaf

    c__pdaf_incremental (&steps,
                         c__PDAFcython.c__dist_stateinc_pdaf
                        )

def add_increment (state_p
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ 

    Parameters
    ----------
    state_p : ndarray[float]
        state vector

    Returns
    -------
    state_p : ndarray[float]
        state vector
    """
    cdef double[::1] state_p_view = np.array(state_p).ravel(order='F')
    cdef int dim_p
    dim_p,  = state_p.shape


    c__pdaf_add_increment (&dim_p,
                           &state_p_view[0]
                          )

    return np.asarray(state_p_view).reshape((dim_p), order='F')

