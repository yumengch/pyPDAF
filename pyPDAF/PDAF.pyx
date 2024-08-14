import sys

import pyPDAF.UserFunc as PDAFcython
cimport pyPDAF.UserFunc as c__PDAFcython

import numpy as np
cimport numpy as cnp

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
                sys.stderr.write("Uncaught exception was detected on rank {}. \n".format(
                    mpi4py.MPI.COMM_WORLD.Get_rank()))

                print_exception(exctype, value, traceback)
                sys.stderr.write("\n")
                sys.stderr.flush()
            finally:
                try:
                    mpi4py.MPI.COMM_WORLD.Abort(1)
                except Exception as e:
                    sys.stderr.write("MPI Abort failed, this process will hang.\n")
                    sys.stderr.flush()
                    raise e
        else:
            sys.__excepthook__(exctype, value, traceback)
    except ImportError:
        sys.__excepthook__(exctype, value, traceback)

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_3dvar or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_en3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_en3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_enkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_etkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_hyb3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_hyb3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lenkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_letkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lnetf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lknetf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_lseik or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_netf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_pf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_seek or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_seik or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_assimilate_prepost or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int flag

    c__pdaf_assimilate_prepost (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__distribute_state_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__next_observation_pdaf,
                                &flag
                               )

    return flag

def deallocate ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_deallocate or PDAF source files 

    """
    c__pdaf_deallocate ()

def diag_effsample (double[::1] weights
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_effsample or PDAF source files 

    Parameters
    ----------
    weights : ndarray[float]
        weights of the samples; Dimension: dim_sample

    Returns
    -------
    effsample : float
        effecfive sample size
    """
    cdef int dim_sample
    dim_sample = weights.shape[0]


    cdef double effsample

    c__pdaf_diag_effsample (&dim_sample,
                            &weights[0],
                            &effsample
                           )

    return effsample

def diag_ensstats (int element,
                   double[::1] state,
                   double[:,:] ens
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_ensstats or PDAF source files 

    Parameters
    ----------
    element : int
        id of element to be used. if element=0, mean values over all elements are computed
    state : ndarray[float]
        state vector; Dimension: dim
    ens : ndarray[float]
        state ensemble; Dimension: dim,dim_ens

    Returns
    -------
    skewness : float
        skewness of ensemble
    kurtosis : float
        kurtosis of ensemble
    status : int
        status flag (0=success)
    """
    cdef double[::1] ens_f = np.asfortranarray(ens).ravel(order="F")
    cdef int dim, dim_ens
    dim = ens.shape[0]
    dim_ens = ens.shape[1]


    cdef double skewness
    cdef double kurtosis
    cdef int status

    c__pdaf_diag_ensstats (&dim,
                           &dim_ens,
                           &element,
                           &state[0],
                           &ens_f[0],
                           &skewness,
                           &kurtosis,
                           &status
                          )

    return skewness, kurtosis, status

def diag_histogram (int ncall,
                    int element,
                    double[::1] state,
                    double[:,:] ens,
                    int[::1] hist
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_histogram or PDAF source files 

    Parameters
    ----------
    ncall : int
        number of calls to routine
    element : int
        element of vector used for histogram
    state : ndarray[float]
        if element=0, all elements are usedstate vector; Dimension: dim
    ens : ndarray[float]
        state ensemble; Dimension: dim,dim_ens
    hist : ndarray[int]
        histogram about the state; Dimension: dim_ens+1

    Returns
    -------
    hist : ndarray[int]
        histogram about the state; Dimension: dim_ens+1
    delta : float
        deviation measure from flat histogram
    status : int
        status flag (0=success)
    """
    cdef double[::1] ens_f = np.asfortranarray(ens).ravel(order="F")
    cdef int dim, dim_ens
    dim = ens.shape[0]
    dim_ens = ens.shape[1]


    cdef double delta
    cdef int status

    c__pdaf_diag_histogram (&ncall,
                            &dim,
                            &dim_ens,
                            &element,
                            &state[0],
                            &ens_f[0],
                            &hist[0],
                            &delta,
                            &status
                           )

    return np.asarray(hist).reshape((dim_ens+1), order='F'), delta, status

def eofcovar (int[::1] dim_fields,
              int[::1] offsets,
              int remove_mstate,
              int do_mv,
              double[:,:] states,
              double[::1] meanstate,
              int verbose
             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_eofcovar or PDAF source files 

    Parameters
    ----------
    dim_fields : ndarray[int]
        size of each field; Dimension: nfields
    offsets : ndarray[int]
        start position of each field; Dimension: nfields
    remove_mstate : int
        1: subtract mean state from statesbefore computing eofs; 0: don't remove
    do_mv : int
        1: do multivariate scaling; 0: no scalingnfields, dim_fields and offsets are only used if do_mv=1
    states : ndarray[float]
        state perturbations; Dimension: dim_state,nstates
    meanstate : ndarray[float]
        mean state (only changed if remove_mstate=1); Dimension: dim_state
    verbose : int
        verbosity flag

    Returns
    -------
    states : ndarray[float]
        state perturbations; Dimension: dim_state,nstates
    stddev : ndarray[float]
        standard deviation of field variabilitywithout multivariate scaling (do_mv=0), it is stddev = 1.0; Dimension: nfields
    svals : ndarray[float]
        singular values divided by sqrt(nstates-1); Dimension: nstates
    svec : ndarray[float]
        singular vectors; Dimension: dim_state,nstates
    meanstate : ndarray[float]
        mean state (only changed if remove_mstate=1); Dimension: dim_state
    status : int
        status flag
    """
    cdef double[::1] states_f = np.asfortranarray(states).ravel(order="F")
    cdef int nfields, dim_state, nstates
    dim_state = states.shape[0]
    nstates = states.shape[1]
    nfields = dim_fields.shape[0]


    cdef double [::1] stddev = np.zeros((nfields), dtype=np.float64).ravel()
    cdef double [::1] svals = np.zeros((nstates), dtype=np.float64).ravel()
    cdef double [::1] svec = np.zeros((dim_state, nstates), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_eofcovar (&dim_state,
                      &nstates,
                      &nfields,
                      &dim_fields[0],
                      &offsets[0],
                      &remove_mstate,
                      &do_mv,
                      &states_f[0],
                      &stddev[0],
                      &svals[0],
                      &svec[0],
                      &meanstate[0],
                      &verbose,
                      &status
                     )

    return np.asarray(states).reshape((dim_state, nstates), order='F'), np.asarray(stddev).reshape((nfields), order='F'), np.asarray(svals).reshape((nstates), order='F'), np.asarray(svec).reshape((dim_state, nstates), order='F'), np.asarray(meanstate).reshape((dim_state), order='F'), status

def gather_dim_obs_f (int dim_obs_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_dim_obs_f or PDAF source files 

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

def gather_obs_f (double[::1] obs_p,
                  int dimobs_f
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f or PDAF source files 

    Parameters
    ----------
    obs_p : ndarray[float]
        pe-local vector; Dimension: dimobs_p
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    obs_f : ndarray[float]
        full gathered vector; Dimension: dimobs_f
    status : int
        status flag:(0) no error(1) when pdaf_gather_dim_obs_f not executed before
    """
    cdef int dimobs_p
    dimobs_p = obs_p.shape[0]


    cdef double [::1] obs_f = np.zeros((dimobs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f (&obs_p[0],
                          &dimobs_p,
                          &obs_f[0],
                          &dimobs_f,
                          &status
                         )

    return np.asarray(obs_f).reshape((dimobs_f), order='F'), status

def gather_obs_f2 (double[:,:] coords_p,
                   int dimobs_f
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f2 or PDAF source files 

    Parameters
    ----------
    coords_p : ndarray[float]
        pe-local array; Dimension: nrows,dimobs_p
    dimobs_f : int
        dimension of full gathered obs

    Returns
    -------
    coords_f : ndarray[float]
        full gathered array; Dimension: nrows,dimobs_f
    status : int
        status flag:(0) no error(1) when pdaf_gather dim_obs_f not executed before
    """
    cdef double[::1] coords_p_f = np.asfortranarray(coords_p).ravel(order="F")
    cdef int nrows, dimobs_p
    nrows = coords_p.shape[0]
    dimobs_p = coords_p.shape[1]


    cdef double [::1] coords_f = np.zeros((nrows, dimobs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f2 (&coords_p_f[0],
                           &dimobs_p,
                           &coords_f[0],
                           &dimobs_f,
                           &nrows,
                           &status
                          )

    return np.asarray(coords_f).reshape((nrows, dimobs_f), order='F'), status

def generate_obs (py__collect_state_pdaf,
                  py__distribute_state_pdaf,
                  py__init_dim_obs_f_pdaf,
                  py__obs_op_f_pdaf,
                  py__get_obs_f_pdaf,
                  py__init_obserr_f_pdaf,
                  py__prepoststep_pdaf,
                  py__next_observation_pdaf
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_generate_obs or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_assim_flag or PDAF source files 


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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_ensstats or PDAF source files 


    Returns
    -------
    dims : ndarray[int]
        dimension of pointer; Dimension: 1
    c_skew_ptr : ndarray[float]
        pointer to skewness array
    c_kurt_ptr : ndarray[float]
        pointer to kurtosis array
    status : int
        status flag
    """

    cdef int [::1] dims = np.zeros((1), dtype=np.intc).ravel()
    cdef double* c_skew_ptr
    cdef double* c_kurt_ptr
    cdef int status

    c__pdaf_get_ensstats (&dims[0],
                          &c_skew_ptr,
                          &c_kurt_ptr,
                          &status
                         )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_skew_ptr).reshape(dims, order='F'), \
           np.asarray(<double[:np.prod(dims)]> c_kurt_ptr).reshape(dims, order='F'), \
           status

def get_localfilter ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_localfilter or PDAF source files 


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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_memberid or PDAF source files 

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_obsmemberid or PDAF source files 

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_smootherens or PDAF source files 


    Returns
    -------
    c_sens_point : ndarray[float]
        pointer to smoother array
    maxlag : int
        number of past timesteps processed in sens
    dims : ndarray[int]
        dimension of pointer; Dimension: 3
    status : int
        status flag
    """

    cdef double* c_sens_point
    cdef int maxlag
    cdef int [::1] dims = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_get_smootherens (&c_sens_point,
                             &maxlag,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           maxlag, status

def get_state (int steps,
               int doexit,
               py__next_observation_pdaf,
               py__distribute_state_pdaf,
               py__prepoststep_pdaf,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_get_state or PDAF source files 

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
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
          int[::1] param_int,
          double[::1] param_real,
          int comm_model,
          int comm_filter,
          int comm_couple,
          int task_id,
          int n_modeltasks,
          bint in_filterpe,
          py__init_ens_pdaf,
          int in_screen
         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_init or PDAF source files 

    Parameters
    ----------
    filtertype : int
        type of filter
    subtype : int
        sub-type of filter
    stepnull : int
        initial time step of assimilation
    param_int : ndarray[int]
        integer parameter array; Dimension: dim_pint
    param_real : ndarray[float]
        real parameter array; Dimension: dim_preal
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
        integer parameter array; Dimension: dim_pint
    param_real : ndarray[float]
        real parameter array; Dimension: dim_preal
    flag : int
        status flag, 0: no error, error codes:
    """
    cdef int dim_pint, dim_preal
    dim_pint = param_int.shape[0]
    dim_preal = param_real.shape[0]

    c__PDAFcython.init_ens_pdaf = <void*>py__init_ens_pdaf

    cdef int flag

    if (filtertype == 0) or (filtertype == 200 and subtype == 0):
        c__pdaf_init (&filtertype, &subtype, &stepnull,
                      &param_int[0], &dim_pint,
                      &param_real[0], &dim_preal,
                      &comm_model, &comm_filter, &comm_couple,
                      &task_id, &n_modeltasks, &in_filterpe,
                      c__PDAFcython.c__init_ens_pdaf_single_member,
                      &in_screen, &flag)
    else:
        c__pdaf_init (&filtertype,
                      &subtype,
                      &stepnull,
                      &param_int[0],
                      &dim_pint,
                      &param_real[0],
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

    return np.asarray(param_int).reshape((dim_pint), order='F'), np.asarray(param_real).reshape((dim_preal), order='F'), flag

def local_weight (int wtype,
                  int rtype,
                  double cradius,
                  double sradius,
                  double distance,
                  double[:,:] a,
                  double var_obs,
                  int verbose
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_local_weight or PDAF source files 

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
        input matrix; Dimension: nrows,ncols
    var_obs : float
        observation variance
    verbose : int
        verbosity flag

    Returns
    -------
    weight : float
        weights
    """
    cdef double[::1] a_f = np.asfortranarray(a).ravel(order="F")
    cdef int nrows, ncols
    nrows = a.shape[0]
    ncols = a.shape[1]


    cdef double weight

    c__pdaf_local_weight (&wtype,
                          &rtype,
                          &cradius,
                          &sradius,
                          &distance,
                          &nrows,
                          &ncols,
                          &a_f[0],
                          &var_obs,
                          &weight,
                          &verbose
                         )

    return weight

def print_info (int printtype
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_info or PDAF source files 

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_3dvar or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_en3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_en3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_enkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_etkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_generate_obs or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.init_obserr_f_pdaf = <void*>py__init_obserr_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_hyb3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_hyb3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.init_obs_f_pdaf = <void*>py__init_obs_f_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lenkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    c__PDAFcython.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_letkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lnetf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lknetf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf
    c__PDAFcython.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    c__PDAFcython.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_lseik or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.init_obs_l_pdaf = <void*>py__init_obs_l_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.g2l_obs_pdaf = <void*>py__g2l_obs_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf
    c__PDAFcython.init_obsvar_l_pdaf = <void*>py__init_obsvar_l_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_netf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_pf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.likelihood_pdaf = <void*>py__likelihood_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_prepost or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_seek or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_put_state_seik or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    c__PDAFcython.init_obsvar_pdaf = <void*>py__init_obsvar_pdaf

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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_reset_forget or PDAF source files 

    Parameters
    ----------
    forget_in : float
        new value of forgetting factor
    """

    c__pdaf_reset_forget (&forget_in
                         )

def sampleens (double[:,:] modes,
               double[::1] svals,
               double[::1] state,
               int verbose,
               int flag
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_SampleEns or PDAF source files 

    Parameters
    ----------
    modes : ndarray[float]
        array of eof modes; Dimension: dim,dim_ens-1
    svals : ndarray[float]
        vector of singular values; Dimension: dim_ens-1
    state : ndarray[float]
        pe-local model state; Dimension: dim
    verbose : int
        verbosity flag
    flag : int
        status flag

    Returns
    -------
    modes : ndarray[float]
        array of eof modes; Dimension: dim,dim_ens-1
    state : ndarray[float]
        pe-local model state; Dimension: dim
    ens : ndarray[float]
        state ensemble; Dimension: dim,dim_ens
    flag : int
        status flag
    """
    cdef double[::1] modes_f = np.asfortranarray(modes).ravel(order="F")
    cdef int dim, dim_ens
    dim = modes.shape[0]
    dim_ens = modes.shape[1]
    dim_ens = dim_ens + 1


    cdef double [::1] ens = np.zeros((dim, dim_ens), dtype=np.float64).ravel()

    c__pdaf_sampleens (&dim,
                       &dim_ens,
                       &modes_f[0],
                       &svals[0],
                       &state[0],
                       &ens[0],
                       &verbose,
                       &flag
                      )

    return np.asarray(modes).reshape((dim, dim_ens-1), order='F'), np.asarray(state).reshape((dim), order='F'), np.asarray(ens).reshape((dim, dim_ens), order='F'), flag

def set_debug_flag (int debugval
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_debug_flag or PDAF source files 

    Parameters
    ----------
    debugval : int
        value of debugging flag; print debug information for >0
    """

    c__pdaf_set_debug_flag (&debugval
                           )

def set_ens_pointer ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_ens_pointer or PDAF source files 


    Returns
    -------
    c_ens_point : ndarray[float]
        pointer to smoother array
    dims : ndarray[int]
        dimension of the pointer; Dimension: 2
    status : int
        status flag
    """

    cdef double* c_ens_point
    cdef int [::1] dims = np.zeros((2), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_ens_pointer (&c_ens_point,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_ens_point).reshape(dims, order='F'), \
           status

def set_smootherens (int maxlag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_smootherens or PDAF source files 

    Parameters
    ----------
    maxlag : int
        number of past timesteps processed in sens

    Returns
    -------
    c_sens_point : ndarray[float]
        pointer to smoother array
    dims : ndarray[int]
        dimension of the pointer; Dimension: 3
    status : int
        status flag
    """

    cdef double* c_sens_point
    cdef int [::1] dims = np.zeros((3), dtype=np.intc).ravel()
    cdef int status

    c__pdaf_set_smootherens (&c_sens_point,
                             &maxlag,
                             &dims[0],
                             &status
                            )

    dims = np.asarray(dims)
    return np.asarray(<double[:np.prod(dims)]> c_sens_point).reshape(dims, order='F'), \
           status

def seik_ttimesa (double[:,:] a
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_seik_TtimesA or PDAF source files 

    Parameters
    ----------
    a : ndarray[float]
        input matrix; Dimension: rank,dim_col

    Returns
    -------
    b : ndarray[float]
        output matrix (ta); Dimension: rank+1,dim_col
    """
    cdef double[::1] a_f = np.asfortranarray(a).ravel(order="F")
    cdef int rank, dim_col
    rank = a.shape[0]
    dim_col = a.shape[1]


    cdef double [::1] b = np.zeros((rank+1, dim_col), dtype=np.float64).ravel()

    c__pdaf_seik_ttimesa (&rank,
                          &dim_col,
                          &a_f[0],
                          &b[0]
                         )

    return np.asarray(b).reshape((rank+1, dim_col), order='F')

def etkf_tleft (double[:,:] a
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_etkf_Tleft or PDAF source files 

    Parameters
    ----------
    a : ndarray[float]
        input/output matrix; Dimension: dim_ens,dim

    Returns
    -------
    a : ndarray[float]
        input/output matrix; Dimension: dim_ens,dim
    """
    cdef double[::1] a_f = np.asfortranarray(a).ravel(order="F")
    cdef int dim_ens, dim
    dim_ens = a.shape[0]
    dim = a.shape[1]


    c__pdaf_etkf_tleft (&dim_ens,
                        &dim,
                        &a_f[0]
                       )

    return np.asarray(a).reshape((dim_ens, dim), order='F')

def estkf_omegaa (double[:,:] a
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_estkf_OmegaA or PDAF source files 

    Parameters
    ----------
    a : ndarray[float]
        input matrix; Dimension: rank,dim_col

    Returns
    -------
    b : ndarray[float]
        output matrix (ta); Dimension: rank+1,dim_col
    """
    cdef double[::1] a_f = np.asfortranarray(a).ravel(order="F")
    cdef int rank, dim_col
    rank = a.shape[0]
    dim_col = a.shape[1]


    cdef double [::1] b = np.zeros((rank+1, dim_col), dtype=np.float64).ravel()

    c__pdaf_estkf_omegaa (&rank,
                          &dim_col,
                          &a_f[0],
                          &b[0]
                         )

    return np.asarray(b).reshape((rank+1, dim_col), order='F')

def enkf_omega (int[::1] seed,
                double[:,:] omega,
                double norm,
                int otype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_enkf_omega or PDAF source files 

    Parameters
    ----------
    seed : ndarray[int]
        seed for random number generation; Dimension: 4
    omega : ndarray[float]
        random matrix; Dimension: dim_ens,r
    norm : float
        norm for ensemble transformation
    otype : int
        type of omega
    screen : int
        verbosity flag

    Returns
    -------
    omega : ndarray[float]
        random matrix; Dimension: dim_ens,r
    norm : float
        norm for ensemble transformation
    """
    cdef double[::1] omega_f = np.asfortranarray(omega).ravel(order="F")
    cdef int dim_ens, r
    dim_ens = omega.shape[0]
    r = omega.shape[1]


    c__pdaf_enkf_omega (&seed[0],
                        &r,
                        &dim_ens,
                        &omega_f[0],
                        &norm,
                        &otype,
                        &screen
                       )

    return np.asarray(omega).reshape((dim_ens, r), order='F'), norm

def seik_omega (double[:,:] omega,
                int omegatype,
                int screen
               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_seik_omega or PDAF source files 

    Parameters
    ----------
    omega : ndarray[float]
        matrix omega; Dimension: rank+1,rank
    omegatype : int
        select type of omega
    screen : int
        verbosity flag

    Returns
    -------
    omega : ndarray[float]
        matrix omega; Dimension: rank+1,rank
    """
    cdef double[::1] omega_f = np.asfortranarray(omega).ravel(order="F")
    cdef int rank
    rank = omega.shape[0]
    _ = omega.shape[1]
    rank = rank - 1


    c__pdaf_seik_omega (&rank,
                        &omega_f[0],
                        &omegatype,
                        &screen
                       )

    return np.asarray(omega).reshape((rank+1, rank), order='F')

def incremental (int steps,
                 py__dist_stateinc_pdaf
                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_incremental or PDAF source files 

    Parameters
    ----------
    steps : int
        time steps over which increment is distributed
    py__dist_stateinc_pdaf : func
        add state increment during integration
    """
    c__PDAFcython.dist_stateinc_pdaf = <void*>py__dist_stateinc_pdaf

    c__pdaf_incremental (&steps,
                         c__PDAFcython.c__dist_stateinc_pdaf
                        )

def add_increment (double[::1] state_p
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_add_increment or PDAF source files 

    Parameters
    ----------
    state_p : ndarray[float]
        state vector; Dimension: dim_p

    Returns
    -------
    state_p : ndarray[float]
        state vector; Dimension: dim_p
    """
    cdef int dim_p
    dim_p = state_p.shape[0]


    c__pdaf_add_increment (&dim_p,
                           &state_p[0]
                          )

    return np.asarray(state_p).reshape((dim_p), order='F')

def local_weights (int wtype,
                   double cradius,
                   double sradius,
                   double[::1] distance,
                   int verbose
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_local_weights or PDAF source files 

    Parameters
    ----------
    wtype : int
        type of weight function(0): unit weight (=1 up to distance=cradius)(1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)(2): 5th order polynomial (gaspari&cohn 1999; 0 for distance>cradius)
    cradius : float
        parameter for cut-off
    sradius : float
        support radius
    distance : ndarray[float]
        array holding distances; Dimension: dim
    verbose : int
        verbosity flag

    Returns
    -------
    weight : ndarray[float]
        array for weights; Dimension: dim
    """
    cdef int dim
    dim = distance.shape[0]


    cdef double [::1] weight = np.zeros((dim), dtype=np.float64).ravel()

    c__pdaf_local_weights (&wtype,
                           &cradius,
                           &sradius,
                           &dim,
                           &distance[0],
                           &weight[0],
                           &verbose
                          )

    return np.asarray(weight).reshape((dim), order='F')

def diag_crps (int element,
               double[:,:] oens,
               double[::1] obs
              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_diag_CRPS or PDAF source files 

    Parameters
    ----------
    element : int
        id of element to be used. if element=0, mean values over all elements are computed
    oens : ndarray[float]
        state ensemble; Dimension: dim,dim_ens
    obs : ndarray[float]
        state ensemble; Dimension: dim

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
    cdef double[::1] oens_f = np.asfortranarray(oens).ravel(order="F")
    cdef int dim, dim_ens
    dim = oens.shape[0]
    dim_ens = oens.shape[1]


    cdef double crps
    cdef double reli
    cdef double resol
    cdef double uncert
    cdef int status

    c__pdaf_diag_crps (&dim,
                       &dim_ens,
                       &element,
                       &oens_f[0],
                       &obs[0],
                       &crps,
                       &reli,
                       &resol,
                       &uncert,
                       &status
                      )

    return crps, reli, resol, uncert, status

def force_analysis ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_force_analysis or PDAF source files 

    """
    c__pdaf_force_analysis ()

def gather_obs_f2_flex (int dim_obs_f,
                        double[:,:] coords_p
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f2_flex or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        full observation dimension
    coords_p : ndarray[float]
        pe-local array; Dimension: nrows,dim_obs_p

    Returns
    -------
    coords_f : ndarray[float]
        full gathered array; Dimension: nrows,dim_obs_f
    status : int
        status flag: (0) no error
    """
    cdef double[::1] coords_p_f = np.asfortranarray(coords_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = coords_p.shape[0]
    dim_obs_p = coords_p.shape[1]


    cdef double [::1] coords_f = np.zeros((nrows, dim_obs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f2_flex (&dim_obs_p,
                                &dim_obs_f,
                                &coords_p_f[0],
                                &coords_f[0],
                                &nrows,
                                &status
                               )

    return np.asarray(coords_f).reshape((nrows, dim_obs_f), order='F'), status

def gather_obs_f_flex (int dim_obs_f,
                       double[::1] obs_p
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_gather_obs_f_flex or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        full observation dimension
    obs_p : ndarray[float]
        pe-local vector; Dimension: dim_obs_p

    Returns
    -------
    obs_f : ndarray[float]
        full gathered vector; Dimension: dim_obs_f
    status : int
        status flag: (0) no error
    """
    cdef int dim_obs_p
    dim_obs_p = obs_p.shape[0]


    cdef double [::1] obs_f = np.zeros((dim_obs_f), dtype=np.float64).ravel()
    cdef int status

    c__pdaf_gather_obs_f_flex (&dim_obs_p,
                               &dim_obs_f,
                               &obs_p[0],
                               &obs_f[0],
                               &status
                              )

    return np.asarray(obs_f).reshape((dim_obs_f), order='F'), status

def prepost (py__collect_state_pdaf,
             py__distribute_state_pdaf,
             py__prepoststep_pdaf,
             py__next_observation_pdaf
            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_prepost or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

    cdef int outflag

    c__pdaf_prepost (c__PDAFcython.c__collect_state_pdaf,
                     c__PDAFcython.c__distribute_state_pdaf,
                     c__PDAFcython.c__prepoststep_pdaf,
                     c__PDAFcython.c__next_observation_pdaf,
                     &outflag
                    )

    return outflag

def set_memberid (int memberid
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_memberid or PDAF source files 

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

def set_comm_pdaf (int in_comm_pdaf
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_comm_pdaf or PDAF source files 

    Parameters
    ----------
    in_comm_pdaf : int
        mpi communicator for pdaf
    """

    c__pdaf_set_comm_pdaf (&in_comm_pdaf
                          )

def set_offline_mode (int screen
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_set_offline_mode or PDAF source files 

    Parameters
    ----------
    screen : int
        verbosity flag
    """

    c__pdaf_set_offline_mode (&screen
                             )

def print_domain_stats (int n_domains_p
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_domain_stats or PDAF source files 

    Parameters
    ----------
    n_domains_p : int
        number of pe-local analysis domains
    """

    c__pdaf_print_domain_stats (&n_domains_p
                               )

def init_local_obsstats ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_init_local_obsstats or PDAF source files 

    """
    c__pdaf_init_local_obsstats ()

def incr_local_obsstats (int dim_obs_l
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_incr_local_obsstats or PDAF source files 

    Parameters
    ----------
    dim_obs_l : int
        number of locally assimilated observations
    """

    c__pdaf_incr_local_obsstats (&dim_obs_l
                                )

def print_local_obsstats (int screen
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_print_local_obsstats or PDAF source files 

    Parameters
    ----------
    screen : int
        verbosity flag

    Returns
    -------
    n_domains_with_obs : int
        number of domains with observations
    """

    cdef int n_domains_with_obs

    c__pdaf_print_local_obsstats (&screen,
                                  &n_domains_with_obs
                                 )

    return n_domains_with_obs

def omit_obs_omi (double[::1] state_p,
                  double[:,:] ens_p,
                  double[::1] obs_p,
                  py__init_obs_pdaf,
                  py__obs_op_pdaf,
                  int compute_mean,
                  int screen
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAF_omit_obs_omi or PDAF source files 

    Parameters
    ----------
    state_p : ndarray[float]
        on exit: pe-local forecast mean state; Dimension: dim_p
    ens_p : ndarray[float]
        pe-local state ensemble; Dimension: dim_p,dim_ens
    obs_p : ndarray[float]
        pe-local observation vector; Dimension: dim_obs_p
    py__init_obs_pdaf : func
        initialize observation vector
    py__obs_op_pdaf : func
        observation operator
    compute_mean : int
        (1) compute mean; (0) state_p holds mean
    screen : int
        verbosity flag

    Returns
    -------
    state_p : ndarray[float]
        on exit: pe-local forecast mean state; Dimension: dim_p
    obs_p : ndarray[float]
        pe-local observation vector; Dimension: dim_obs_p
    """
    cdef double[::1] ens_p_f = np.asfortranarray(ens_p).ravel(order="F")
    cdef int dim_p, dim_ens, dim_obs_p
    dim_p = ens_p.shape[0]
    dim_ens = ens_p.shape[1]
    dim_obs_p = obs_p.shape[0]

    c__PDAFcython.init_obs_pdaf = <void*>py__init_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf

    c__pdaf_omit_obs_omi (&dim_p,
                          &dim_obs_p,
                          &dim_ens,
                          &state_p[0],
                          &ens_p_f[0],
                          &obs_p[0],
                          c__PDAFcython.c__init_obs_pdaf,
                          c__PDAFcython.c__obs_op_pdaf,
                          &compute_mean,
                          &screen
                         )

    return np.asarray(state_p).reshape((dim_p), order='F'), np.asarray(obs_p).reshape((dim_obs_p), order='F')

def omi_init (int n_obs
             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init or PDAF source files 

    Parameters
    ----------
    n_obs : int
        number of observations
    """

    c__pdafomi_init (&n_obs
                    )

def omi_set_doassim (int i_obs,
                     int doassim
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_doassim or PDAF source files 

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

def omi_set_disttype (int i_obs,
                      int disttype
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_disttype or PDAF source files 

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

def omi_set_ncoord (int i_obs,
                    int ncoord
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_ncoord or PDAF source files 

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

def omi_set_id_obs_p (int i_obs,
                      int[:,:] id_obs_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_id_obs_p or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    id_obs_p : ndarray[int]
        setter value; Dimension: nrows,dim_obs_p
    """
    cdef int[::1] id_obs_p_f = np.asfortranarray(id_obs_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = id_obs_p.shape[0]
    dim_obs_p = id_obs_p.shape[1]


    c__pdafomi_set_id_obs_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &id_obs_p_f[0]
                            )

def omi_set_icoeff_p (int i_obs,
                      double[:,:] icoeff_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_icoeff_p or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    icoeff_p : ndarray[float]
        setter value; Dimension: nrows,dim_obs_p
    """
    cdef double[::1] icoeff_p_f = np.asfortranarray(icoeff_p).ravel(order="F")
    cdef int nrows, dim_obs_p
    nrows = icoeff_p.shape[0]
    dim_obs_p = icoeff_p.shape[1]


    c__pdafomi_set_icoeff_p (&i_obs,
                             &nrows,
                             &dim_obs_p,
                             &icoeff_p_f[0]
                            )

def omi_set_domainsize (int i_obs,
                        double[::1] domainsize
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_domainsize or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    domainsize : ndarray[float]
        setter value; Dimension: ncoord
    """
    cdef int ncoord
    ncoord = domainsize.shape[0]


    c__pdafomi_set_domainsize (&i_obs,
                               &ncoord,
                               &domainsize[0]
                              )

def omi_set_obs_err_type (int i_obs,
                          int obs_err_type
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_obs_err_type or PDAF source files 

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

def omi_set_use_global_obs (int i_obs,
                            int use_global_obs
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_use_global_obs or PDAF source files 

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

def omi_set_inno_omit (int i_obs,
                       double inno_omit
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_inno_omit or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit : float
        setter value
    """

    c__pdafomi_set_inno_omit (&i_obs,
                              &inno_omit
                             )

def omi_set_inno_omit_ivar (int i_obs,
                            double inno_omit_ivar
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_inno_omit_ivar or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    inno_omit_ivar : float
        setter value
    """

    c__pdafomi_set_inno_omit_ivar (&i_obs,
                                   &inno_omit_ivar
                                  )

def omi_gather_obs (int i_obs,
                    double[::1] obs_p,
                    double[::1] ivar_obs_p,
                    double[:,:] ocoord_p,
                    double cradius
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/pdafomi_gather_obs or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    obs_p : ndarray[float]
        pe-local observation vector; Dimension: dim_obs_p
    ivar_obs_p : ndarray[float]
        pe-local inverse observation error variance; Dimension: dim_obs_p
    ocoord_p : ndarray[float]
        pe-local observation coordinates; Dimension: thisobs(i_obs)%ncoord,dim_obs_p
    cradius : float
        localization radius

    Returns
    -------
    dim_obs : int
        full number of observations
    """
    cdef double[::1] ocoord_p_f = np.asfortranarray(ocoord_p).ravel(order="F")
    cdef int dim_obs_p
    _ = ocoord_p.shape[0]
    dim_obs_p = ocoord_p.shape[1]


    cdef int dim_obs

    c__pdafomi_gather_obs (&i_obs,
                           &dim_obs_p,
                           &obs_p[0],
                           &ivar_obs_p[0],
                           &ocoord_p_f[0],
                           &cradius,
                           &dim_obs
                          )

    return dim_obs

def omi_gather_obsstate (int i_obs,
                         double[::1] obsstate_p,
                         double[::1] obsstate_f
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_gather_obsstate or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    obsstate_p : ndarray[float]
        vector of process-local observed state; Dimension: thisobs(i_obs)%dim_obs_p
    obsstate_f : ndarray[float]
        full observed vector for all types; Dimension: nobs_f_all

    Returns
    -------
    obsstate_f : ndarray[float]
        full observed vector for all types; Dimension: nobs_f_all
    """
    cdef int nobs_f_all
    nobs_f_all = obsstate_f.shape[0]


    c__pdafomi_gather_obsstate (&i_obs,
                                &obsstate_p[0],
                                &obsstate_f[0],
                                &nobs_f_all
                               )

    return np.asarray(obsstate_f).reshape((nobs_f_all), order='F')

def omi_set_domain_limits (double[:,:] lim_coords
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_domain_limits or PDAF source files 

    Parameters
    ----------
    lim_coords : ndarray[float]
        geographic coordinate array (1: longitude, 2: latitude); Dimension: 2,2
    """
    cdef double[::1] lim_coords_f = np.asfortranarray(lim_coords).ravel(order="F")

    c__pdafomi_set_domain_limits (&lim_coords_f[0]
                                 )

def omi_set_debug_flag (int debugval
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_set_debug_flag or PDAF source files 

    Parameters
    ----------
    debugval : int
        value for debugging flag
    """

    c__pdafomi_set_debug_flag (&debugval
                              )

def omi_deallocate_obs (int i_obs
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_deallocate_obs or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    """

    c__pdafomi_deallocate_obs (&i_obs
                              )

def omi_obs_op_gridpoint (int i_obs,
                          double[::1] state_p,
                          double[::1] obs_f_all
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_gridpoint or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_gridpoint (&i_obs,
                                 &state_p[0],
                                 &dim_p,
                                 &obs_f_all[0],
                                 &nobs_f_all
                                )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_gridavg (int i_obs,
                        int nrows,
                        double[::1] state_p,
                        double[::1] obs_f_all
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_gridavg or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_gridavg (&i_obs,
                               &nrows,
                               &state_p[0],
                               &dim_p,
                               &obs_f_all[0],
                               &nobs_f_all
                              )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_interp_lin (int i_obs,
                           int nrows,
                           double[::1] state_p,
                           double[::1] obs_f_all
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_interp_lin or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_interp_lin (&i_obs,
                                  &nrows,
                                  &state_p[0],
                                  &dim_p,
                                  &obs_f_all[0],
                                  &nobs_f_all
                                 )

    return np.asarray(obs_f_all).reshape((nobs_f_all), order='F')

def omi_obs_op_adj_gridavg (int i_obs,
                            int nrows,
                            double[::1] state_p,
                            double[::1] obs_f_all
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_gridavg or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_gridavg (&i_obs,
                                   &nrows,
                                   &state_p[0],
                                   &dim_p,
                                   &obs_f_all[0],
                                   &nobs_f_all
                                  )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_obs_op_adj_gridpoint (int i_obs,
                              double[::1] state_p,
                              double[::1] obs_f_all
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_gridpoint or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_gridpoint (&i_obs,
                                     &state_p[0],
                                     &dim_p,
                                     &obs_f_all[0],
                                     &nobs_f_all
                                    )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_obs_op_adj_interp_lin (int i_obs,
                               int nrows,
                               double[::1] state_p,
                               double[::1] obs_f_all
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obs_op_adj_interp_lin or PDAF source files 

    Parameters
    ----------
    i_obs : int
        index of observations
    nrows : int
        number of values to be averaged
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    obs_f_all : ndarray[float]
        full observed state for all observation types (nobs_f_all); Dimension: nobs_f_all

    Returns
    -------
    state_p : ndarray[float]
        pe-local model state (dim_p); Dimension: dim_p
    """
    cdef int dim_p, nobs_f_all
    dim_p = state_p.shape[0]
    nobs_f_all = obs_f_all.shape[0]


    c__pdafomi_obs_op_adj_interp_lin (&i_obs,
                                      &nrows,
                                      &state_p[0],
                                      &dim_p,
                                      &obs_f_all[0],
                                      &nobs_f_all
                                     )

    return np.asarray(state_p).reshape((dim_p), order='F')

def omi_get_interp_coeff_tri (double[:,:] gpc,
                              double[::1] oc,
                              double[::1] icoeff
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_tri or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points; dim(3,2); Dimension: 3,2
    oc : ndarray[float]
        3 rows; each containing lon and lat coordinatescoordinates of observation; dim(2); Dimension: 2
    icoeff : ndarray[float]
        interpolation coefficients; dim(3); Dimension: 3

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients; dim(3); Dimension: 3
    """
    cdef double[::1] gpc_f = np.asfortranarray(gpc).ravel(order="F")

    c__pdafomi_get_interp_coeff_tri (&gpc_f[0],
                                     &oc[0],
                                     &icoeff[0]
                                    )

    return np.asarray(icoeff).reshape((3), order='F')

def omi_get_interp_coeff_lin1d (double[::1] gpc,
                                double oc,
                                double[::1] icoeff
                               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_lin1D or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points (dim=2); Dimension: 2
    oc : float
        coordinates of observation
    icoeff : ndarray[float]
        interpolation coefficients (dim=2); Dimension: 2

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients (dim=2); Dimension: 2
    """

    c__pdafomi_get_interp_coeff_lin1d (&gpc[0],
                                       &oc,
                                       &icoeff[0]
                                      )

    return np.asarray(icoeff).reshape((2), order='F')

def omi_get_interp_coeff_lin (double[:,:] gpc,
                              double[::1] oc,
                              double[::1] icoeff
                             ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_get_interp_coeff_lin or PDAF source files 

    Parameters
    ----------
    gpc : ndarray[float]
        coordinates of grid points; Dimension: num_gp,n_dim
    oc : ndarray[float]
        coordinates of observation; Dimension: n_dim
    icoeff : ndarray[float]
        interpolation coefficients (num_gp); Dimension: num_gp

    Returns
    -------
    icoeff : ndarray[float]
        interpolation coefficients (num_gp); Dimension: num_gp
    """
    cdef double[::1] gpc_f = np.asfortranarray(gpc).ravel(order="F")
    cdef int num_gp, n_dim
    num_gp = gpc.shape[0]
    n_dim = gpc.shape[1]


    c__pdafomi_get_interp_coeff_lin (&num_gp,
                                     &n_dim,
                                     &gpc_f[0],
                                     &oc[0],
                                     &icoeff[0]
                                    )

    return np.asarray(icoeff).reshape((num_gp), order='F')

def omi_assimilate_3dvar (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_3dvar or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_en3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_en3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_en3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_en3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_global (py__collect_state_pdaf,
                           py__distribute_state_pdaf,
                           py__init_dim_obs_pdaf,
                           py__obs_op_pdaf,
                           py__prepoststep_pdaf,
                           py__next_observation_pdaf
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_global or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_hyb3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_hyb3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_hyb3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_hyb3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_lenkf (py__collect_state_pdaf,
                          py__distribute_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__prepoststep_pdaf,
                          py__localize_covar_pdaf,
                          py__next_observation_pdaf
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_lenkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_assimilate_local (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_assimilate_local or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_generate_obs (py__collect_state_pdaf,
                      py__distribute_state_pdaf,
                      py__init_dim_obs_f_pdaf,
                      py__obs_op_f_pdaf,
                      py__get_obs_f_pdaf,
                      py__prepoststep_pdaf,
                      py__next_observation_pdaf
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_generate_obs or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.next_observation_pdaf = <void*>py__next_observation_pdaf

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

def omi_put_state_3dvar (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__cvt_pdaf,
                         py__cvt_adj_pdaf,
                         py__obs_op_lin_pdaf,
                         py__obs_op_adj_pdaf,
                         py__prepoststep_pdaf,
                         int outflag
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_3dvar or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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

def omi_put_state_en3dvar_estkf (py__collect_state_pdaf,
                                 py__init_dim_obs_pdaf,
                                 py__obs_op_pdaf,
                                 py__cvt_ens_pdaf,
                                 py__cvt_adj_ens_pdaf,
                                 py__obs_op_lin_pdaf,
                                 py__obs_op_adj_pdaf,
                                 py__prepoststep_pdaf,
                                 int outflag
                                ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_en3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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

def omi_put_state_en3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_en3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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

def omi_put_state_generate_obs (py__collect_state_pdaf,
                                py__init_dim_obs_f_pdaf,
                                py__obs_op_f_pdaf,
                                py__get_obs_f_pdaf,
                                py__prepoststep_pdaf
                               ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_generate_obs or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_generate_obs (c__PDAFcython.c__collect_state_pdaf,
                                       c__PDAFcython.c__init_dim_obs_f_pdaf,
                                       c__PDAFcython.c__obs_op_f_pdaf,
                                       c__PDAFcython.c__get_obs_f_pdaf,
                                       c__PDAFcython.c__prepoststep_pdaf,
                                       &flag
                                      )

    return flag

def omi_put_state_global (py__collect_state_pdaf,
                          py__init_dim_obs_pdaf,
                          py__obs_op_pdaf,
                          py__prepoststep_pdaf
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_global or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

    cdef int flag

    c__pdafomi_put_state_global (c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__init_dim_obs_pdaf,
                                 c__PDAFcython.c__obs_op_pdaf,
                                 c__PDAFcython.c__prepoststep_pdaf,
                                 &flag
                                )

    return flag

def omi_put_state_hyb3dvar_estkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_hyb3dvar_estkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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

def omi_put_state_hyb3dvar_lestkf (py__collect_state_pdaf,
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
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_hyb3dvar_lestkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_f_pdaf = <void*>py__init_dim_obs_f_pdaf
    c__PDAFcython.obs_op_f_pdaf = <void*>py__obs_op_f_pdaf
    c__PDAFcython.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    c__PDAFcython.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    c__PDAFcython.cvt_pdaf = <void*>py__cvt_pdaf
    c__PDAFcython.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    c__PDAFcython.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    c__PDAFcython.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf

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

def omi_put_state_lenkf (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__prepoststep_pdaf,
                         py__localize_covar_pdaf
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_lenkf or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.localize_covar_pdaf = <void*>py__localize_covar_pdaf

    cdef int flag

    c__pdafomi_put_state_lenkf (c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__init_dim_obs_pdaf,
                                c__PDAFcython.c__obs_op_pdaf,
                                c__PDAFcython.c__prepoststep_pdaf,
                                c__PDAFcython.c__localize_covar_pdaf,
                                &flag
                               )

    return flag

def omi_put_state_local (py__collect_state_pdaf,
                         py__init_dim_obs_pdaf,
                         py__obs_op_pdaf,
                         py__prepoststep_pdaf,
                         py__init_n_domains_p_pdaf,
                         py__init_dim_l_pdaf,
                         py__init_dim_obs_l_pdaf,
                         py__g2l_state_pdaf,
                         py__l2g_state_pdaf
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_put_state_local or PDAF source files 

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
    c__PDAFcython.collect_state_pdaf = <void*>py__collect_state_pdaf
    c__PDAFcython.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    c__PDAFcython.obs_op_pdaf = <void*>py__obs_op_pdaf
    c__PDAFcython.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    c__PDAFcython.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    c__PDAFcython.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    c__PDAFcython.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    c__PDAFcython.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    c__PDAFcython.l2g_state_pdaf = <void*>py__l2g_state_pdaf

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

def omi_init_obs_f_cb (int step,
                       int dim_obs_f,
                       double[::1] observation_f
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obs_f_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        dimension of full observation vector
    observation_f : ndarray[float]
        full observation vector; Dimension: dim_obs_f

    Returns
    -------
    observation_f : ndarray[float]
        full observation vector; Dimension: dim_obs_f
    """
    c__pdafomi_init_obs_f_cb (&step,
                              &dim_obs_f,
                              &observation_f[0]
                             )

    return np.asarray(observation_f).reshape((dim_obs_f), order='F')

def omi_init_obsvar_cb (int step,
                        int dim_obs_p,
                        double[::1] obs_p,
                        double meanvar
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obsvar_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        pe-local dimension of observation vector
    obs_p : ndarray[float]
        pe-local observation vector; Dimension: dim_obs_p
    meanvar : float
        mean observation error variance

    Returns
    -------
    meanvar : float
        mean observation error variance
    """
    c__pdafomi_init_obsvar_cb (&step,
                               &dim_obs_p,
                               &obs_p[0],
                               &meanvar
                              )

    return meanvar

def omi_g2l_obs_cb (int domain_p,
                    int step,
                    int dim_obs_f,
                    int dim_obs_l,
                    double[::1] ostate_f,
                    double[::1] ostate_l
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_g2l_obs_cb or PDAF source files 

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
        full pe-local obs.ervation vector; Dimension: dim_obs_f
    ostate_l : ndarray[float]
        observation vector on local domain; Dimension: dim_obs_l

    Returns
    -------
    ostate_l : ndarray[float]
        observation vector on local domain; Dimension: dim_obs_l
    """
    c__pdafomi_g2l_obs_cb (&domain_p,
                           &step,
                           &dim_obs_f,
                           &dim_obs_l,
                           &ostate_f[0],
                           &ostate_l[0]
                          )

    return np.asarray(ostate_l).reshape((dim_obs_l), order='F')

def omi_init_obs_l_cb (int domain_p,
                       int step,
                       int dim_obs_l,
                       double[::1] observation_l
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obs_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain index
    step : int
        current time step
    dim_obs_l : int
        local dimension of observation vector
    observation_l : ndarray[float]
        local observation vector; Dimension: dim_obs_l

    Returns
    -------
    observation_l : ndarray[float]
        local observation vector; Dimension: dim_obs_l
    """
    c__pdafomi_init_obs_l_cb (&domain_p,
                              &step,
                              &dim_obs_l,
                              &observation_l[0]
                             )

    return np.asarray(observation_l).reshape((dim_obs_l), order='F')

def omi_init_obsvar_l_cb (int domain_p,
                          int step,
                          int dim_obs_l,
                          double[::1] obs_l,
                          double meanvar_l
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obsvar_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        local dimension of observation vector
    obs_l : ndarray[float]
        local observation vector; Dimension: dim_obs_l
    meanvar_l : float
        mean local observation error variance

    Returns
    -------
    meanvar_l : float
        mean local observation error variance
    """
    c__pdafomi_init_obsvar_l_cb (&domain_p,
                                 &step,
                                 &dim_obs_l,
                                 &obs_l[0],
                                 &meanvar_l
                                )

    return meanvar_l

def omi_prodrinva_l_cb (int domain_p,
                        int step,
                        int dim_obs_l,
                        int rank,
                        double[::1] obs_l,
                        double[::1,:] a_l,
                        double[::1,:] c_l
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_l_cb or PDAF source files 

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
        local vector of observations; Dimension: dim_obs_l
    a_l : ndarray[float]
        input matrix; Dimension: dim_obs_l,rank
    c_l : ndarray[float]
        output matrix; Dimension: dim_obs_l,rank

    Returns
    -------
    a_l : ndarray[float]
        input matrix; Dimension: dim_obs_l,rank
    c_l : ndarray[float]
        output matrix; Dimension: dim_obs_l,rank
    """
    cdef double[::1] a_l_f = np.asfortranarray(a_l).ravel(order="F")
    cdef double[::1] c_l_f = np.asfortranarray(c_l).ravel(order="F")
    c__pdafomi_prodrinva_l_cb (&domain_p,
                               &step,
                               &dim_obs_l,
                               &rank,
                               &obs_l[0],
                               &a_l_f[0],
                               &c_l_f[0]
                              )

    return np.asarray(a_l).reshape((dim_obs_l, rank), order='F'), np.asarray(c_l).reshape((dim_obs_l, rank), order='F')

def omi_likelihood_l_cb (int domain_p,
                         int step,
                         int dim_obs_l,
                         double[::1] obs_l,
                         double[::1] resid_l,
                         double lhood_l
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        pe-local dimension of obs. vector
    obs_l : ndarray[float]
        pe-local vector of observations; Dimension: dim_obs_l
    resid_l : ndarray[float]
        input vector of residuum; Dimension: dim_obs_l
    lhood_l : float
        output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[float]
        input vector of residuum; Dimension: dim_obs_l
    lhood_l : float
        output vector - log likelihood
    """
    c__pdafomi_likelihood_l_cb (&domain_p,
                                &step,
                                &dim_obs_l,
                                &obs_l[0],
                                &resid_l[0],
                                &lhood_l
                               )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), lhood_l

def omi_prodrinva_cb (int step,
                      int dim_obs_p,
                      int ncol,
                      double[::1] obs_p,
                      double[::1,:] a_p,
                      double[::1,:] c_p
                     ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of pe-local observation vector
    ncol : int
        number of columns in a_p and c_p
    obs_p : ndarray[float]
        pe-local vector of observations; Dimension: dim_obs_p
    a_p : ndarray[float]
        input matrix; Dimension: dim_obs_p,ncol
    c_p : ndarray[float]
        output matrix; Dimension: dim_obs_p,ncol

    Returns
    -------
    c_p : ndarray[float]
        output matrix; Dimension: dim_obs_p,ncol
    """
    cdef double[::1] a_p_f = np.asfortranarray(a_p).ravel(order="F")
    cdef double[::1] c_p_f = np.asfortranarray(c_p).ravel(order="F")
    c__pdafomi_prodrinva_cb (&step,
                             &dim_obs_p,
                             &ncol,
                             &obs_p[0],
                             &a_p_f[0],
                             &c_p_f[0]
                            )

    return np.asarray(c_p).reshape((dim_obs_p, ncol), order='F')

def omi_likelihood_cb (int step,
                       int dim_obs,
                       double[::1] obs,
                       double[::1] resid,
                       double lhood
                      ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        pe-local dimension of obs. vector
    obs : ndarray[float]
        pe-local vector of observations; Dimension: dim_obs
    resid : ndarray[float]
        input vector of residuum; Dimension: dim_obs
    lhood : float
        output vector - log likelihood

    Returns
    -------
    lhood : float
        output vector - log likelihood
    """
    c__pdafomi_likelihood_cb (&step,
                              &dim_obs,
                              &obs[0],
                              &resid[0],
                              &lhood
                             )

    return lhood

def omi_add_obs_error_cb (int step,
                          int dim_obs_p,
                          double[::1,:] c_p
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_add_obs_error_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_p : int
        dimension of pe-local observation vector
    c_p : ndarray[float]
        matrix to which r is added; Dimension: dim_obs_p,dim_obs_p

    Returns
    -------
    c_p : ndarray[float]
        matrix to which r is added; Dimension: dim_obs_p,dim_obs_p
    """
    cdef double[::1] c_p_f = np.asfortranarray(c_p).ravel(order="F")
    c__pdafomi_add_obs_error_cb (&step,
                                 &dim_obs_p,
                                 &c_p_f[0]
                                )

    return np.asarray(c_p).reshape((dim_obs_p, dim_obs_p), order='F')

def omi_init_obscovar_cb (int step,
                          int dim_obs,
                          int dim_obs_p,
                          double[::1,:] covar,
                          double[::1] m_state_p,
                          bint isdiag
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obscovar_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    dim_obs_p : int
        pe-local dimension of obs. vector
    covar : ndarray[float]
        observation error covar. matrix; Dimension: dim_obs,dim_obs
    m_state_p : ndarray[float]
        observation vector; Dimension: dim_obs_p
    isdiag : bool
        whether matrix r is diagonal

    Returns
    -------
    covar : ndarray[float]
        observation error covar. matrix; Dimension: dim_obs,dim_obs
    isdiag : bool
        whether matrix r is diagonal
    """
    cdef double[::1] covar_f = np.asfortranarray(covar).ravel(order="F")
    c__pdafomi_init_obscovar_cb (&step,
                                 &dim_obs,
                                 &dim_obs_p,
                                 &covar_f[0],
                                 &m_state_p[0],
                                 &isdiag
                                )

    return np.asarray(covar).reshape((dim_obs, dim_obs), order='F'), isdiag

def omi_init_obserr_f_cb (int step,
                          int dim_obs_f,
                          double[::1] obs_f,
                          double[::1] obserr_f
                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_obserr_f_cb or PDAF source files 

    Parameters
    ----------
    step : int
        current time step
    dim_obs_f : int
        full dimension of observation vector
    obs_f : ndarray[float]
        full observation vector; Dimension: dim_obs_f
    obserr_f : ndarray[float]
        full observation error stddev; Dimension: dim_obs_f

    Returns
    -------
    obserr_f : ndarray[float]
        full observation error stddev; Dimension: dim_obs_f
    """
    c__pdafomi_init_obserr_f_cb (&step,
                                 &dim_obs_f,
                                 &obs_f[0],
                                 &obserr_f[0]
                                )

    return np.asarray(obserr_f).reshape((dim_obs_f), order='F')

def omi_prodrinva_hyb_l_cb (int domain_p,
                            int step,
                            int dim_obs_l,
                            int rank,
                            double[::1] obs_l,
                            double alpha,
                            double[::1,:] a_l,
                            double[::1,:] c_l
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_prodRinvA_hyb_l_cb or PDAF source files 

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
        local vector of observations; Dimension: dim_obs_l
    alpha : float
        hybrid weight
    a_l : ndarray[float]
        input matrix; Dimension: dim_obs_l,rank
    c_l : ndarray[float]
        output matrix; Dimension: dim_obs_l,rank

    Returns
    -------
    a_l : ndarray[float]
        input matrix; Dimension: dim_obs_l,rank
    c_l : ndarray[float]
        output matrix; Dimension: dim_obs_l,rank
    """
    cdef double[::1] a_l_f = np.asfortranarray(a_l).ravel(order="F")
    cdef double[::1] c_l_f = np.asfortranarray(c_l).ravel(order="F")
    c__pdafomi_prodrinva_hyb_l_cb (&domain_p,
                                   &step,
                                   &dim_obs_l,
                                   &rank,
                                   &obs_l[0],
                                   &alpha,
                                   &a_l_f[0],
                                   &c_l_f[0]
                                  )

    return np.asarray(a_l).reshape((dim_obs_l, rank), order='F'), np.asarray(c_l).reshape((dim_obs_l, rank), order='F')

def omi_likelihood_hyb_l_cb (int domain_p,
                             int step,
                             int dim_obs_l,
                             double[::1] obs_l,
                             double[::1] resid_l,
                             double alpha,
                             double lhood_l
                            ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_likelihood_hyb_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        current local analysis domain
    step : int
        current time step
    dim_obs_l : int
        pe-local dimension of obs. vector
    obs_l : ndarray[float]
        pe-local vector of observations; Dimension: dim_obs_l
    resid_l : ndarray[float]
        input vector of residuum; Dimension: dim_obs_l
    alpha : float
        hybrid weight
    lhood_l : float
        output vector - log likelihood

    Returns
    -------
    resid_l : ndarray[float]
        input vector of residuum; Dimension: dim_obs_l
    lhood_l : float
        output vector - log likelihood
    """
    c__pdafomi_likelihood_hyb_l_cb (&domain_p,
                                    &step,
                                    &dim_obs_l,
                                    &obs_l[0],
                                    &resid_l[0],
                                    &alpha,
                                    &lhood_l
                                   )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), lhood_l

def omi_obsstats_l (int screen
                   ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obsstats_l or PDAF source files 

    Parameters
    ----------
    screen : int
        < verbosity flag
    """

    c__pdafomi_obsstats_l (&screen
                          )

def omi_weights_l (int verbose,
                   int locweight,
                   double[::1] cradius,
                   double[::1] sradius,
                   double[:,:] mata,
                   double[::1] ivar_obs_l,
                   double[::1] dist_l
                  ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_weights_l or PDAF source files 

    Parameters
    ----------
    verbose : int
        < verbosity flag
    locweight : int
        < localization weight type
    cradius : ndarray[float]
        < localization cut-off radius; Dimension: nobs_l
    sradius : ndarray[float]
        < support radius for weight functions; Dimension: nobs_l
    mata : ndarray[float]
        <; Dimension: nobs_l,ncols
    ivar_obs_l : ndarray[float]
        < local vector of inverse obs. variances (nobs_l); Dimension: nobs_l
    dist_l : ndarray[float]
        < local vector of obs. distances (nobs_l); Dimension: nobs_l

    Returns
    -------
    weight_l : ndarray[float]
        < output: vector of weights; Dimension: nobs_l
    """
    cdef double[::1] mata_f = np.asfortranarray(mata).ravel(order="F")
    cdef int nobs_l, ncols
    nobs_l = mata.shape[0]
    ncols = mata.shape[1]


    cdef double [::1] weight_l = np.zeros((nobs_l), dtype=np.float64).ravel()

    c__pdafomi_weights_l (&verbose,
                          &nobs_l,
                          &ncols,
                          &locweight,
                          &cradius[0],
                          &sradius[0],
                          &mata_f[0],
                          &ivar_obs_l[0],
                          &dist_l[0],
                          &weight_l[0]
                         )

    return np.asarray(weight_l).reshape((nobs_l), order='F')

def omi_weights_l_sgnl (int verbose,
                        int locweight,
                        double cradius,
                        double sradius,
                        double[:,:] mata,
                        double[::1] ivar_obs_l,
                        double[::1] dist_l
                       ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_weights_l_sgnl or PDAF source files 

    Parameters
    ----------
    verbose : int
        < verbosity flag
    locweight : int
        < localization weight type
    cradius : float
        < localization cut-off radius
    sradius : float
        < support radius for weight functions
    mata : ndarray[float]
        <; Dimension: nobs_l,ncols
    ivar_obs_l : ndarray[float]
        < local vector of inverse obs. variances (nobs_l); Dimension: nobs_l
    dist_l : ndarray[float]
        < local vector of obs. distances (nobs_l); Dimension: nobs_l

    Returns
    -------
    weight_l : ndarray[float]
        < output: vector of weights; Dimension: nobs_l
    """
    cdef double[::1] mata_f = np.asfortranarray(mata).ravel(order="F")
    cdef int nobs_l, ncols
    nobs_l = mata.shape[0]
    ncols = mata.shape[1]


    cdef double [::1] weight_l = np.zeros((nobs_l), dtype=np.float64).ravel()

    c__pdafomi_weights_l_sgnl (&verbose,
                               &nobs_l,
                               &ncols,
                               &locweight,
                               &cradius,
                               &sradius,
                               &mata_f[0],
                               &ivar_obs_l[0],
                               &dist_l[0],
                               &weight_l[0]
                              )

    return np.asarray(weight_l).reshape((nobs_l), order='F')

def omi_check_error (int flag
                    ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_check_error or PDAF source files 

    Parameters
    ----------
    flag : int
        < error flag

    Returns
    -------
    flag : int
        < error flag
    """

    c__pdafomi_check_error (&flag
                           )

    return flag

def omi_gather_obsdims ():
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_gather_obsdims or PDAF source files 

    """
    c__pdafomi_gather_obsdims ()

def omi_obsstats (int screen
                 ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_obsstats or PDAF source files 

    Parameters
    ----------
    screen : int
        < verbosity flag
    """

    c__pdafomi_obsstats (&screen
                        )

def omi_init_dim_obs_l_iso (int i_obs,
                            double[::1] coords_l,
                            int locweight,
                            double cradius,
                            double sradius,
                            int cnt_obs_l
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_iso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < index of observation type
    coords_l : ndarray[float]
        < coordinates of current analysis domain; Dimension: ncoord
    locweight : int
        < type of localization function
    cradius : float
        < localization cut-off radius (single or vector)
    sradius : float
        < support radius of localization function (single or vector)
    cnt_obs_l : int
        < local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        < local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_iso (&i_obs,
                                   &ncoord,
                                   &coords_l[0],
                                   &locweight,
                                   &cradius,
                                   &sradius,
                                   &cnt_obs_l
                                  )

    return cnt_obs_l

def omi_init_dim_obs_l_noniso (int i_obs,
                               double[::1] coords_l,
                               int locweight,
                               double[::1] cradius,
                               double[::1] sradius,
                               int cnt_obs_l
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_noniso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < index of observation type
    coords_l : ndarray[float]
        < coordinates of current analysis domain; Dimension: ncoord
    locweight : int
        < type of localization function
    cradius : ndarray[float]
        < vector of localization cut-off radii; Dimension: ncoord
    sradius : ndarray[float]
        < vector of support radii of localization function; Dimension: ncoord
    cnt_obs_l : int
        < local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        < local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_noniso (&i_obs,
                                      &ncoord,
                                      &coords_l[0],
                                      &locweight,
                                      &cradius[0],
                                      &sradius[0],
                                      &cnt_obs_l
                                     )

    return cnt_obs_l

def omi_init_dim_obs_l_noniso_locweights (int i_obs,
                                          double[::1] coords_l,
                                          int[::1] locweights,
                                          double[::1] cradius,
                                          double[::1] sradius,
                                          int cnt_obs_l
                                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l_noniso_locweights or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < index of observation type
    coords_l : ndarray[float]
        < coordinates of current analysis domain; Dimension: ncoord
    locweights : ndarray[int]
        < types of localization function; Dimension: 2
    cradius : ndarray[float]
        < vector of localization cut-off radii; Dimension: ncoord
    sradius : ndarray[float]
        < vector of support radii of localization function; Dimension: ncoord
    cnt_obs_l : int
        < local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        < local dimension of current observation vector
    """
    cdef int ncoord
    ncoord = coords_l.shape[0]


    c__pdafomi_init_dim_obs_l_noniso_locweights (&i_obs,
                                                 &ncoord,
                                                 &coords_l[0],
                                                 &locweights[0],
                                                 &cradius[0],
                                                 &sradius[0],
                                                 &cnt_obs_l
                                                )

    return cnt_obs_l

def omi_localize_covar_iso (int i_obs,
                            int locweight,
                            double cradius,
                            double sradius,
                            double[:,:] coords,
                            double[:,:] hp,
                            double[:,:] hph
                           ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_iso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < index of observation type
    locweight : int
        < localization weight type
    cradius : float
        < localization radius
    sradius : float
        < support radius for weight functions
    coords : ndarray[float]
        < coordinates of state vector elements; Dimension: ncoord,dim_p
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs

    Returns
    -------
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] hp_f = np.asfortranarray(hp).ravel(order="F")
    cdef double[::1] hph_f = np.asfortranarray(hph).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = hp.shape[0]
    _ = hp.shape[1]


    c__pdafomi_localize_covar_iso (&i_obs,
                                   &dim_p,
                                   &dim_obs,
                                   &ncoord,
                                   &locweight,
                                   &cradius,
                                   &sradius,
                                   &coords_f[0],
                                   &hp_f[0],
                                   &hph_f[0]
                                  )

    return np.asarray(hp).reshape((dim_obs, dim_p), order='F'), np.asarray(hph).reshape((dim_obs, dim_obs), order='F')

def omi_localize_covar_noniso (int i_obs,
                               int locweight,
                               double[::1] cradius,
                               double[::1] sradius,
                               double[:,:] coords,
                               double[:,:] hp,
                               double[:,:] hph
                              ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_noniso or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < data type with full observation
    locweight : int
        < localization weight type
    cradius : ndarray[float]
        < vector of localization cut-off radii; Dimension: ncoord
    sradius : ndarray[float]
        < vector of support radii of localization function; Dimension: ncoord
    coords : ndarray[float]
        < coordinates of state vector elements; Dimension: ncoord,dim_p
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs

    Returns
    -------
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] hp_f = np.asfortranarray(hp).ravel(order="F")
    cdef double[::1] hph_f = np.asfortranarray(hph).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = hp.shape[0]
    _ = hp.shape[1]


    c__pdafomi_localize_covar_noniso (&i_obs,
                                      &dim_p,
                                      &dim_obs,
                                      &ncoord,
                                      &locweight,
                                      &cradius[0],
                                      &sradius[0],
                                      &coords_f[0],
                                      &hp_f[0],
                                      &hph_f[0]
                                     )

    return np.asarray(hp).reshape((dim_obs, dim_p), order='F'), np.asarray(hph).reshape((dim_obs, dim_obs), order='F')

def omi_localize_covar_noniso_locweights (int i_obs,
                                          int[::1] locweights,
                                          double[::1] cradius,
                                          double[::1] sradius,
                                          double[:,:] coords,
                                          double[:,:] hp,
                                          double[:,:] hph
                                         ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar_noniso_locweights or PDAF source files 

    Parameters
    ----------
    i_obs : int
        < index of observation type
    locweights : ndarray[int]
        < types of localization function; Dimension: 2
    cradius : ndarray[float]
        < vector of localization cut-off radii; Dimension: ncoord
    sradius : ndarray[float]
        < vector of support radii of localization function; Dimension: ncoord
    coords : ndarray[float]
        < coordinates of state vector elements; Dimension: ncoord,dim_p
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs

    Returns
    -------
    hp : ndarray[float]
        < matrix hp, dimension (nobs, dim); Dimension: dim_obs,dim_p
    hph : ndarray[float]
        < matrix hph, dimension (nobs, nobs); Dimension: dim_obs,dim_obs
    """
    cdef double[::1] coords_f = np.asfortranarray(coords).ravel(order="F")
    cdef double[::1] hp_f = np.asfortranarray(hp).ravel(order="F")
    cdef double[::1] hph_f = np.asfortranarray(hph).ravel(order="F")
    cdef int ncoord, dim_p, dim_obs
    ncoord = coords.shape[0]
    dim_p = coords.shape[1]
    dim_obs = hp.shape[0]
    _ = hp.shape[1]


    c__pdafomi_localize_covar_noniso_locweights (&i_obs,
                                                 &dim_p,
                                                 &dim_obs,
                                                 &ncoord,
                                                 &locweights[0],
                                                 &cradius[0],
                                                 &sradius[0],
                                                 &coords_f[0],
                                                 &hp_f[0],
                                                 &hph_f[0]
                                                )

    return np.asarray(hp).reshape((dim_obs, dim_p), order='F'), np.asarray(hph).reshape((dim_obs, dim_obs), order='F')

def omi_omit_by_inno_l_cb (int domain_p,
                           int dim_obs_l,
                           double[::1] resid_l,
                           double[::1] obs_l
                          ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_omit_by_inno_l_cb or PDAF source files 

    Parameters
    ----------
    domain_p : int
        < current local analysis domain
    dim_obs_l : int
        < pe-local dimension of obs. vector
    resid_l : ndarray[float]
        < input vector of residuum; Dimension: dim_obs_l
    obs_l : ndarray[float]
        < input vector of local observations; Dimension: dim_obs_l

    Returns
    -------
    resid_l : ndarray[float]
        < input vector of residuum; Dimension: dim_obs_l
    obs_l : ndarray[float]
        < input vector of local observations; Dimension: dim_obs_l
    """
    c__pdafomi_omit_by_inno_l_cb (&domain_p,
                                  &dim_obs_l,
                                  &resid_l[0],
                                  &obs_l[0]
                                 )

    return np.asarray(resid_l).reshape((dim_obs_l), order='F'), np.asarray(obs_l).reshape((dim_obs_l), order='F')

def omi_omit_by_inno_cb (int dim_obs_f,
                         double[::1] resid_f,
                         double[::1] obs_f
                        ):
    """See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/PDAFomi_omit_by_inno_cb or PDAF source files 

    Parameters
    ----------
    dim_obs_f : int
        < full dimension of obs. vector
    resid_f : ndarray[float]
        < input vector of residuum; Dimension: dim_obs_f
    obs_f : ndarray[float]
        < input vector of full observations; Dimension: dim_obs_f

    Returns
    -------
    resid_f : ndarray[float]
        < input vector of residuum; Dimension: dim_obs_f
    obs_f : ndarray[float]
        < input vector of full observations; Dimension: dim_obs_f
    """
    c__pdafomi_omit_by_inno_cb (&dim_obs_f,
                                &resid_f[0],
                                &obs_f[0]
                               )

    return np.asarray(resid_f).reshape((dim_obs_f), order='F'), np.asarray(obs_f).reshape((dim_obs_f), order='F')

