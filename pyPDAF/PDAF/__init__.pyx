import pyPDAF.Cython as PDAFcython
import pyPDAF.Cython.c__localization as localization
import pyPDAF.Cython.c__PDAFomi as cPDAFomi

cimport pyPDAF.Cython as c__PDAFcython
cimport pyPDAF.Cython.c__localization as c__localization
cimport pyPDAF.Cython.c__PDAFomi as c__PDAFomi

import numpy as np

def init(int filtertype, int subtype, 
         filter_param_i,  filter_param_r, 
         int COMM_model, int COMM_filter, 
         int COMM_couple, int task_id, 
         int n_modeltasks, bint filterpe, 
         py__init_ens_pdaf, int screen):
    cdef int[::1] filter_param_i_view = np.array(
                                    filter_param_i, 
                                    dtype=np.intc)
    cdef int dim_pint = len(filter_param_i)
    cdef double[::1] filter_param_r_view = filter_param_r
    cdef int dim_preal = len(filter_param_r)

    cdef int status_pdaf

    PDAFcython.py__init_ens_pdaf = py__init_ens_pdaf

    c__pdaf_init(&filtertype, &subtype,
        &filter_param_i_view[0], &dim_pint,
        &filter_param_r_view[0], &dim_preal,
        &COMM_model, &COMM_filter,
        &COMM_couple,
        &task_id, &n_modeltasks,
        &filterpe, c__PDAFcython.c__init_ens_pdaf,
        &screen,
        &status_pdaf)

    return status_pdaf

def get_state(py__next_observation_pdaf,
                    py__distribute_state_pdaf,
                    py__prepoststep_ens_pdaf):
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf

    cdef int steps, doexit, status_pdaf
    cdef double timenow

    c__pdaf_get_state(
                &steps, &timenow, &doexit, 
                c__PDAFcython.c__next_observation_pdaf,
                c__PDAFcython.c__distribute_state_pdaf,
                c__PDAFcython.c__prepoststep_ens_pdaf,
                &status_pdaf)
    return steps, timenow, doexit, status_pdaf

def get_localfilter():
    cdef int localfilter
    c__pdaf_get_localfilter(&localfilter)
    return localfilter

def PDAFomi_assimilate_local(py__collect_state_pdaf,
                             py__distribute_state_pdaf,
                             py__init_dim_obs_PDAFomi,
                             py__obs_op_PDAFomi,
                             py__prepoststep_ens_pdaf,
                             py__init_n_domains_pdaf,
                             py__init_dim_l_pdaf,
                             py__init_dim_obs_l_PDAFomi,
                             py__g2l_state_pdaf,
                             py__l2g_state_pdaf,
                             py__next_observation_pdaf):
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    cPDAFomi.py__init_dim_obs_PDAFomi = py__init_dim_obs_PDAFomi
    cPDAFomi.py__obs_op_PDAFomi = py__obs_op_PDAFomi
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf
    localization.py__init_n_domains_pdaf = py__init_n_domains_pdaf
    localization.py__init_dim_l_pdaf = py__init_dim_l_pdaf
    cPDAFomi.py__init_dim_obs_l_PDAFomi = py__init_dim_obs_l_PDAFomi
    localization.py__g2l_state_pdaf = py__g2l_state_pdaf
    localization.py__l2g_state_pdaf = py__l2g_state_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int status_pdaf
    c__pdafomi_assimilate_local(
             c__PDAFcython.c__collect_state_pdaf,
             c__PDAFcython.c__distribute_state_pdaf,
             c__PDAFomi.c__init_dim_obs_PDAFomi,
             c__PDAFomi.c__obs_op_PDAFomi,
             c__PDAFcython.c__prepoststep_ens_pdaf,
             c__localization.c__init_n_domains_pdaf,
             c__localization.c__init_dim_l_pdaf,
             c__PDAFomi.c__init_dim_obs_l_PDAFomi,
             c__localization.c__g2l_state_pdaf,
             c__localization.c__l2g_state_pdaf,
             c__PDAFcython.c__next_observation_pdaf, &status_pdaf)
    return status_pdaf

def PDAFomi_assimilate_global(py__collect_state_pdaf,
                              py__distribute_state_pdaf,
                              py__init_dim_obs_PDAFomi,
                              py__obs_op_PDAFomi,
                              py__prepoststep_ens_pdaf,
                              py__next_observation_pdaf):
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    cPDAFomi.py__init_dim_obs_PDAFomi = py__init_dim_obs_PDAFomi
    cPDAFomi.py__obs_op_PDAFomi = py__obs_op_PDAFomi
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int status_pdaf
    c__pdafomi_assimilate_global(
        c__PDAFcython.c__collect_state_pdaf,
        c__PDAFcython.c__distribute_state_pdaf,
        c__PDAFomi.c__init_dim_obs_PDAFomi,
        c__PDAFomi.c__obs_op_PDAFomi,
        c__PDAFcython.c__prepoststep_ens_pdaf,
        c__PDAFcython.c__next_observation_pdaf, &status_pdaf)
    return status_pdaf

def PDAFomi_assimilate_lenkf(py__collect_state_pdaf,
                             py__distribute_state_pdaf,
                             py__init_dim_obs_PDAFomi,
                             py__obs_op_PDAFomi,
                             py__prepoststep_ens_pdaf,
                             py__localize_covar_PDAFomi,
                             py__next_observation_pdaf):
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    cPDAFomi.py__init_dim_obs_PDAFomi = py__init_dim_obs_PDAFomi
    cPDAFomi.py__obs_op_PDAFomi = py__obs_op_PDAFomi
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf
    cPDAFomi.py__localize_covar_PDAFomi = py__localize_covar_PDAFomi
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
    cdef int status_pdaf
    c__pdafomi_assimilate_lenkf(
        c__PDAFcython.c__collect_state_pdaf,
        c__PDAFcython.c__distribute_state_pdaf,
        c__PDAFomi.c__init_dim_obs_PDAFomi,
        c__PDAFomi.c__obs_op_PDAFomi,
        c__PDAFcython.c__prepoststep_ens_pdaf,
        c__PDAFomi.c__localize_covar_PDAFomi,
        c__PDAFcython.c__next_observation_pdaf, &status_pdaf)
    return status_pdaf