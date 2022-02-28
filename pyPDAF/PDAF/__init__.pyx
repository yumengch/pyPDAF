"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and National Centre for Earth Observation

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This module calls PDAF iso_C_binding wrappers.

Function calls convert Python objects for Fortran routines.
"""
import pyPDAF.UserFunc as PDAFcython
import pyPDAF.UserFunc.c__localization as localization
import pyPDAF.UserFunc.c__PDAFomi as cPDAFomi

cimport pyPDAF.UserFunc as c__PDAFcython
cimport pyPDAF.UserFunc.c__localization as c__localization
cimport pyPDAF.UserFunc.c__PDAFomi as c__PDAFomi

import numpy as np


def init(int filtertype, int subtype, 
         filter_param_i,  filter_param_r, 
         int COMM_model, int COMM_filter, 
         int COMM_couple, int task_id, 
         int n_modeltasks, bint filterpe, 
         py__init_ens_pdaf, int screen):
    """Python wrapper for PDAFinit in PDAF.

        Parameters
        ----------
        filtertype : int
            Filter index
        subtype : int
            Subtype of selected filter
        filter_param_i : ndarray
            an array of integer filter parameters
        filter_param_r : ndarray
            an array of float filter parameters
        COMM_model : int
            Fortran communicator for model
        COMM_filter : int
            Fortran communicator for filter
        n_modeltasks : int
            number of model tasks
        filterpe : bool
            Whether the current PE is filter PE
        py__init_ens_pdaf : func
            user-supplied functions
        screen : int
            verbosity of PDAF screen output


        Returns
        -------
        status_pdaf : int
            PDAF status code; 0 means pass
    """
    cdef int[::1] filter_param_i_view = np.array(filter_param_i, 
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
    """Python wrapper for PDAF_get_state.

        Parameters
        ----------
        py__next_observation_pdaf : func
            user-supplied Python functions
        py__distribute_state_pdaf : func
            user-supplied Python functions
        py__prepoststep_ens_pdaf : func
            user-supplied Python functions


        Returns
        -------
        steps : int
            flag and number of time steps
        timenow : int
            current time step
        doexit : int
            Whether to exit from forecasts. Exit if 1
        status_pdaf : int
            PDAF status code; 0 means pass
    """
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf

    cdef int steps, doexit, status_pdaf
    cdef double timenow

    c__pdaf_get_state(&steps, &timenow, &doexit, 
                      c__PDAFcython.c__next_observation_pdaf,
                      c__PDAFcython.c__distribute_state_pdaf,
                      c__PDAFcython.c__prepoststep_ens_pdaf,
                      &status_pdaf)

    return steps, timenow, doexit, status_pdaf

def get_localfilter():
    """Python wrapper for PDAF_get_localfilter.

        Returns
        -------
        localfilter : int
            Whether PDAF uses a local filter
    """
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
    """Python wrapper for PDAFomi_assimilate_local.

        Parameters
        ----------
        py__collect_state_pdaf : func
            user-supplied functions
        py__distribute_state_pdaf : func
            user-supplied functions
        py__init_dim_obs_PDAFomi : func
            user-supplied functions
        py__obs_op_PDAFomi : func
            user-supplied functions
        py__prepoststep_ens_pdaf : func
            user-supplied functions
        py__init_n_domains_pdaf : func
            user-supplied functions
        py__init_dim_l_pdaf : func
            user-supplied functions
        py__init_dim_obs_l_PDAFom : func
            user-supplied functions
        py__g2l_state_pdaf : func
            user-supplied functions
        py__l2g_state_pdaf : func
            user-supplied functions
        py__next_observation_pdaf : func
            user-supplied functions

        Returns
        -------
        status_pdaf : int
            PDAF status code; 0 means pass
    """
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
    c__pdafomi_assimilate_local(c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__distribute_state_pdaf,
                                c__PDAFomi.c__init_dim_obs_PDAFomi,
                                c__PDAFomi.c__obs_op_PDAFomi,
                                c__PDAFcython.c__prepoststep_ens_pdaf,
                                c__localization.c__init_n_domains_pdaf,
                                c__localization.c__init_dim_l_pdaf,
                                c__PDAFomi.c__init_dim_obs_l_PDAFomi,
                                c__localization.c__g2l_state_pdaf,
                                c__localization.c__l2g_state_pdaf,
                                c__PDAFcython.c__next_observation_pdaf,
                                &status_pdaf)
    return status_pdaf

def PDAFomi_assimilate_global(py__collect_state_pdaf,
                              py__distribute_state_pdaf,
                              py__init_dim_obs_PDAFomi,
                              py__obs_op_PDAFomi,
                              py__prepoststep_ens_pdaf,
                              py__next_observation_pdaf):
    """Python wrapper for PDAFomi_assimilate_global.

        Parameters
        ----------
        py__collect_state_pdaf : func
            user-supplied functions
        py__distribute_state_pdaf : func
            user-supplied functions
        py__init_dim_obs_PDAFomi : func
            user-supplied functions
        py__obs_op_PDAFomi : func
            user-supplied functions
        py__prepoststep_ens_pdaf : func
            user-supplied functions
        py__next_observation_pdaf : func
            user-supplied functions

        Returns
        -------
        status_pdaf : int
            PDAF status code; 0 means pass
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    cPDAFomi.py__init_dim_obs_PDAFomi = py__init_dim_obs_PDAFomi
    cPDAFomi.py__obs_op_PDAFomi = py__obs_op_PDAFomi
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf

    cdef int status_pdaf
    c__pdafomi_assimilate_global(c__PDAFcython.c__collect_state_pdaf,
                                 c__PDAFcython.c__distribute_state_pdaf,
                                 c__PDAFomi.c__init_dim_obs_PDAFomi,
                                 c__PDAFomi.c__obs_op_PDAFomi,
                                 c__PDAFcython.c__prepoststep_ens_pdaf,
                                 c__PDAFcython.c__next_observation_pdaf,
                                 &status_pdaf)
    return status_pdaf

def PDAFomi_assimilate_lenkf(py__collect_state_pdaf,
                             py__distribute_state_pdaf,
                             py__init_dim_obs_PDAFomi,
                             py__obs_op_PDAFomi,
                             py__prepoststep_ens_pdaf,
                             py__localize_covar_PDAFomi,
                             py__next_observation_pdaf):
    """Python wrapper for PDAFomi_assimilate_local.

        Parameters
        ----------
        py__collect_state_pdaf : func
            user-supplied functions
        py__distribute_state_pdaf : func
            user-supplied functions
        py__init_dim_obs_PDAFomi : func
            user-supplied functions
        py__obs_op_PDAFomi : func
            user-supplied functions
        py__prepoststep_ens_pdaf : func
            user-supplied functions
        py__init_n_domains_pdaf : func
            user-supplied functions
        py__init_dim_l_pdaf : func
            user-supplied functions
        py__init_dim_obs_l_PDAFom : func
            user-supplied functions
        py__g2l_state_pdaf : func
            user-supplied functions
        py__l2g_state_pdaf : func
            user-supplied functions
        py__next_observation_pdaf : func
            user-supplied functions

        Returns
        -------
        status_pdaf : int
            PDAF status code; 0 means pass
    """
    PDAFcython.py__collect_state_pdaf = py__collect_state_pdaf
    PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
    cPDAFomi.py__init_dim_obs_PDAFomi = py__init_dim_obs_PDAFomi
    cPDAFomi.py__obs_op_PDAFomi = py__obs_op_PDAFomi
    PDAFcython.py__prepoststep_ens_pdaf = py__prepoststep_ens_pdaf
    cPDAFomi.py__localize_covar_PDAFomi = py__localize_covar_PDAFomi
    PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
    cdef int status_pdaf
    c__pdafomi_assimilate_lenkf(c__PDAFcython.c__collect_state_pdaf,
                                c__PDAFcython.c__distribute_state_pdaf,
                                c__PDAFomi.c__init_dim_obs_PDAFomi,
                                c__PDAFomi.c__obs_op_PDAFomi,
                                c__PDAFcython.c__prepoststep_ens_pdaf,
                                c__PDAFomi.c__localize_covar_PDAFomi,
                                c__PDAFcython.c__next_observation_pdaf,
                                &status_pdaf)
    return status_pdaf