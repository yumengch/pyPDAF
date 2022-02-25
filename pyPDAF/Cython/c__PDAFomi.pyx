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

Implementation of the user-defined PDAF subroutines
in Cython syntax. 
"""
import numpy as np


def py__init_dim_obs_PDAFomi(step, dim_obs):
    """default user-supplied init_dim_obs_PDAFOmi
    
    Parameters
    ----------
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong init_dim_obs_PDAFomi is called!!!...')

def py__init_dim_obs_l_PDAFomi(domain_p, step, dim_obs, dim_obs_l):
    """default user-supplied init_dim_obs_l_PDAFomi
    
    Parameters
    ----------
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    dim_obs_l : int
        dimension of local observation vector
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong init_dim_obs_l_PDAFomi is called!!!...')

def py__obs_op_PDAFomi(step, state_p, ostate):
    """default user-supplied obs_op_PDAFomi
    
    Parameters
    ----------
    step : int
        current time step
    state_p : ndarray
        local PE state vector
    ostate : ndarray
        state vector in obs space
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong obs_op_PDAFomi is called!!!...')

def py__localize_covar_PDAFomi(HP_p, HPH):
    """Summary
    
    Parameters
    ----------
    HP_p : ndarray
        matrix HP
    HPH : ndarray
        matrix HPH
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong localize_covar_PDAFomi is called!!!...')


cdef void c__init_dim_obs_PDAFomi(int* step, int* dim_obs):
    dim_obs[0] = py__init_dim_obs_PDAFomi(step[0], dim_obs[0])


cdef void c__init_dim_obs_l_PDAFomi(int* domain_p, int* step,
                                    int* dim_obs, int* dim_obs_l):
    dim_obs_l[0] = py__init_dim_obs_l_PDAFomi(domain_p[0], step[0], 
                               dim_obs[0], dim_obs_l[0]);


cdef void c__obs_op_PDAFomi(int* step, int* dim_p, int* dim_obs, 
                            double* state_p, double* ostate):
    state_p_numpy = np.asarray(<double[:dim_p[0]]> state_p)
    ostate_numpy = np.asarray(<double[:dim_obs[0]]> ostate)
    py__obs_op_PDAFomi(step[0], state_p_numpy, ostate_numpy)


cdef void c__localize_covar_PDAFomi(int* dim_p, int* dim_obs, 
                                    double* hp_p, double* hph):
    if (dim_p[0] != dim_obs[0]):
        hp_p_numpy = np.asarray(<double[:dim_p[0], :dim_obs[0]]> hp_p)
    else:
        hp_p_numpy = np.asarray(<double[:dim_p[0], :dim_p[0]]> hp_p).T
    hph_numpy = np.asarray(<double[:dim_obs[0], :dim_obs[0]]> hph).T
    py__localize_covar_PDAFomi(hp_p_numpy, hph_numpy)
