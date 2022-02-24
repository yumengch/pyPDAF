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

Implementation of the user-defined routines 
for localization in PDAF
"""
import numpy as np


def py__init_dim_l_pdaf(step, domain_p, dim_l):
    """Initialize local state dimension
    
    Parameters
    ----------
    step : int
        Current time step
    domain_p : int
        index of current local analysis domain
    dim_l : int
        local state dimension 
    
    Returns
    -------
    dim_l : int
        local state dimension 

    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong init_dim_l_pdaf is called!!!...')


def py__init_n_domains_pdaf(step, n_domains_p):
    """Initialise the number of local domains
    
    Parameters
    ----------
    step : int
        current time step
    n_domains_p : int
        pe-local number of analysis domains  

    Returns
    -------
    n_domains_p : int
        pe-local number of analysis domains  
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong init_n_domains_pdaf is called!!!...')


def py__g2l_state_pdaf(step, domain_p, state_p, state_l):
    """Convert from global state vector to local
    
    Parameters
    ----------
    step : int
        current time step
    domain_p : int
        current local analysis domain
    state_p : ndarray
        pe-local full state vector (shape: (dim_p))
    state_l : ndarray
        state vector on local analysis domain (shape: (dim_l))
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong distribute_state_pdaf is called!!!...')


def py__l2g_state_pdaf(step, domain_p, state_l, state_p):
    """Convert from local state vector to global
    
    Parameters
    ----------
    step : int
        current time step
    domain_p : int
        current local analysis domain
    state_l : ndarray
        state vector on local analysis domain (shape: (dim_l))
    state_p : ndarray
        pe-local full state vector (shape: (dim_p))
    
    Raises
    ------
    RuntimeError
        No user-supplied function
    """
    raise RuntimeError('...Wrong l2g_state_pdaf is called!!!...')


cdef void c__init_dim_l_pdaf(int* step, int* domain_p, int* dim_l):
    dim_l[0] = py__init_dim_l_pdaf(step[0], domain_p[0], dim_l[0])


cdef void c__init_n_domains_pdaf(int* step, int* n_domains_p):
    n_domains_p[0] = py__init_n_domains_pdaf(step[0], n_domains_p[0])


cdef void c__g2l_state_pdaf(int* step, int* domain_p, int* dim_p, 
                            double* state_p, int* dim_l, 
                            double* state_l):
    state_p_numpy = np.asarray(<double[:dim_p[0]]> state_p)
    state_l_numpy = np.asarray(<double[:dim_l[0]]> state_l)

    py__g2l_state_pdaf(step[0], domain_p[0], 
                       state_p_numpy, state_l_numpy)


cdef void c__l2g_state_pdaf(int* step, int* domain_p, int* dim_l,
                            double* state_l, int* dim_p, 
                            double* state_p):
    state_l_numpy = np.asarray(<double[:dim_l[0]]> state_l)
    state_p_numpy = np.asarray(<double[:dim_p[0]]> state_p)
    py__l2g_state_pdaf(step[0], domain_p[0], 
                       state_l_numpy, state_p_numpy)
