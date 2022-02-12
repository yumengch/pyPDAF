"""Implementation file for the user-defined routines 
for localization in PDAF

The main goal of this file is to convert PDAF given Fortran
array/pointers to numpy arrays
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
        Description
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
        Description
    """
    raise RuntimeError('...Wrong init_n_domains_pdaf is called!!!...')

def py__g2l_state_pdaf(step, domain_p, state_p, state_l):
    """Summary
    
    Parameters
    ----------
    step : TYPE
        Description
    domain_p : TYPE
        Description
    state_p : TYPE
        Description
    state_l : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
    """
    raise RuntimeError('...Wrong distribute_state_pdaf is called!!!...')

def py__l2g_state_pdaf(step, domain_p, state_l, state_p):
    """Summary
    
    Parameters
    ----------
    step : TYPE
        Description
    domain_p : TYPE
        Description
    state_l : TYPE
        Description
    state_p : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
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
