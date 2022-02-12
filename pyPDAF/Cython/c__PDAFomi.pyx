"""Implementation file for the user-defined PDAFomi routines

The main goal of this file is to convert PDAF given Fortran
array/pointers to numpy arrays
"""
import numpy as np


def py__init_dim_obs_PDAFomi(step, dim_obs):
    """Summary
    
    Parameters
    ----------
    step : TYPE
        Description
    dim_obs : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
    """
    raise RuntimeError('...Wrong init_dim_obs_PDAFomi is called!!!...')

def py__init_dim_obs_l_PDAFomi(domain_p, step, dim_obs, dim_obs_l):
    """Summary
    
    Parameters
    ----------
    domain_p : TYPE
        Description
    step : TYPE
        Description
    dim_obs : TYPE
        Description
    dim_obs_l : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
    """
    raise RuntimeError('...Wrong init_dim_obs_l_PDAFomi is called!!!...')

def py__obs_op_PDAFomi(step, state_p, ostate):
    """Summary
    
    Parameters
    ----------
    step : TYPE
        Description
    state_p : TYPE
        Description
    ostate : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
    """
    raise RuntimeError('...Wrong obs_op_PDAFomi is called!!!...')

def py__localize_covar_PDAFomi(hp_p_numpy, hph_numpy):
    """Summary
    
    Parameters
    ----------
    hp_p_numpy : TYPE
        Description
    hph_numpy : TYPE
        Description
    
    Raises
    ------
    RuntimeError
        Description
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
