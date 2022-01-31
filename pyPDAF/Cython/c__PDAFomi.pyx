import numpy as np


def py__init_dim_obs_PDAFomi(step, dim_obs):
    raise RuntimeError('...Wrong init_dim_obs_PDAFomi is called!!!...')

def py__init_dim_obs_l_PDAFomi(domain_p, step, dim_obs, dim_obs_l):
    raise RuntimeError('...Wrong init_dim_obs_l_PDAFomi is called!!!...')

def py__obs_op_PDAFomi(step, state_p, ostate):
    raise RuntimeError('...Wrong obs_op_PDAFomi is called!!!...')

def py__localize_covar_PDAFomi(hp_p_numpy, hph_numpy):
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
