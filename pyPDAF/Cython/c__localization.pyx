import numpy as np


def py__init_dim_l_pdaf(step, domain_p, dim_l):
    raise RuntimeError('...Wrong init_dim_l_pdaf is called!!!...')

def py__init_n_domains_pdaf(step, n_domains_p):
    raise RuntimeError('...Wrong init_n_domains_pdaf is called!!!...')

def py__g2l_state_pdaf(step, domain_p, state_p, state_l):
    raise RuntimeError('...Wrong distribute_state_pdaf is called!!!...')

def py__l2g_state_pdaf(step, domain_p, state_l, state_p):
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
