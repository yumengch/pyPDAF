"""Declaration file of pyPDAF.

The declaration file declares the user-defined PDAFomi subroutines
in C language format. 
"""
cdef void c__init_dim_obs_PDAFomi(int* step, int* dim_obs);
cdef void c__obs_op_PDAFomi(int* step, int* dim_p, int* dim_obs, 
                            double* state_p, double* ostate);
cdef void c__init_dim_obs_l_PDAFomi(int* domain_p, int* step,
                                    int* dim_obs, int* dim_obs_l);
cdef void c__localize_covar_PDAFomi(int* dim_p, int* dim_obs, 
                                    double* hp_p, double* hph);