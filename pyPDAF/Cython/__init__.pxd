"""Declaration file of pyPDAF.

The declaration file declares the user-defined PDAF subroutines
in C language format. 
"""
cdef void c__init_ens_pdaf(int* filtertype, int* dim_p, int* dim_ens, 
	                       double* state_p, double* uinv, 
                           double* ens_p, int* flag)

cdef void c__distribute_state_pdaf(int* dim_p, double* state_p);

cdef void c__collect_state_pdaf(int* dim_p, double* state_p);

cdef void c__next_observation_pdaf(int* stepnow, 
    int* nsteps, int* doexit, double* time);

cdef void c__prepoststep_ens_pdaf(int* step, int* dim_p, int* dim_ens,
            int* dim_ens_p, int* dim_obs_p, 
            double* state_p, double* uinv, double* ens_p, int* flag);
