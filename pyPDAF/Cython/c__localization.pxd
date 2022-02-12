"""Declaration file of pyPDAF.

The declaration file declares the user-defined PDAF subroutines
for localization in C language format. 
"""
cdef void c__init_dim_l_pdaf(int* step, int* domain_p, int* dim_l);

cdef void c__init_n_domains_pdaf(int* step, int* n_domains_p);

cdef void c__g2l_state_pdaf(int* step, int* domain_p, int* dim_p, 
                            double* state_p, int* dim_l, 
                            double* state_l);

cdef void c__l2g_state_pdaf(int* step, int* domain_p, int* dim_l,
                            double* state_l, int* dim_p, 
                            double* state_p);