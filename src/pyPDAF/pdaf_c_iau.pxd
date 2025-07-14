from .cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_iau_init(int* type_iau_in, int* nsteps_iau_in, 
    int* flag) noexcept nogil;

cdef extern void c__pdaf_iau_reset(int* type_iau_in, int* nsteps_iau_in, 
    int* flag) noexcept nogil;

cdef extern void c__pdaf_iau_set_weights(int* iweights, 
    double* weights) noexcept nogil;

cdef extern void c__pdaf_iau_set_pointer(CFI_cdesc_t* iau_ptr, 
    int* flag) noexcept nogil;

cdef extern void c__pdaf_iau_init_inc(int* dim_p, int* dim_ens_l, 
    double* ens_inc, int* flag) noexcept nogil;

cdef extern void c__pdaf_iau_add_inc(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , 
                                     double* )) noexcept nogil;
    

