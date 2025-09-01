from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_iau_init_weights(int* type_iau,
    int* nsteps_iau) noexcept nogil;

cdef extern void c__pdaf_iau_update_inc(
    CFI_cdesc_t* ens_ana, CFI_cdesc_t* state_ana) noexcept nogil;

cdef extern void c__pdaf_iau_add_inc_ens(int* step, int* dim_p,
    int* dim_ens_task, CFI_cdesc_t* ens, CFI_cdesc_t* state,
    void (*c__collect_state_pdaf)(int* , double* ),
    void (*c__distribute_state_pdaf)(int* ,
                                     double* )) noexcept nogil;


cdef extern void c__pdaf_iau_update_ens(
    CFI_cdesc_t* ens, CFI_cdesc_t* state) noexcept nogil;

cdef extern void c__pdaf_iau_dealloc() noexcept nogil;

