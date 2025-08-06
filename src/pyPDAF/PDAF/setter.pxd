from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_set_iparam_filters(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_set_rparam_filters(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_set_comm_pdaf(
    int* in_comm_pdaf) noexcept nogil;

cdef extern void c__pdaf_set_debug_flag(
    int* debugval) noexcept nogil;

cdef extern void c__pdaf_set_ens_pointer(CFI_cdesc_t* ens_ptr,
    int* status) noexcept nogil;

cdef extern void c__pdaf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_set_memberid(
    int* memberid) noexcept nogil;

cdef extern void c__pdaf_set_offline_mode(
    int* screen) noexcept nogil;

cdef extern void c__pdaf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_set_seedset(
    int* seedset_in) noexcept nogil;

cdef extern void c__pdaf_set_smootherens(CFI_cdesc_t* sens_point,
    int* maxlag, int* status) noexcept nogil;

cdef extern void c__pdaf_set_forget(int* step, int* localfilter,
    int* dim_obs_p, int* dim_ens, double* mens_p, double* mstate_p,
    double* obs_p,
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    double* forget_in, double* forget_out,
    int* screen) noexcept nogil;

