from pyPDAF.cfi_binding cimport CFI_cdesc_t

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


