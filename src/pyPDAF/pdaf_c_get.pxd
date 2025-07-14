from .cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_get_assim_flag(
    int* did_assim) noexcept nogil;

cdef extern void c__pdaf_get_localfilter(
    int* localfilter_out) noexcept nogil;

cdef extern void c__pdaf_get_local_type(
    int* localtype) noexcept nogil;

cdef extern void c__pdaf_get_memberid(
    int* memberid) noexcept nogil;

cdef extern void c__pdaf_get_obsmemberid(
    int* memberid) noexcept nogil;

cdef extern void c__pdaf_get_smootherens(CFI_cdesc_t* sens_point, 
    int* maxlag, int* status) noexcept nogil;

