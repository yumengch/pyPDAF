from .cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaflocal_set_indices(int* dim_l, 
    int* map) noexcept nogil;

cdef extern void c__pdaflocal_set_increment_weights(int* dim_l, 
    double* weights) noexcept nogil;

cdef extern void c__pdaflocal_clear_increment_weights() noexcept nogil;

