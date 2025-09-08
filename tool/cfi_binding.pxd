from libc.stdint cimport int8_t, int16_t
from libc.stddef cimport ptrdiff_t

cdef extern from "ISO_Fortran_binding.h":
    """
    typedef CFI_CDESC_T(1) CFI_cdesc_rank1;
    typedef CFI_CDESC_T(2) CFI_cdesc_rank2;
    typedef CFI_CDESC_T(3) CFI_cdesc_rank3;
    """
    ctypedef ptrdiff_t CFI_index_t;
    ctypedef int8_t CFI_rank_t;
    ctypedef int8_t CFI_attribute_t;
    ctypedef int16_t CFI_type_t;

    ctypedef struct CFI_dim_t:
        CFI_index_t extent;

    ctypedef struct CFI_cdesc_t:
        void *base_addr;
        CFI_rank_t rank;
        CFI_attribute_t attribute;
        CFI_dim_t* dim;

    ctypedef CFI_cdesc_t CFI_cdesc_rank1;
    ctypedef CFI_cdesc_t CFI_cdesc_rank2;
    ctypedef CFI_cdesc_t CFI_cdesc_rank3;

    cdef int CFI_attribute_pointer;
    cdef int CFI_attribute_allocatable;
    cdef int CFI_attribute_other;

    cdef CFI_type_t CFI_type_double;
    cdef CFI_type_t CFI_type_int;

    cdef extern void *CFI_address (const CFI_cdesc_t *, const CFI_index_t*) noexcept nogil;
    cdef extern int CFI_establish(CFI_cdesc_t *, void *, CFI_attribute_t,
                                  CFI_type_t, size_t, CFI_rank_t, const CFI_index_t*) noexcept nogil;

