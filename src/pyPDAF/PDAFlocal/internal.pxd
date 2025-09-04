
cdef extern void c__pdaflocal_g2l_cb(int* step, int* domain_p, int* dim_p,
    double* state_p, int* dim_l,
    double* state_l) noexcept nogil;

cdef extern void c__pdaflocal_l2g_cb(int* step, int* domain_p, int* dim_l,
    double* state_l, int* dim_p,
    double* state_p) noexcept nogil;