from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_localize_covar_iso(int* i_obs, int* dim,
    int* locweight, double* cradius, double* sradius, CFI_cdesc_t* coords,
    CFI_cdesc_t* hp, CFI_cdesc_t* hph) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_noniso_locweights(int* i_obs,
    int* dim, CFI_cdesc_t* locweights, CFI_cdesc_t* cradius,
    CFI_cdesc_t* sradius, CFI_cdesc_t* coords, CFI_cdesc_t* hp,
    CFI_cdesc_t* hph) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_noniso(int* i_obs, int* dim,
    int* locweight, CFI_cdesc_t* cradius, CFI_cdesc_t* sradius,
    CFI_cdesc_t* coords, CFI_cdesc_t* hp,
    CFI_cdesc_t* hph) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_serial_iso(int* i_obs,
    int* iobs_all, int* dim, int* dim_obs, int* locweight, double* cradius,
    double* sradius, CFI_cdesc_t* coords, CFI_cdesc_t* hp,
    CFI_cdesc_t* hxy) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_serial_noniso_locweights(
    int* i_obs, int* iobs_all, int* dim, int* dim_obs,
    CFI_cdesc_t* locweights, CFI_cdesc_t* cradius, CFI_cdesc_t* sradius,
    CFI_cdesc_t* coords, CFI_cdesc_t* hp,
    CFI_cdesc_t* hxy) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_serial_noniso(int* i_obs,
    int* iobs_all, int* dim, int* dim_obs, int* locweight,
    CFI_cdesc_t* cradius, CFI_cdesc_t* sradius, CFI_cdesc_t* coords,
    CFI_cdesc_t* hp, CFI_cdesc_t* hxy) noexcept nogil;

cdef extern void c__pdafomi_deallocate_obs(
    int* i_obs) noexcept nogil;

