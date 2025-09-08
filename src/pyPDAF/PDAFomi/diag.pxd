from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_diag_dimobs(
    CFI_cdesc_t* dim_obs_ptr) noexcept nogil;

cdef extern void c__pdafomi_diag_get_hx(int* id_obs, int* dim_obs_diag,
    CFI_cdesc_t* hx_p_ptr) noexcept nogil;

cdef extern void c__pdafomi_diag_get_hxmean(int* id_obs, int* dim_obs_diag,
    CFI_cdesc_t* hxmean_p_ptr) noexcept nogil;

cdef extern void c__pdafomi_diag_get_ivar(int* id_obs, int* dim_obs_diag,
    CFI_cdesc_t* ivar_ptr) noexcept nogil;

cdef extern void c__pdafomi_diag_get_obs(int* id_obs, int* dim_obs_diag,
    int* ncoord, CFI_cdesc_t* obs_p_ptr,
    CFI_cdesc_t* ocoord_p_ptr) noexcept nogil;

cdef extern void c__pdafomi_diag_nobstypes(
    int* nobs) noexcept nogil;

cdef extern void c__pdafomi_diag_obs_rmsd(int* nobs,
    CFI_cdesc_t* rmsd_pointer, int* verbose) noexcept nogil;

cdef extern void c__pdafomi_diag_stats(int* nobs,
    CFI_cdesc_t* obsstats_ptr, int* verbose) noexcept nogil;

