from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_set_doassim(int* i_obs,
    int* doassim) noexcept nogil;

cdef extern void c__pdafomi_set_disttype(int* i_obs,
    int* disttype) noexcept nogil;

cdef extern void c__pdafomi_set_ncoord(int* i_obs,
    int* ncoord) noexcept nogil;

cdef extern void c__pdafomi_set_obs_err_type(int* i_obs,
    int* obs_err_type) noexcept nogil;

cdef extern void c__pdafomi_set_use_global_obs(int* i_obs,
    int* use_global_obs) noexcept nogil;

cdef extern void c__pdafomi_set_inno_omit(int* i_obs,
    double* inno_omit) noexcept nogil;

cdef extern void c__pdafomi_set_inno_omit_ivar(int* i_obs,
    double* inno_omit_ivar) noexcept nogil;

cdef extern void c__pdafomi_set_id_obs_p(int* i_obs, int* nrows,
    int* dim_obs_p, int* id_obs_p) noexcept nogil;

cdef extern void c__pdafomi_set_icoeff_p(int* i_obs, int* nrows,
    int* dim_obs_p, double* icoeff_p) noexcept nogil;

cdef extern void c__pdafomi_set_domainsize(int* i_obs, int* ncoord,
    double* domainsize) noexcept nogil;

cdef extern void c__pdafomi_set_name(int* i_obs, char* obsname) noexcept nogil;