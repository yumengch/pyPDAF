from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_init_obs_f_cb(int* step, int* dim_obs_f,
    double* observation_f) noexcept nogil;

cdef extern void c__pdafomi_init_obsvar_cb(int* step, int* dim_obs_p,
    double* obs_p, double* meanvar) noexcept nogil;

cdef extern void c__pdafomi_init_obsvars_f_cb(int* step, int* dim_obs_f,
    double* var_f) noexcept nogil;

cdef extern void c__pdafomi_g2l_obs_cb(int* domain_p, int* step,
    int* dim_obs_f, int* dim_obs_l, double* ostate_f,
    double* ostate_l) noexcept nogil;

cdef extern void c__pdafomi_init_obs_l_cb(int* domain_p, int* step,
    int* dim_obs_l, double* observation_l) noexcept nogil;

cdef extern void c__pdafomi_init_obsvar_l_cb(int* domain_p, int* step,
    int* dim_obs_l, double* obs_l,
    double* meanvar_l) noexcept nogil;

cdef extern void c__pdafomi_prodrinva_l_cb(int* domain_p, int* step,
    int* dim_obs_l, int* rank, double* obs_l, double* a_l,
    double* c_l) noexcept nogil;

cdef extern void c__pdafomi_likelihood_l_cb(int* domain_p, int* step,
    int* dim_obs_l, double* obs_l, double* resid_l,
    double* lhood_l) noexcept nogil;

cdef extern void c__pdafomi_prodrinva_hyb_l_cb(int* domain_p, int* step,
    int* dim_obs_l, int* rank, double* obs_l, double* alpha, double* a_l,
    double* c_l) noexcept nogil;

cdef extern void c__pdafomi_likelihood_hyb_l_cb(int* domain_p, int* step,
    int* dim_obs_l, double* obs_l, double* resid_l, double* alpha,
    double* lhood_l) noexcept nogil;

cdef extern void c__pdafomi_prodrinva_cb(int* step, int* dim_obs_p,
    int* ncol, double* obs_p, double* a_p,
    double* c_p) noexcept nogil;

cdef extern void c__pdafomi_likelihood_cb(int* step, int* dim_obs,
    double* obs, double* resid, double* lhood) noexcept nogil;

cdef extern void c__pdafomi_add_obs_error_cb(int* step, int* dim_obs_p,
    double* c_p) noexcept nogil;

cdef extern void c__pdafomi_init_obscovar_cb(int* step, int* dim_obs,
    int* dim_obs_p, double* covar, double* m_state_p,
    bint* isdiag) noexcept nogil;

cdef extern void c__pdafomi_init_obserr_f_cb(int* step, int* dim_obs_f,
    double* obs_f, double* obserr_f) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_cb(int* dim_p, int* dim_obs,
    double* hp_p, double* hph) noexcept nogil;

cdef extern void c__pdafomi_localize_covar_serial_cb(int* iobs, int* dim_p,
    int* dim_obs, double* hp_p, double* hxy_p) noexcept nogil;

cdef extern void c__pdafomi_omit_by_inno_l_cb(int* domain_p,
    int* dim_obs_l, double* resid_l,
    double* obs_l) noexcept nogil;

cdef extern void c__pdafomi_omit_by_inno_cb(int* dim_obs_f,
    double* resid_f, double* obs_f) noexcept nogil;
