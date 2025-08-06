from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_init(int* n_obs) noexcept nogil;

cdef extern void c__pdafomi_init_local() noexcept nogil;

cdef extern void c__pdafomi_check_error(
    int* flag) noexcept nogil;

cdef extern void c__pdafomi_gather_obs(int* i_obs, int* dim_obs_p,
    CFI_cdesc_t* obs_p, CFI_cdesc_t* ivar_obs_p, CFI_cdesc_t* ocoord_p,
    int* ncoord, double* lradius,
    int* dim_obs_f) noexcept nogil;

cdef extern void c__pdafomi_gather_obsstate(int* i_obs,
    CFI_cdesc_t* obsstate_p,
    CFI_cdesc_t* obsstate_f) noexcept nogil;

cdef extern void c__pdafomi_get_interp_coeff_tri(CFI_cdesc_t* gpc,
    CFI_cdesc_t* oc, CFI_cdesc_t* icoeff) noexcept nogil;

cdef extern void c__pdafomi_get_interp_coeff_lin1d(CFI_cdesc_t* gpc,
    double* oc, CFI_cdesc_t* icoeff) noexcept nogil;

cdef extern void c__pdafomi_get_interp_coeff_lin(int* num_gp, int* n_dim,
    CFI_cdesc_t* gpc, CFI_cdesc_t* oc,
    CFI_cdesc_t* icoeff) noexcept nogil;

cdef extern void c__pdafomi_init_dim_obs_l_iso(int* i_obs,
    CFI_cdesc_t* coords_l, int* locweight, double* cradius,
    double* sradius, int* cnt_obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_init_dim_obs_l_noniso(int* i_obs,
    CFI_cdesc_t* coords_l, int* locweight, CFI_cdesc_t* cradius,
    CFI_cdesc_t* sradius, int* cnt_obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_init_dim_obs_l_noniso_locweights(int* i_obs,
    CFI_cdesc_t* coords_l, CFI_cdesc_t* locweights, CFI_cdesc_t* cradius,
    CFI_cdesc_t* sradius, int* cnt_obs_l) noexcept nogil;

cdef extern void c__pdafomi_obs_op_gridpoint(int* i_obs,
    CFI_cdesc_t* state_p,
    CFI_cdesc_t* obs_f_all) noexcept nogil;

cdef extern void c__pdafomi_obs_op_gridavg(int* i_obs, int* nrows,
    CFI_cdesc_t* state_p,
    CFI_cdesc_t* obs_f_all) noexcept nogil;

cdef extern void c__pdafomi_obs_op_interp_lin(int* i_obs, int* nrows,
    CFI_cdesc_t* state_p,
    CFI_cdesc_t* obs_f_all) noexcept nogil;

cdef extern void c__pdafomi_obs_op_adj_gridpoint(int* i_obs,
    CFI_cdesc_t* obs_f_all,
    CFI_cdesc_t* state_p) noexcept nogil;

cdef extern void c__pdafomi_obs_op_adj_gridavg(int* i_obs, int* nrows,
    CFI_cdesc_t* obs_f_all,
    CFI_cdesc_t* state_p) noexcept nogil;

cdef extern void c__pdafomi_obs_op_adj_interp_lin(int* i_obs, int* nrows,
    CFI_cdesc_t* obs_f_all,
    CFI_cdesc_t* state_p) noexcept nogil;

cdef extern void c__pdafomi_observation_localization_weights(int* i_obs,
    int* ncols, CFI_cdesc_t* a_l, double* weight,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_set_debug_flag(
    int* debugval) noexcept nogil;

cdef extern void c__pdafomi_set_dim_obs_l(int* i_obs, int* cnt_obs_l_all,
    int* cnt_obs_l) noexcept nogil;

cdef extern void c__pdafomi_set_localization(int* i_obs, double* cradius,
    double* sradius, int* locweight) noexcept nogil;

cdef extern void c__pdafomi_set_localization_noniso(int* i_obs,
    int* nradii, double* cradius, double* sradius, int* locweight,
    int* locweight_v) noexcept nogil;

cdef extern void c__pdafomi_set_localize_covar_iso(int* i_obs, int* dim,
    int* ncoords, CFI_cdesc_t* coords, int* locweight, double* cradius,
    double* sradius) noexcept nogil;

cdef extern void c__pdafomi_set_localize_covar_noniso(int* i_obs, int* dim,
    int* ncoords, CFI_cdesc_t* coords, int* locweight,
    CFI_cdesc_t* cradius, CFI_cdesc_t* sradius) noexcept nogil;

cdef extern void c__pdafomi_set_localize_covar_noniso_locweights(
    int* i_obs, int* dim, int* ncoords, CFI_cdesc_t* coords,
    CFI_cdesc_t* locweights, CFI_cdesc_t* cradius,
    CFI_cdesc_t* sradius) noexcept nogil;

cdef extern void c__pdafomi_set_obs_diag(
    int* diag) noexcept nogil;

cdef extern void c__pdafomi_set_domain_limits(
    double* lim_coords) noexcept nogil;

cdef extern void c__pdafomi_get_domain_limits_unstr(int* npoints_p,
    CFI_cdesc_t* coords_p) noexcept nogil;

cdef extern void c__pdafomi_store_obs_l_index(int* i_obs, int* idx,
    int* id_obs_l, double* distance, double* cradius_l,
    double* sradius_l) noexcept nogil;

cdef extern void c__pdafomi_store_obs_l_index_vdist(int* i_obs, int* idx,
    int* id_obs_l, double* distance, double* cradius_l, double* sradius_l,
    double* vdist) noexcept nogil;

