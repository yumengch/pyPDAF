from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdafomi_set_globalobs(
    int* globalobs_in) noexcept nogil;

cdef extern void c__pdafomi_diag_omit_by_inno() noexcept nogil;

cdef extern void c__pdafomi_cnt_dim_obs_l(int* i_obs,
    CFI_cdesc_t* coords_l) noexcept nogil;

cdef extern void c__pdafomi_cnt_dim_obs_l_noniso(int* i_obs,
    CFI_cdesc_t* coords_l) noexcept nogil;

cdef extern void c__pdafomi_init_obsarrays_l(int* i_obs,
    CFI_cdesc_t* coords_l, int* off_obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_init_obsarrays_l_noniso(int* i_obs,
    CFI_cdesc_t* coords_l, int* off_obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_g2l_obs(int* i_obs, CFI_cdesc_t* obs_f_all,
    CFI_cdesc_t* obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_init_obs_l(int* i_obs,
    CFI_cdesc_t* obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_init_obsvar_l(int* i_obs, double* meanvar_l,
    int* cnt_obs_l) noexcept nogil;

cdef extern void c__pdafomi_prodrinva_l(int* i_obs, int* nobs_all,
    int* ncols, CFI_cdesc_t* a_l, CFI_cdesc_t* c_l,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_prodrinva_hyb_l(int* i_obs, int* nobs_all,
    int* ncols, double* gamma, CFI_cdesc_t* a_l, CFI_cdesc_t* c_l,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_likelihood_l(int* i_obs,
    CFI_cdesc_t* resid_l_all, double* lhood_l,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_likelihood_hyb_l(int* i_obs,
    CFI_cdesc_t* resid_l_all, double* gamma, double* lhood_l,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_g2l_obs_internal(int* i_obs,
    CFI_cdesc_t* obs_f_one, int* offset_obs_l_all,
    CFI_cdesc_t* obs_l_all) noexcept nogil;

cdef extern void c__pdafomi_comp_dist2(int* i_obs, CFI_cdesc_t* coordsa,
    CFI_cdesc_t* coordsb, double* distance2,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_check_dist2(int* i_obs, CFI_cdesc_t* coordsa,
    CFI_cdesc_t* coordsb, double* distance2, bint* checkdist, int* verbose,
    int* cnt_obs) noexcept nogil;

cdef extern void c__pdafomi_check_dist2_noniso(int* i_obs,
    CFI_cdesc_t* coordsa, CFI_cdesc_t* coordsb, double* distance2,
    CFI_cdesc_t* dists, double* cradius, double* sradius, bint* checkdist,
    int* verbose, int* cnt_obs) noexcept nogil;

cdef extern void c__pdafomi_weights_l(int* verbose, int* nobs_l,
    int* ncols, int* locweight, CFI_cdesc_t* cradius, CFI_cdesc_t* sradius,
    CFI_cdesc_t* mata, CFI_cdesc_t* ivar_obs_l, CFI_cdesc_t* dist_l,
    CFI_cdesc_t* weight_l) noexcept nogil;

cdef extern void c__pdafomi_weights_l_sgnl(int* verbose, int* nobs_l,
    int* ncols, int* locweight, double* cradius, double* sradius,
    CFI_cdesc_t* mata, CFI_cdesc_t* ivar_obs_l, CFI_cdesc_t* dist_l,
    CFI_cdesc_t* weight_l) noexcept nogil;

cdef extern void c__pdafomi_omit_by_inno_l(int* i_obs, CFI_cdesc_t* inno_l,
    CFI_cdesc_t* obs_l_all, int* obsid, int* cnt_all,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_obsstats_l(
    int* screen) noexcept nogil;

cdef extern void c__pdafomi_dealloc() noexcept nogil;

cdef extern void c__pdafomi_ocoord_all(int* ncoord,
    CFI_cdesc_t* oc_all) noexcept nogil;

cdef extern void c__pdafomi_local_weight(int* wtype, int* rtype,
    double* cradius, double* sradius, double* distance, int* nrows,
    int* ncols, double* a, double* var_obs, double* weight,
    int* verbose) noexcept nogil;

cdef extern void c__pdafomi_check_dist2_loop(int* i_obs,
    CFI_cdesc_t* coordsa, int* cnt_obs,
    int* mode) noexcept nogil;

cdef extern void c__pdafomi_check_dist2_noniso_loop(int* i_obs,
    CFI_cdesc_t* coordsa, int* cnt_obs,
    int* mode) noexcept nogil;

cdef extern void c__pdafomi_obs_op_gatheronly(int* i_obs,
    CFI_cdesc_t* state_p,
    CFI_cdesc_t* obs_f_all) noexcept nogil;

cdef extern void c__pdafomi_obs_op_adj_gatheronly(int* i_obs,
    CFI_cdesc_t* obs_f_all,
    CFI_cdesc_t* state_p) noexcept nogil;

cdef extern void c__pdafomi_init_obs_f(int* i_obs, int* dim_obs_f,
    CFI_cdesc_t* obsstate_f, int* offset) noexcept nogil;

cdef extern void c__pdafomi_init_obsvars_f(int* i_obs, int* dim_obs_f,
    CFI_cdesc_t* var_f, int* offset) noexcept nogil;

cdef extern void c__pdafomi_init_obsvar_f(int* i_obs, double* meanvar,
    int* cnt_obs) noexcept nogil;

cdef extern void c__pdafomi_prodrinva(int* i_obs, int* ncols,
    CFI_cdesc_t* a_p, CFI_cdesc_t* c_p) noexcept nogil;

cdef extern void c__pdafomi_likelihood(int* i_obs, CFI_cdesc_t* resid,
    double* lhood) noexcept nogil;

cdef extern void c__pdafomi_add_obs_error(int* i_obs, int* nobs_all,
    CFI_cdesc_t* matc) noexcept nogil;

cdef extern void c__pdafomi_init_obscovar(int* i_obs, int* nobs_all,
    CFI_cdesc_t* covar, bint* isdiag) noexcept nogil;

cdef extern void c__pdafomi_init_obserr_f(int* i_obs,
    CFI_cdesc_t* obserr_f) noexcept nogil;

cdef extern void c__pdafomi_get_local_ids_obs_f(int* dim_obs_g,
    double* lradius, CFI_cdesc_t* oc_f, int* cnt_lim, CFI_cdesc_t* id_lim,
    int* disttype, CFI_cdesc_t* domainsize) noexcept nogil;

cdef extern void c__pdafomi_limit_obs_f(int* i_obs, int* offset,
    CFI_cdesc_t* obs_f_one,
    CFI_cdesc_t* obs_f_lim) noexcept nogil;

cdef extern void c__pdafomi_gather_dim_obs_f(int* dim_obs_p,
    int* dim_obs_f) noexcept nogil;

cdef extern void c__pdafomi_gather_obs_f_flex(int* dim_obs_p,
    CFI_cdesc_t* obs_p, CFI_cdesc_t* obs_f,
    int* status) noexcept nogil;

cdef extern void c__pdafomi_gather_obs_f2_flex(int* dim_obs_p,
    CFI_cdesc_t* coords_p, CFI_cdesc_t* coords_f, int* nrows,
    int* status) noexcept nogil;

cdef extern void c__pdafomi_omit_by_inno(int* i_obs, CFI_cdesc_t* inno_f,
    CFI_cdesc_t* obs_f_all, int* obsid,
    int* cnt_all) noexcept nogil;

cdef extern void c__pdafomi_obsstats(
    int* screen) noexcept nogil;

cdef extern void c__pdafomi_gather_obsdims() noexcept nogil;

