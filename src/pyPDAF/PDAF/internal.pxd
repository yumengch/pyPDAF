from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_set_forget_local(int* domain, int* step,
    int* dim_obs_l, int* dim_ens, double* hx_l, double* hxbar_l,
    double* obs_l,
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    double* forget, double* aforget) noexcept nogil;

cdef extern void c__pdaf_fcst_operations(int* step,
    void (*c__collect_state_pdaf)(int* , double* ),
    void (*c__distribute_state_pdaf)(int* , double* ),
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_letkf_ana_t(int* domain_p, int* step, int* dim_l,
    int* dim_obs_l, int* dim_ens, double* state_l, double* ainv_l,
    double* ens_l, double* hz_l, double* hxbar_l, double* obs_l,
    double* rndmat, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* type_trans, int* screen, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdafseik_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, int* rank, double* state_p, double* uinv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf3dvar_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, int* dim_cvec, double* state_p, double* ainv,
    double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdafen3dvar_update_estkf(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec_ens, double* state_p,
    double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdafen3dvar_update_lestkf(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec_ens, double* state_p,
    double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_dim_obs_f_pdaf)(int* , int* ),
    void (*c__obs_op_f_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_f_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdafetkf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* dim_lag, double* sens_p,
    int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_netf_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p, double* rndmat,
    double* t, int* type_forget, double* forget, int* type_winf,
    double* limit_winf, int* type_noise, double* noise_amp, double* hz_p,
    double* obs_p,
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_netf_smoothert(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, double* ens_p, double* rndmat, double* t,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ),
    int* screen, int* flag) noexcept nogil;

cdef extern void c__pdaf_smoother_netf(int* dim_p, int* dim_ens,
    int* dim_lag, double* ainv, double* sens_p, int* cnt_maxlag,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_lnetf_ana(int* domain_p, int* step, int* dim_l,
    int* dim_obs_l, int* dim_ens, double* ens_l, double* hx_l,
    double* obs_l, double* rndmat,
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    int* type_forget, double* forget, int* type_winf, double* limit_winf,
    int* cnt_small_svals, double* eff_dimens, double* t, int* screen,
    int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_lnetf_smoothert(int* domain_p, int* step,
    int* dim_obs_f, int* dim_obs_l, int* dim_ens, double* hx_f,
    double* rndmat,
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    int* screen, double* t, int* flag) noexcept nogil;

cdef extern void c__pdaf_smoother_lnetf(int* domain_p, int* step,
    int* dim_p, int* dim_l, int* dim_ens, int* dim_lag, double* ainv,
    double* ens_l, double* sens_p, int* cnt_maxlag,
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    int* screen) noexcept nogil;

cdef extern void c__pdaf_memcount_ini(
    int* ncounters) noexcept nogil;

cdef extern void c__pdaf_memcount_define(char* stortype,
    int* wordlength) noexcept nogil;

cdef extern void c__pdaf_memcount(int* id, char* stortype,
    int* dim) noexcept nogil;

cdef extern void c__pdaf_init_filters(int* type_filter, int* subtype,
    int* param_int, int* dim_pint, double* param_real, int* dim_preal,
    char* filterstr, bint* ensemblefilter, bint* fixedbasis, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_alloc_filters(char* filterstr, int* subtype,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_configinfo_filters(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_options_filters(
    int* type_filter) noexcept nogil;

cdef extern void c__pdaf_print_info_filters(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_allreduce(int* val_p, int* val_g, int* mpitype,
    int* mpiop, int* status) noexcept nogil;

cdef extern void c__pdaflseik_update(int* step, int* dim_p, int* dim_obs_f,
    int* dim_ens, int* rank, double* state_p, double* uinv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_ensrf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_ensrf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_ensrf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_ensrf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_ensrf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_ensrf_options() noexcept nogil;

cdef extern void c__pdaf_ensrf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_estkf_ana_fixed(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank, double* state_p, double* ainv,
    double* ens_p, double* hl_p, double* hxbar_p, double* obs_p,
    double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* type_sqrt, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_etkf_ana_fixed(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, double* state_p, double* ainv,
    double* ens_p, double* hz_p, double* hxbar_p, double* obs_p,
    double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdafestkf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* envar_mode, int* dim_lag,
    double* sens_p, int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaflknetf_update_step(int* step, int* dim_p,
    int* dim_obs_f, int* dim_ens, double* state_p, double* ainv,
    double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_hyb_l_pdaf)(int* , int* , int* , int* , double* ,
                                    double* , double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    void (*c__likelihood_hyb_l_pdaf)(int* , int* , int* , double* ,
                                     double* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdafletkf_update(int* step, int* dim_p, int* dim_obs_f,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* dim_lag, double* sens_p,
    int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_lseik_ana_trans(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, int* rank, double* state_l,
    double* uinv_l, double* ens_l, double* hl_l, double* hxbar_l,
    double* obs_l, double* omegat_in, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* nm1vsn, int* type_sqrt, int* screen, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_en3dvar_optim_lbfgs(int* step, int* dim_p,
    int* dim_ens, int* dim_cvec_p, int* dim_obs_p, double* ens_p,
    double* obs_p, double* dy_p, double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_en3dvar_optim_cgplus(int* step, int* dim_p,
    int* dim_ens, int* dim_cvec_p, int* dim_obs_p, double* ens_p,
    double* obs_p, double* dy_p, double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_en3dvar_optim_cg(int* step, int* dim_p,
    int* dim_ens, int* dim_cvec_p, int* dim_obs_p, double* ens_p,
    double* obs_p, double* dy_p, double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_en3dvar_costf_cvt(int* step, int* iter,
    int* dim_p, int* dim_ens, int* dim_cvec_p, int* dim_obs_p,
    double* ens_p, double* obs_p, double* dy_p, double* v_p, double* j_tot,
    double* gradj,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel) noexcept nogil;

cdef extern void c__pdaf_en3dvar_costf_cg_cvt(int* step, int* iter,
    int* dim_p, int* dim_ens, int* dim_cvec_p, int* dim_obs_p,
    double* ens_p, double* obs_p, double* dy_p, double* v_p, double* d_p,
    double* j_tot, double* gradj, double* hessjd,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel) noexcept nogil;

cdef extern void c__pdaf_gather_ens(int* dim_p, int* dim_ens_p,
    CFI_cdesc_t* ens, int* screen) noexcept nogil;

cdef extern void c__pdaf_scatter_ens(int* dim_p, int* dim_ens_p,
    CFI_cdesc_t* ens, CFI_cdesc_t* state,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_hyb3dvar_optim_lbfgs(int* step, int* dim_p,
    int* dim_ens, int* dim_cv_par_p, int* dim_cv_ens_p, int* dim_obs_p,
    double* ens_p, double* obs_p, double* dy_p, double* v_par_p,
    double* v_ens_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, double* beta_3dvar,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_hyb3dvar_optim_cgplus(int* step, int* dim_p,
    int* dim_ens, int* dim_cv_par_p, int* dim_cv_ens_p, int* dim_obs_p,
    double* ens_p, double* obs_p, double* dy_p, double* v_par_p,
    double* v_ens_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, double* beta_3dvar,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_hyb3dvar_optim_cg(int* step, int* dim_p,
    int* dim_ens, int* dim_cv_par_p, int* dim_cv_ens_p, int* dim_obs_p,
    double* ens_p, double* obs_p, double* dy_p, double* v_par_p,
    double* v_ens_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, double* beta_3dvar,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_hyb3dvar_costf_cvt(int* step, int* iter,
    int* dim_p, int* dim_ens, int* dim_cv_p, int* dim_cv_par_p,
    int* dim_cv_ens_p, int* dim_obs_p, double* ens_p, double* obs_p,
    double* dy_p, double* v_par_p, double* v_ens_p, double* v_p,
    double* j_tot, double* gradj,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, double* beta) noexcept nogil;

cdef extern void c__pdaf_hyb3dvar_costf_cg_cvt(int* step, int* iter,
    int* dim_p, int* dim_ens, int* dim_cv_par_p, int* dim_cv_ens_p,
    int* dim_obs_p, double* ens_p, double* obs_p, double* dy_p,
    double* v_par_p, double* v_ens_p, double* d_par_p, double* d_ens_p,
    double* j_tot, double* gradj_par, double* gradj_ens,
    double* hessjd_par, double* hessjd_ens,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, double* beta) noexcept nogil;

cdef extern void c__pdaf_print_version() noexcept nogil;

cdef extern void c__pdafen3dvar_analysis_cvt(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec_ens, double* state_p,
    double* ens_p, double* state_inc_p, double* hxbar_p, double* obs_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* screen, int* type_opt, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_sisort(int* n,
    double* veca) noexcept nogil;

cdef extern void c__pdaf_enkf_ana_rlm(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank_ana, double* state_p,
    double* ens_p, double* hzb, double* hx_p, double* hxbar_p,
    double* obs_p,
    void (*c__add_obs_err_pdaf)(int* , int* , double* ),
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_smoother_enkf(int* dim_p, int* dim_ens,
    int* dim_lag, double* ainv, double* sens_p, int* cnt_maxlag,
    double* forget, int* screen) noexcept nogil;

cdef extern void c__pdafensrf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__localize_covar_serial_pdaf)(int* , int* , int* , double* ,
                                          double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_pf_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p, int* type_resample,
    int* type_winf, double* limit_winf, int* type_noise, double* noise_amp,
    double* hz_p, double* obs_p,
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_pf_resampling(int* method, int* nin, int* nout,
    double* weights, int* ids, int* screen) noexcept nogil;

cdef extern void c__pdaf_mvnormalize(int* mode, int* dim_state,
    int* dim_field, int* offset, int* ncol, double* states, double* stddev,
    int* status) noexcept nogil;

cdef extern void c__pdaf_3dvar_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_3dvar_alloc(int* subtype,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_3dvar_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_3dvar_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_3dvar_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_3dvar_options() noexcept nogil;

cdef extern void c__pdaf_3dvar_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_reset_dim_ens(int* dim_ens_in,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_reset_dim_p(int* dim_p_in,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_3dvar_optim_lbfgs(int* step, int* dim_p,
    int* dim_cvec_p, int* dim_obs_p, double* obs_p, double* dy_p,
    double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_3dvar_optim_cgplus(int* step, int* dim_p,
    int* dim_cvec_p, int* dim_obs_p, double* obs_p, double* dy_p,
    double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_3dvar_optim_cg(int* step, int* dim_p,
    int* dim_cvec_p, int* dim_obs_p, double* obs_p, double* dy_p,
    double* v_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel, int* screen) noexcept nogil;

cdef extern void c__pdaf_3dvar_costf_cvt(int* step, int* iter, int* dim_p,
    int* dim_cvec_p, int* dim_obs_p, double* obs_p, double* dy_p,
    double* v_p, double* j_tot, double* gradj,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel) noexcept nogil;

cdef extern void c__pdaf_3dvar_costf_cg_cvt(int* step, int* iter,
    int* dim_p, int* dim_cvec_p, int* dim_obs_p, double* obs_p,
    double* dy_p, double* v_p, double* d_p, double* j_tot, double* gradj,
    double* hessjd,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* opt_parallel) noexcept nogil;

cdef extern void c__pdaf_lknetf_analysis_t(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, double* state_l,
    double* ainv_l, double* ens_l, double* hx_l, double* hxbar_l,
    double* obs_l, double* rndmat, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    int* screen, int* type_forget, double* eff_dimens, int* type_hyb,
    double* hyb_g, double* hyb_k, double* gamma, double* skew_mabs,
    double* kurt_mabs, int* flag) noexcept nogil;

cdef extern void c__pdaf_get_ensstats(CFI_cdesc_t* skew_ptr,
    CFI_cdesc_t* kurt_ptr, int* status) noexcept nogil;

cdef extern void c__pdaf_estkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_estkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_estkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_estkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_estkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_estkf_options() noexcept nogil;

cdef extern void c__pdaf_estkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_gen_obs(int* step, int* dim_p, int* dim_obs_f,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_f_pdaf)(int* , int* ),
    void (*c__obs_op_f_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__get_obs_f_pdaf)(int* , int* , double* ),
    void (*c__init_obserr_f_pdaf)(int* , int* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* flag) noexcept nogil;

cdef extern void c__pdafobs_init(int* step, int* dim_p, int* dim_ens,
    int* dim_obs_p, double* state_p, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    int* screen, int* debug, bint* do_ens_mean, bint* do_init_dim,
    bint* do_hx, bint* do_hxbar,
    bint* do_init_obs) noexcept nogil;

cdef extern void c__pdafobs_init_local(int* domain_p, int* step,
    int* dim_obs_l, int* dim_obs_f, int* dim_ens,
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    int* debug) noexcept nogil;

cdef extern void c__pdafobs_init_obsvars(int* step, int* dim_obs_p,
    void (*c__init_obsvars_pdaf)(int* , int* ,
                                 double* )) noexcept nogil;


cdef extern void c__pdafobs_dealloc() noexcept nogil;

cdef extern void c__pdafobs_dealloc_local() noexcept nogil;

cdef extern void c__pdaf_netf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_netf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_netf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_netf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_netf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_netf_options() noexcept nogil;

cdef extern void c__pdaf_netf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_lenkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lenkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lenkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_lenkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lenkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lenkf_options() noexcept nogil;

cdef extern void c__pdaf_lenkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_lseik_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lseik_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lseik_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_lseik_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lseik_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lseik_options() noexcept nogil;

cdef extern void c__pdaf_lseik_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_etkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_etkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_etkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_etkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_etkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_etkf_options() noexcept nogil;

cdef extern void c__pdaf_etkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaflenkf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__add_obs_err_pdaf)(int* , int* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_pf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_pf_alloc(int* outflag) noexcept nogil;

cdef extern void c__pdaf_pf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_pf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_pf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_pf_options() noexcept nogil;

cdef extern void c__pdaf_pf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_lknetf_ana_letkft(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, double* state_l,
    double* ainv_l, double* ens_l, double* hz_l, double* hxbar_l,
    double* obs_l, double* rndmat, double* forget,
    void (*c__prodrinva_hyb_l_pdaf)(int* , int* , int* , int* , double* ,
                                    double* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    double* gamma, int* screen, int* type_forget,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_ana_lnetf(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, double* ens_l, double* hx_l,
    double* rndmat, double* obs_l,
    void (*c__likelihood_hyb_l_pdaf)(int* , int* , int* , double* ,
                                     double* , double* , double* ),
    int* cnt_small_svals, double* n_eff_all, double* gamma, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_enkf_ana_rsm(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank_ana, double* state_p,
    double* ens_p, double* hx_p, double* hxbar_p, double* obs_p,
    void (*c__add_obs_err_pdaf)(int* , int* , double* ),
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lknetf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lknetf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_lknetf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_options() noexcept nogil;

cdef extern void c__pdaf_lknetf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_lknetf_alpha_neff(int* dim_ens, double* weights,
    double* hlimit, double* alpha) noexcept nogil;

cdef extern void c__pdaf_lknetf_compute_gamma(int* domain_p, int* step,
    int* dim_obs_l, int* dim_ens, double* hx_l, double* hxbar_l,
    double* obs_l, int* type_hyb, double* hyb_g, double* hyb_k,
    double* gamma, double* n_eff_out, double* skew_mabs, double* kurt_mabs,
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    int* screen, int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_set_gamma(int* domain_p, int* dim_obs_l,
    int* dim_ens, double* hx_l, double* hxbar_l, double* weights,
    int* type_hyb, double* hyb_g, double* hyb_k, double* gamma,
    double* n_eff_out, double* maskew, double* makurt, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lknetf_reset_gamma(
    double* gamma_in) noexcept nogil;

cdef extern void c__pdafhyb3dvar_analysis_cvt(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec, int* dim_cvec_ens,
    double* beta_3dvar, double* state_p, double* ens_p,
    double* state_inc_p, double* hxbar_p, double* obs_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* screen, int* type_opt, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf3dvar_analysis_cvt(int* step, int* dim_p,
    int* dim_obs_p, int* dim_cvec, double* state_p, double* hxbar_p,
    double* obs_p,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    int* screen, int* type_opt, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lestkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lestkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lestkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_lestkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lestkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lestkf_options() noexcept nogil;

cdef extern void c__pdaf_lestkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_seik_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, int* rank, double* state_p, double* uinv, double* ens_p,
    double* hl_p, double* hxbar_p, double* obs_p, double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_resample(int* subtype, int* dim_p,
    int* dim_ens, int* rank, double* uinv, double* state_p, double* enst_p,
    int* type_sqrt, int* type_trans, int* nm1vsn, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lseik_ana(int* domain_p, int* step, int* dim_l,
    int* dim_obs_l, int* dim_ens, int* rank, double* state_l,
    double* uinv_l, double* ens_l, double* hl_l, double* hxbar_l,
    double* obs_l, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_lseik_resample(int* domain_p, int* subtype,
    int* dim_l, int* dim_ens, int* rank, double* uinv_l, double* state_l,
    double* ens_l, double* omegat_in, int* type_sqrt, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_prepost(
    void (*c__collect_state_pdaf)(int* , double* ),
    void (*c__distribute_state_pdaf)(int* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ),
    int* outflag) noexcept nogil;

cdef extern void c__pdafenkf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__add_obs_err_pdaf)(int* , int* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* dim_lag, double* sens_p,
    int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_init_parallel(int* dim_ens, bint* ensemblefilter,
    bint* fixedbasis, int* comm_model, int* in_comm_filter,
    int* in_comm_couple, int* in_n_modeltasks, int* in_task_id,
    int* screen, int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_seik_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_seik_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_seik_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_options() noexcept nogil;

cdef extern void c__pdaf_seik_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdafnetf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* dim_lag, double* sens_p,
    int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_ana_newt(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank, double* state_p, double* uinv,
    double* ens_p, double* hl_p, double* hxbar_p, double* obs_p,
    double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_resample_newt(int* subtype, int* dim_p,
    int* dim_ens, int* rank, double* uinv, double* state_p, double* ens_p,
    int* type_sqrt, int* type_trans, int* nm1vsn, int* screen,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lenkf_ana_rsm(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank_ana, double* state_p,
    double* ens_p, double* hx_p, double* hxbar_p, double* obs_p,
    void (*c__add_obs_err_pdaf)(int* , int* , double* ),
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_lestkf_ana(int* domain_p, int* step, int* dim_l,
    int* dim_obs_l, int* dim_ens, int* rank, double* state_l,
    double* ainv_l, double* ens_l, double* hl_l, double* hxbar_l,
    double* obs_l, double* omegat_in, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* envar_mode, int* type_sqrt, double* ta, int* screen, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaflestkf_update(int* step, int* dim_p,
    int* dim_obs_f, int* dim_ens, int* rank, double* state_p, double* ainv,
    double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* envar_mode, int* dim_lag,
    double* sens_p, int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_lnetf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lnetf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_lnetf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_lnetf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lnetf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_lnetf_options() noexcept nogil;

cdef extern void c__pdaf_lnetf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_enkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_enkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_enkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_enkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_enkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_enkf_options() noexcept nogil;

cdef extern void c__pdaf_enkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_enkf_gather_resid(int* dim_obs, int* dim_obs_p,
    int* dim_ens, double* resid_p,
    double* resid) noexcept nogil;

cdef extern void c__pdaf_enkf_obs_ensemble(int* step, int* dim_obs_p,
    int* dim_obs, int* dim_ens, double* obsens_p, double* obs_p,
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* ,
                                   bint* ),
    int* screen, int* flag) noexcept nogil;

cdef extern void c__pdafpf_update(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_generate_rndmat(int* dim, double* rndmat,
    int* mattype) noexcept nogil;

cdef extern void c__pdaf_print_domain_stats(
    int* n_domains_p) noexcept nogil;

cdef extern void c__pdaf_init_local_obsstats() noexcept nogil;

cdef extern void c__pdaf_incr_local_obsstats(
    int* dim_obs_l) noexcept nogil;

cdef extern void c__pdaf_print_local_obsstats(int* screen,
    int* n_domains_with_obs) noexcept nogil;

cdef extern void c__pdaf_seik_matrixt(int* dim, int* dim_ens,
    double* a) noexcept nogil;

cdef extern void c__pdaf_seik_ttimesa(int* rank, int* dim_col, double* a,
    double* b) noexcept nogil;

cdef extern void c__pdaf_seik_omega(int* rank, double* omega,
    int* omegatype, int* screen) noexcept nogil;

cdef extern void c__pdaf_seik_uinv(int* rank,
    double* uinv) noexcept nogil;

cdef extern void c__pdaf_ens_omega(int* seed, int* r, int* dim_ens,
    double* omega, double* norm, int* otype,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_estkf_omegaa(int* rank, int* dim_col, double* a,
    double* b) noexcept nogil;

cdef extern void c__pdaf_estkf_aomega(int* dim, int* dim_ens,
    double* a) noexcept nogil;

cdef extern void c__pdaf_subtract_rowmean(int* dim, int* dim_ens,
    double* a) noexcept nogil;

cdef extern void c__pdaf_subtract_colmean(int* dim_ens, int* dim,
    double* a) noexcept nogil;

cdef extern void c__pdaf_add_particle_noise(int* dim_p, int* dim_ens,
    double* state_p, double* ens_p, int* type_noise, double* noise_amp,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_inflate_weights(int* screen, int* dim_ens,
    double* alpha, double* weights) noexcept nogil;

cdef extern void c__pdaf_inflate_ens(int* dim, int* dim_ens,
    double* meanstate, double* ens, double* forget,
    bint* do_ensmean) noexcept nogil;

cdef extern void c__pdaf_alloc(int* dim_p, int* dim_ens, int* dim_ens_task,
    int* dim_es, int* dim_bias_p, int* dim_lag, int* statetask,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_smoothing(int* dim_p, int* dim_ens, int* dim_lag,
    double* ainv, double* sens_p, int* cnt_maxlag, double* forget,
    int* screen) noexcept nogil;

cdef extern void c__pdaf_smoothing_local(int* domain_p, int* step,
    int* dim_p, int* dim_l, int* dim_ens, int* dim_lag, double* ainv,
    double* ens_l, double* sens_p, int* cnt_maxlag,
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    double* forget, int* screen) noexcept nogil;

cdef extern void c__pdaf_smoother_shift(int* dim_p, int* dim_ens,
    int* dim_lag, double* ens_p, double* sens_p, int* cnt_maxlag,
    int* screen) noexcept nogil;

cdef extern void c__pdaflknetf_update_sync(int* step, int* dim_p,
    int* dim_obs_f, int* dim_ens, double* state_p, double* ainv,
    double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_etkf_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    double* hz_p, double* hxbar_p, double* obs_p, double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* type_trans, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_letkf_ana(int* domain_p, int* step, int* dim_l,
    int* dim_obs_l, int* dim_ens, double* state_l, double* ainv_l,
    double* ens_l, double* hz_l, double* hxbar_l, double* obs_l,
    double* rndmat, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* type_trans, int* screen, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_letkf_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_letkf_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_letkf_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_letkf_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_letkf_set_rparam(int* id, double* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_letkf_options() noexcept nogil;

cdef extern void c__pdaf_letkf_memtime(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_estkf_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, int* rank, double* state_p, double* ainv, double* ens_p,
    double* hl_p, double* hxbar_p, double* obs_p, double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* envar_mode, int* type_sqrt, int* type_trans,
    double* ta, int* debug, int* flag) noexcept nogil;

cdef extern void c__pdaf_ensrf_ana(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ens_p, double* hx_p,
    double* hxbar_p, double* obs_p, double* var_obs_p,
    void (*c__localize_covar_serial_pdaf)(int* , int* , int* , double* ,
                                          double* ),
    int* screen, int* debug) noexcept nogil;

cdef extern void c__pdaf_ensrf_ana_2step(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, double* state_p, double* ens_p,
    double* hx_p, double* hxbar_p, double* obs_p, double* var_obs_p,
    void (*c__localize_covar_serial_pdaf)(int* , int* , int* , double* ,
                                          double* ),
    int* screen, int* debug) noexcept nogil;

cdef extern void c__pdaflnetf_update(int* step, int* dim_p, int* dim_obs_f,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* ,
                                 double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* screen, int* subtype, int* dim_lag, double* sens_p,
    int* cnt_maxlag, int* flag) noexcept nogil;

cdef extern void c__pdaf_seik_ana_trans(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* rank, double* state_p, double* uinv,
    double* ens_p, double* hl_p, double* hxbar_p, double* obs_p,
    double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* type_sqrt, int* type_trans, int* nm1vsn, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdafhyb3dvar_update_estkf(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec, int* dim_cvec_ens,
    double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdafhyb3dvar_update_lestkf(int* step, int* dim_p,
    int* dim_obs_p, int* dim_ens, int* dim_cvec, int* dim_cvec_ens,
    double* state_p, double* ainv, double* ens_p,
    void (*c__init_dim_obs_pdaf)(int* , int* ),
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_pdaf)(int* , int* , double* ),
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* ,
                            double* ),
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_dim_obs_f_pdaf)(int* , int* ),
    void (*c__obs_op_f_pdaf)(int* , int* , int* , double* , double* ),
    void (*c__init_obs_f_pdaf)(int* , int* , double* ),
    void (*c__init_obs_l_pdaf)(int* , int* , int* , double* ),
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    void (*c__init_n_domains_p_pdaf)(int* , int* ),
    void (*c__init_dim_l_pdaf)(int* , int* , int* ),
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ),
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* ,
                              double* ),
    void (*c__g2l_obs_pdaf)(int* , int* , int* , int* , int* , int* ,
                            int* , int* ),
    void (*c__init_obsvar_pdaf)(int* , int* , double* , double* ),
    void (*c__init_obsvar_l_pdaf)(int* , int* , int* , double* , int* ,
                                  double* ),
    int* screen, int* subtype, int* flag) noexcept nogil;

cdef extern void c__pdaf_lestkf_ana_fixed(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, int* rank, double* state_l,
    double* ainv_l, double* ens_l, double* hl_l, double* hxbar_l,
    double* obs_l, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* type_sqrt, int* screen, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_genobs_init(int* subtype, int* param_int,
    int* dim_pint, double* param_real, int* dim_preal,
    bint* ensemblefilter, bint* fixedbasis, int* verbose,
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_genobs_alloc(
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_genobs_config(int* subtype,
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_genobs_set_iparam(int* id, int* value,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_genobs_options() noexcept nogil;

cdef extern void c__pdaf_etkf_ana_t(int* step, int* dim_p, int* dim_obs_p,
    int* dim_ens, double* state_p, double* ainv, double* ens_p,
    double* hz_p, double* hxbar_p, double* obs_p, double* forget,
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* ,
                              double* ),
    int* screen, int* type_trans, int* debug,
    int* flag) noexcept nogil;

cdef extern void c__pdaf_letkf_ana_fixed(int* domain_p, int* step,
    int* dim_l, int* dim_obs_l, int* dim_ens, double* state_l,
    double* ainv_l, double* ens_l, double* hz_l, double* hxbar_l,
    double* obs_l, double* forget,
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* ,
                                double* , double* ),
    int* screen, int* debug, int* flag) noexcept nogil;

