from .cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_get_fcst_info(int* steps, double* time, 
    int* doexit) noexcept nogil;

cdef extern void c__pdaf3_assimilate_3dvar_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_en3dvar_estkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_en3dvar_lestkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_hyb3dvar_estkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_hyb3dvar_lestkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_3dvar_all(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_3dvar(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_en3dvar(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_en3dvar_estkf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_en3dvar_lestkf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_hyb3dvar(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_hyb3dvar_estkf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_hyb3dvar_lestkf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_3dvar_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_en3dvar_estkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_en3dvar_lestkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_hyb3dvar_estkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_hyb3dvar_lestkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_local(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* , 
                              double* ), 
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* , 
                              double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_global(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_lenkf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_ensrf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__localize_covar_serial_pdaf)(int* , int* , int* , double* , 
                                          double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_local(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__g2l_state_pdaf)(int* , int* , int* , double* , int* , 
                              double* ), 
    void (*c__l2g_state_pdaf)(int* , int* , int* , double* , int* , 
                              double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_global(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_lenkf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_ensrf(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__localize_covar_serial_pdaf)(int* , int* , int* , double* , 
                                          double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_local_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_global_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_lnetf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* , 
                                 double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_lknetf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__prodrinva_hyb_l_pdaf)(int* , int* , int* , int* , double* , 
                                    double* , double* , double* ), 
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* , 
                                 double* ), 
    void (*c__likelihood_hyb_l_pdaf)(int* , int* , int* , double* , 
                                     double* , double* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_enkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__add_obs_err_pdaf)(int* , int* , double* ), 
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* , 
                                   bint* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_lenkf_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ), 
    void (*c__add_obs_err_pdaf)(int* , int* , double* ), 
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* , 
                                   bint* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assim_offline_nonlin_nondiagr(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_3dvar_all(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_3dvar(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_en3dvar(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_en3dvar_estkf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_en3dvar_lestkf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_hyb3dvar(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_hyb3dvar_estkf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_hyb3dvar_lestkf(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_ens_pdaf)(int* , int* , int* , int* , double* , double* , 
                            double* ), 
    void (*c__cvt_adj_ens_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__cvt_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__cvt_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_lin_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__obs_op_adj_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_local_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_global_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prodrinva_pdaf)(int* , int* , int* , double* , double* , 
                              double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_lnetf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* , 
                                 double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_lknetf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__init_n_domains_p_pdaf)(int* , int* ), 
    void (*c__init_dim_l_pdaf)(int* , int* , int* ), 
    void (*c__init_dim_obs_l_pdaf)(int* , int* , int* , int* ), 
    void (*c__prodrinva_l_pdaf)(int* , int* , int* , int* , double* , 
                                double* , double* ), 
    void (*c__prodrinva_hyb_l_pdaf)(int* , int* , int* , int* , double* , 
                                    double* , double* , double* ), 
    void (*c__likelihood_l_pdaf)(int* , int* , int* , double* , double* , 
                                 double* ), 
    void (*c__likelihood_hyb_l_pdaf)(int* , int* , int* , double* , 
                                     double* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_enkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__add_obs_err_pdaf)(int* , int* , double* ), 
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* , 
                                   bint* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_lenkf_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__localize_covar_pdaf)(int* , int* , double* , double* ), 
    void (*c__add_obs_err_pdaf)(int* , int* , double* ), 
    void (*c__init_obs_covar_pdaf)(int* , int* , int* , double* , double* , 
                                   bint* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_assimilate_nonlin_nondiagr(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__likelihood_pdaf)(int* , int* , double* , double* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_generate_obs(
    void (*c__collect_state_pdaf)(int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__get_obs_f_pdaf)(int* , int* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_generate_obs_offline(
    void (*c__init_dim_obs_pdaf)(int* , int* ), 
    void (*c__obs_op_pdaf)(int* , int* , int* , double* , double* ), 
    void (*c__get_obs_f_pdaf)(int* , int* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

