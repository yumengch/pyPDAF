cdef void c__add_obs_err_pdaf (int* step,
                               int* dim_obs_p,
                               double* c_p
                              );
cdef void c__init_ens_pdaf (int* filtertype,
                            int* dim_p,
                            int* dim_ens,
                            double* state_p,
                            double* uinv,
                            double* ens_p,
                            int* flag
                           );
cdef void c__next_observation_pdaf (int* stepnow,
                                    int* nsteps,
                                    int* doexit,
                                    double* time
                                   );
cdef void c__collect_state_pdaf (int* dim_p,
                                 double* state_p
                                );
cdef void c__distribute_state_pdaf (int* dim_p,
                                    double* state_p
                                   );
cdef void c__prepoststep_pdaf (int* step,
                               int* dim_p,
                               int* dim_ens,
                               int* dim_ens_p,
                               int* dim_obs_p,
                               double* state_p,
                               double* uinv,
                               double* ens_p,
                               int* flag
                              );
cdef void c__init_dim_obs_pdaf (int* step,
                                int* dim_obs_p
                               );
cdef void c__init_obs_pdaf (int* step,
                            int* dim_obs_p,
                            double* observation_p
                           );
cdef void c__init_obs_covar_pdaf (int* step,
                                  int* dim_obs,
                                  int* dim_obs_p,
                                  double* covar,
                                  double* obs_p,
                                  bint* isdiag
                                 );
cdef void c__init_obsvar_pdaf (int* step,
                               int* dim_obs_p,
                               double* obs_p,
                               double* meanvar
                              );
cdef void c__prodrinva_pdaf (int* step,
                             int* dim_obs_p,
                             int* rank,
                             double* obs_p,
                             double* a_p,
                             double* c_p
                            );
cdef void c__obs_op_pdaf (int* step,
                          int* dim_p,
                          int* dim_obs_p,
                          double* state_p,
                          double* m_state_p
                         );
cdef void c__g2l_obs_pdaf (int* domain_p,
                           int* step,
                           int* dim_obs_f,
                           int* dim_obs_l,
                           int* mstate_f,
                           int* dim_p,
                           int* mstate_l,
                           int* dim_l
                          );
cdef void c__g2l_state_pdaf (int* step,
                             int* domain_p,
                             int* dim_p,
                             double* state_p,
                             int* dim_l,
                             double* state_l
                            );
cdef void c__init_dim_l_pdaf (int* step,
                              int* domain_p,
                              int* dim_l
                             );
cdef void c__init_dim_obs_f_pdaf (int* step,
                                  int* dim_obs_f
                                 );
cdef void c__init_dim_obs_l_pdaf (int* domain_p,
                                  int* step,
                                  int* dim_obs_f,
                                  int* dim_obs_l
                                 );
cdef void c__init_n_domains_p_pdaf (int* step,
                                    int* n_domains_p
                                   );
cdef void c__init_obs_f_pdaf (int* step,
                              int* dim_obs_f,
                              double* observation_f
                             );
cdef void c__init_obs_l_pdaf (int* domain_p,
                              int* step,
                              int* dim_obs_l,
                              double* observation_l
                             );
cdef void c__init_obsvar_l_pdaf (int* domain_p,
                                 int* step,
                                 int* dim_obs_l,
                                 double* obs_l,
                                 int* dim_obs_p,
                                 double* meanvar_l
                                );
cdef void c__init_obserr_f_pdaf (int* step,
                                 int* dim_obs_f,
                                 double* obs_f,
                                 double* obserr_f
                                );
cdef void c__l2g_state_pdaf (int* step,
                             int* domain_p,
                             int* dim_l,
                             double* state_l,
                             int* dim_p,
                             double* state_p
                            );
cdef void c__obs_op_f_pdaf (int* step,
                            int* dim_p,
                            int* dim_obs_f,
                            double* state_p,
                            double* m_state_f
                           );
cdef void c__prodrinva_l_pdaf (int* domain_p,
                               int* step,
                               int* dim_obs_l,
                               int* rank,
                               double* obs_l,
                               double* a_l,
                               double* c_l
                              );
cdef void c__localize_covar_pdaf (int* dim_p,
                                  int* dim_obs,
                                  double* hp_p,
                                  double* hph
                                 );
cdef void c__likelihood_pdaf (int* step,
                              int* dim_obs_p,
                              double* obs_p,
                              double* resid,
                              double* likely
                             );
cdef void c__likelihood_l_pdaf (int* domain_p,
                                int* step,
                                int* dim_obs_l,
                                double* obs_l,
                                double* resid_l,
                                double* likely_l
                               );
cdef void c__get_obs_f_pdaf (int* step,
                             int* dim_obs_f,
                             double* observation_f
                            );
cdef void c__cvt_adj_ens_pdaf (int* iter,
                               int* dim_p,
                               int* dim_ens,
                               int* dim_cv_ens_p,
                               double* ens_p,
                               double* vcv_p,
                               double* cv_p
                              );
cdef void c__cvt_adj_pdaf (int* iter,
                           int* dim_p,
                           int* dim_cvec,
                           double* vcv_p,
                           double* cv_p
                          );
cdef void c__cvt_pdaf (int* iter,
                       int* dim_p,
                       int* dim_cvec,
                       double* cv_p,
                       double* vv_p
                      );
cdef void c__cvt_ens_pdaf (int* iter,
                           int* dim_p,
                           int* dim_ens,
                           int* dim_cvec_ens,
                           double* ens_p,
                           double* v_p,
                           double* vv_p
                          );
cdef void c__obs_op_adj_pdaf (int* step,
                              int* dim_p,
                              int* dim_obs_p,
                              double* state_p,
                              double* m_state_p
                             );
cdef void c__obs_op_lin_pdaf (int* step,
                              int* dim_p,
                              int* dim_obs_p,
                              double* state_p,
                              double* m_state_p
                             );
cdef void c__dist_stateinc_pdaf (int* dim_p,
                                 double* state_inc_p,
                                 int* first,
                                 int* steps
                                );
cdef void c__prodrinva_hyb_l_pdaf (int* domain_p,
                                   int* step,
                                   int* dim_obs_l,
                                   double* obs_l,
                                   double* resid_l,
                                   double* gamma,
                                   double* likely_l
                                  );
cdef void c__likelihood_hyb_l_pdaf (int* domain_p,
                                    int* step,
                                    int* dim_obs_l,
                                    int* rank,
                                    double* obs_l,
                                    double* gamma,
                                    double* a_l,
                                    double* c_l
                                   );
