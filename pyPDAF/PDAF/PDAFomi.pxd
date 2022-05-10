cdef extern void c__init (int* n_obs
                         );
cdef extern void c__pdafomi_set_doassim (int* i_obs,
                                         int* doassim
                                        );
cdef extern void c__pdafomi_set_disttype (int* i_obs,
                                          int* disttype
                                         );
cdef extern void c__pdafomi_set_ncoord (int* i_obs,
                                        int* ncoord
                                       );
cdef extern void c__pdafomi_set_id_obs_p (int* i_obs,
                                          int* nrows,
                                          int* dim_obs_p,
                                          int* id_obs_p
                                         );
cdef extern void c__pdafomi_set_icoeff_p (int* i_obs,
                                          int* nrows,
                                          int* dim_obs_p,
                                          double* icoeff_p
                                         );
cdef extern void c__pdafomi_set_domainsize (int* i_obs,
                                            int* ncoord,
                                            double* domainsize
                                           );
cdef extern void c__pdafomi_set_obs_err_type (int* i_obs,
                                              int* obs_err_type
                                             );
cdef extern void c__pdafomi_set_use_global_obs (int* i_obs,
                                                int* use_global_obs
                                               );
cdef extern void c__pdafomi_gather_obs (int* i_obs,
                                        int* dim_obs_p,
                                        double* obs_p,
                                        double* ivar_obs_p,
                                        double* ocoord_p,
                                        double* local_range,
                                        int* dim_obs
                                       );
cdef extern void c__pdafomi_gather_obsstate (int* i_obs,
                                             double* obsstate_p,
                                             double* obsstate_f,
                                             int* nobs_f_all
                                            );
cdef extern void c__pdafomi_localize_covar (int* i_obs,
                                            int* dim_p,
                                            int* dim_obs,
                                            int* dim_coords,
                                            int* locweight,
                                            double* local_range,
                                            double* srange,
                                            double* coords_p,
                                            double* hp_p,
                                            double* hph
                                           );
cdef extern void c__pdafomi_set_domain_limits (double* lim_coords
                                              );
cdef extern void c__pdafomi_set_debug_flag (int* debugval
                                           );
cdef extern void c__pdafomi_deallocate_obs (int* i_obs
                                           );
cdef extern void c__pdafomi_init_dim_obs_l (int* i_obs,
                                            double* coords_l,
                                            int* locweight,
                                            double* local_range,
                                            double* srange,
                                            int* dim_obs_l
                                           );
cdef extern void c__pdafomi_obs_op_gridpoint (int* i_obs,
                                              double* state_p,
                                              int* dim_p,
                                              double* obs_f_all,
                                              int* nobs_f_all
                                             );
cdef extern void c__pdafomi_obs_op_gridavg (int* i_obs,
                                            int* nrows,
                                            double* state_p,
                                            int* dim_p,
                                            double* obs_f_all,
                                            int* nobs_f_all
                                           );
cdef extern void c__pdafomi_obs_op_interp_lin (int* i_obs,
                                               int* nrows,
                                               double* state_p,
                                               int* dim_p,
                                               double* obs_f_all,
                                               int* nobs_f_all
                                              );
cdef extern void c__pdafomi_obs_op_adj_gridavg (int* i_obs,
                                                int* nrows,
                                                double* state_p,
                                                int* dim_p,
                                                double* obs_f_all,
                                                int* nobs_f_all
                                               );
cdef extern void c__pdafomi_obs_op_adj_gridpoint (int* i_obs,
                                                  double* state_p,
                                                  int* dim_p,
                                                  double* obs_f_all,
                                                  int* nobs_f_all
                                                 );
cdef extern void c__pdafomi_obs_op_adj_interp_lin (int* i_obs,
                                                   int* nrows,
                                                   double* state_p,
                                                   int* dim_p,
                                                   double* obs_f_all,
                                                   int* nobs_f_all
                                                  );
cdef extern void c__pdafomi_get_interp_coeff_tri (double* gpc,
                                                  double* oc,
                                                  double* icoeff
                                                 );
cdef extern void c__pdafomi_get_interp_coeff_lin1d (double* gpc,
                                                    double* oc,
                                                    double* icoeff
                                                   );
cdef extern void c__pdafomi_get_interp_coeff_lin (int* num_gp,
                                                  int* n_dim,
                                                  double* gpc,
                                                  double* oc,
                                                  double* icoeff
                                                 );
cdef extern void c__pdafomi_assimilate_3dvar (void (*c__collect_state_pdaf)(int*,
                                                                            double*
                                                                           ),
                                              void (*c__distribute_state_pdaf)(int*,
                                                                               double*
                                                                              ),
                                              void (*c__init_dim_obs_pdaf)(int*,
                                                                           int*
                                                                          ),
                                              void (*c__obs_op_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                              void (*c__cvt_pdaf)(int*,
                                                                  int*,
                                                                  int*,
                                                                  double*,
                                                                  double*
                                                                 ),
                                              void (*c__cvt_adj_pdaf)(int*,
                                                                      int*,
                                                                      int*,
                                                                      double*,
                                                                      double*
                                                                     ),
                                              void (*c__obs_op_lin_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*
                                                                        ),
                                              void (*c__obs_op_adj_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*
                                                                        ),
                                              void (*c__prepoststep_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*,
                                                                          double*,
                                                                          int*
                                                                         ),
                                              void (*c__next_observation_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                              int* outflag
                                             );
cdef extern void c__pdafomi_assimilate_en3dvar_estkf (void (*c__collect_state_pdaf)(int*,
                                                                                    double*
                                                                                   ),
                                                      void (*c__distribute_state_pdaf)(int*,
                                                                                       double*
                                                                                      ),
                                                      void (*c__init_dim_obs_pdaf)(int*,
                                                                                   int*
                                                                                  ),
                                                      void (*c__obs_op_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                      void (*c__cvt_ens_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                      void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                      void (*c__obs_op_lin_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__obs_op_adj_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__prepoststep_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*,
                                                                                  int*
                                                                                 ),
                                                      void (*c__next_observation_pdaf)(int*,
                                                                                       int*,
                                                                                       int*,
                                                                                       double*
                                                                                      ),
                                                      int* outflag
                                                     );
cdef extern void c__pdafomi_assimilate_en3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
                                                                                     double*
                                                                                    ),
                                                       void (*c__distribute_state_pdaf)(int*,
                                                                                        double*
                                                                                       ),
                                                       void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                      int*
                                                                                     ),
                                                       void (*c__obs_op_f_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                       void (*c__cvt_ens_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                       void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*
                                                                                  ),
                                                       void (*c__obs_op_lin_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__obs_op_adj_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__init_n_domains_p_pdaf)(int*,
                                                                                        int*
                                                                                       ),
                                                       void (*c__init_dim_l_pdaf)(int*,
                                                                                  int*,
                                                                                  int*
                                                                                 ),
                                                       void (*c__init_dim_obs_l_pdaf)(int*,
                                                                                      int*,
                                                                                      int*,
                                                                                      int*
                                                                                     ),
                                                       void (*c__g2l_state_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 int*,
                                                                                 double*
                                                                                ),
                                                       void (*c__l2g_state_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 int*,
                                                                                 double*
                                                                                ),
                                                       void (*c__prepoststep_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*,
                                                                                   int*
                                                                                  ),
                                                       void (*c__next_observation_pdaf)(int*,
                                                                                        int*,
                                                                                        int*,
                                                                                        double*
                                                                                       ),
                                                       int* outflag
                                                      );
cdef extern void c__pdafomi_assimilate_global (void (*c__collect_state_pdaf)(int*,
                                                                             double*
                                                                            ),
                                               void (*c__distribute_state_pdaf)(int*,
                                                                                double*
                                                                               ),
                                               void (*c__init_dim_obs_pdaf)(int*,
                                                                            int*
                                                                           ),
                                               void (*c__obs_op_pdaf)(int*,
                                                                      int*,
                                                                      int*,
                                                                      double*,
                                                                      double*
                                                                     ),
                                               void (*c__prepoststep_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           double*,
                                                                           double*,
                                                                           double*,
                                                                           int*
                                                                          ),
                                               void (*c__next_observation_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*
                                                                               ),
                                               int* flag
                                              );
cdef extern void c__pdafomi_assimilate_hyb3dvar_estkf (void (*c__collect_state_pdaf)(int*,
                                                                                     double*
                                                                                    ),
                                                       void (*c__distribute_state_pdaf)(int*,
                                                                                        double*
                                                                                       ),
                                                       void (*c__init_dim_obs_pdaf)(int*,
                                                                                    int*
                                                                                   ),
                                                       void (*c__obs_op_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                       void (*c__cvt_ens_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                       void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*
                                                                                  ),
                                                       void (*c__cvt_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*,
                                                                           double*
                                                                          ),
                                                       void (*c__cvt_adj_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                       void (*c__obs_op_lin_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__obs_op_adj_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__prepoststep_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*,
                                                                                   int*
                                                                                  ),
                                                       void (*c__next_observation_pdaf)(int*,
                                                                                        int*,
                                                                                        int*,
                                                                                        double*
                                                                                       ),
                                                       int* outflag
                                                      );
cdef extern void c__pdafomi_assimilate_hyb3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
                                                                                      double*
                                                                                     ),
                                                        void (*c__distribute_state_pdaf)(int*,
                                                                                         double*
                                                                                        ),
                                                        void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                       int*
                                                                                      ),
                                                        void (*c__obs_op_f_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                        void (*c__cvt_ens_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                        void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                    int*,
                                                                                    int*,
                                                                                    int*,
                                                                                    double*,
                                                                                    double*,
                                                                                    double*
                                                                                   ),
                                                        void (*c__cvt_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
                                                                            double*
                                                                           ),
                                                        void (*c__cvt_adj_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                        void (*c__obs_op_lin_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*
                                                                                  ),
                                                        void (*c__obs_op_adj_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*
                                                                                  ),
                                                        void (*c__init_n_domains_p_pdaf)(int*,
                                                                                         int*
                                                                                        ),
                                                        void (*c__init_dim_l_pdaf)(int*,
                                                                                   int*,
                                                                                   int*
                                                                                  ),
                                                        void (*c__init_dim_obs_l_pdaf)(int*,
                                                                                       int*,
                                                                                       int*,
                                                                                       int*
                                                                                      ),
                                                        void (*c__g2l_state_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  int*,
                                                                                  double*
                                                                                 ),
                                                        void (*c__l2g_state_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  int*,
                                                                                  double*
                                                                                 ),
                                                        void (*c__prepoststep_pdaf)(int*,
                                                                                    int*,
                                                                                    int*,
                                                                                    int*,
                                                                                    int*,
                                                                                    double*,
                                                                                    double*,
                                                                                    double*,
                                                                                    int*
                                                                                   ),
                                                        void (*c__next_observation_pdaf)(int*,
                                                                                         int*,
                                                                                         int*,
                                                                                         double*
                                                                                        ),
                                                        int* outflag
                                                       );
cdef extern void c__pdafomi_assimilate_lenkf (void (*c__collect_state_pdaf)(int*,
                                                                            double*
                                                                           ),
                                              void (*c__distribute_state_pdaf)(int*,
                                                                               double*
                                                                              ),
                                              void (*c__init_dim_obs_pdaf)(int*,
                                                                           int*
                                                                          ),
                                              void (*c__obs_op_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                              void (*c__prepoststep_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*,
                                                                          double*,
                                                                          int*
                                                                         ),
                                              void (*c__localize_covar_pdaf)(int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                              void (*c__next_observation_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                              int* flag
                                             );
cdef extern void c__pdafomi_assimilate_local (void (*c__collect_state_pdaf)(int*,
                                                                            double*
                                                                           ),
                                              void (*c__distribute_state_pdaf)(int*,
                                                                               double*
                                                                              ),
                                              void (*c__init_dim_obs_pdaf)(int*,
                                                                           int*
                                                                          ),
                                              void (*c__obs_op_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                              void (*c__prepoststep_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*,
                                                                          double*,
                                                                          int*
                                                                         ),
                                              void (*c__init_n_domains_p_pdaf)(int*,
                                                                               int*
                                                                              ),
                                              void (*c__init_dim_l_pdaf)(int*,
                                                                         int*,
                                                                         int*
                                                                        ),
                                              void (*c__init_dim_obs_l_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             int*
                                                                            ),
                                              void (*c__g2l_state_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        int*,
                                                                        double*
                                                                       ),
                                              void (*c__l2g_state_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        int*,
                                                                        double*
                                                                       ),
                                              void (*c__next_observation_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                              int* flag
                                             );
cdef extern void c__pdafomi_generate_obs (void (*c__collect_state_pdaf)(int*,
                                                                        double*
                                                                       ),
                                          void (*c__distribute_state_pdaf)(int*,
                                                                           double*
                                                                          ),
                                          void (*c__init_dim_obs_f_pdaf)(int*,
                                                                         int*
                                                                        ),
                                          void (*c__obs_op_f_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   double*,
                                                                   double*
                                                                  ),
                                          void (*c__get_obs_f_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                          void (*c__prepoststep_pdaf)(int*,
                                                                      int*,
                                                                      int*,
                                                                      int*,
                                                                      int*,
                                                                      double*,
                                                                      double*,
                                                                      double*,
                                                                      int*
                                                                     ),
                                          void (*c__next_observation_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                          int* flag
                                         );
cdef extern void c__pdafomi_put_state_3dvar (void (*c__collect_state_pdaf)(int*,
                                                                           double*
                                                                          ),
                                             void (*c__init_dim_obs_pdaf)(int*,
                                                                          int*
                                                                         ),
                                             void (*c__obs_op_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                             void (*c__cvt_pdaf)(int*,
                                                                 int*,
                                                                 int*,
                                                                 double*,
                                                                 double*
                                                                ),
                                             void (*c__cvt_adj_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                             void (*c__obs_op_lin_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                             void (*c__obs_op_adj_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                             void (*c__prepoststep_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         double*,
                                                                         int*
                                                                        ),
                                             int* outflag
                                            );
cdef extern void c__pdafomi_put_state_en3dvar_estkf (void (*c__collect_state_pdaf)(int*,
                                                                                   double*
                                                                                  ),
                                                     void (*c__init_dim_obs_pdaf)(int*,
                                                                                  int*
                                                                                 ),
                                                     void (*c__obs_op_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
                                                                            double*
                                                                           ),
                                                     void (*c__cvt_ens_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                     void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                     void (*c__obs_op_lin_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                     void (*c__obs_op_adj_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                     void (*c__prepoststep_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*,
                                                                                 double*,
                                                                                 int*
                                                                                ),
                                                     int* outflag
                                                    );
cdef extern void c__pdafomi_put_state_en3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
                                                                                    double*
                                                                                   ),
                                                      void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                     int*
                                                                                    ),
                                                      void (*c__obs_op_f_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                      void (*c__cvt_ens_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                      void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                      void (*c__obs_op_lin_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__obs_op_adj_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__init_n_domains_p_pdaf)(int*,
                                                                                       int*
                                                                                      ),
                                                      void (*c__init_dim_l_pdaf)(int*,
                                                                                 int*,
                                                                                 int*
                                                                                ),
                                                      void (*c__init_dim_obs_l_pdaf)(int*,
                                                                                     int*,
                                                                                     int*,
                                                                                     int*
                                                                                    ),
                                                      void (*c__g2l_state_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                int*,
                                                                                double*
                                                                               ),
                                                      void (*c__l2g_state_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                int*,
                                                                                double*
                                                                               ),
                                                      void (*c__prepoststep_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*,
                                                                                  int*
                                                                                 ),
                                                      int* outflag
                                                     );
cdef extern void c__pdafomi_put_state_generate_obs (void (*c__collect_state_pdaf)(int*,
                                                                                  double*
                                                                                 ),
                                                    void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                   int*
                                                                                  ),
                                                    void (*c__obs_op_f_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                    void (*c__get_obs_f_pdaf)(int*,
                                                                              int*,
                                                                              double*
                                                                             ),
                                                    void (*c__prepoststep_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*,
                                                                                double*,
                                                                                int*
                                                                               ),
                                                    int* flag
                                                   );
cdef extern void c__pdafomi_put_state_global (void (*c__collect_state_pdaf)(int*,
                                                                            double*
                                                                           ),
                                              void (*c__init_dim_obs_pdaf)(int*,
                                                                           int*
                                                                          ),
                                              void (*c__obs_op_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                              void (*c__prepoststep_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*,
                                                                          double*,
                                                                          int*
                                                                         ),
                                              int* flag
                                             );
cdef extern void c__pdafomi_put_state_hyb3dvar_estkf (void (*c__collect_state_pdaf)(int*,
                                                                                    double*
                                                                                   ),
                                                      void (*c__init_dim_obs_pdaf)(int*,
                                                                                   int*
                                                                                  ),
                                                      void (*c__obs_op_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                      void (*c__cvt_ens_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                      void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                      void (*c__cvt_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*
                                                                         ),
                                                      void (*c__cvt_adj_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                      void (*c__obs_op_lin_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__obs_op_adj_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                      void (*c__prepoststep_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*,
                                                                                  double*,
                                                                                  int*
                                                                                 ),
                                                      int* outflag
                                                     );
cdef extern void c__pdafomi_put_state_hyb3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
                                                                                     double*
                                                                                    ),
                                                       void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                      int*
                                                                                     ),
                                                       void (*c__obs_op_f_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                       void (*c__cvt_ens_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                       void (*c__cvt_adj_ens_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*
                                                                                  ),
                                                       void (*c__cvt_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*,
                                                                           double*
                                                                          ),
                                                       void (*c__cvt_adj_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                       void (*c__obs_op_lin_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__obs_op_adj_pdaf)(int*,
                                                                                  int*,
                                                                                  int*,
                                                                                  double*,
                                                                                  double*
                                                                                 ),
                                                       void (*c__init_n_domains_p_pdaf)(int*,
                                                                                        int*
                                                                                       ),
                                                       void (*c__init_dim_l_pdaf)(int*,
                                                                                  int*,
                                                                                  int*
                                                                                 ),
                                                       void (*c__init_dim_obs_l_pdaf)(int*,
                                                                                      int*,
                                                                                      int*,
                                                                                      int*
                                                                                     ),
                                                       void (*c__g2l_state_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 int*,
                                                                                 double*
                                                                                ),
                                                       void (*c__l2g_state_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 int*,
                                                                                 double*
                                                                                ),
                                                       void (*c__prepoststep_pdaf)(int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   int*,
                                                                                   double*,
                                                                                   double*,
                                                                                   double*,
                                                                                   int*
                                                                                  ),
                                                       int* outflag
                                                      );
cdef extern void c__pdafomi_put_state_lenkf (void (*c__collect_state_pdaf)(int*,
                                                                           double*
                                                                          ),
                                             void (*c__init_dim_obs_pdaf)(int*,
                                                                          int*
                                                                         ),
                                             void (*c__obs_op_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                             void (*c__prepoststep_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         double*,
                                                                         int*
                                                                        ),
                                             void (*c__localize_covar_pdaf)(int*,
                                                                            int*,
                                                                            double*,
                                                                            double*
                                                                           ),
                                             int* flag
                                            );
cdef extern void c__pdafomi_put_state_local (void (*c__collect_state_pdaf)(int*,
                                                                           double*
                                                                          ),
                                             void (*c__init_dim_obs_pdaf)(int*,
                                                                          int*
                                                                         ),
                                             void (*c__obs_op_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                             void (*c__prepoststep_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         double*,
                                                                         int*
                                                                        ),
                                             void (*c__init_n_domains_p_pdaf)(int*,
                                                                              int*
                                                                             ),
                                             void (*c__init_dim_l_pdaf)(int*,
                                                                        int*,
                                                                        int*
                                                                       ),
                                             void (*c__init_dim_obs_l_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            int*
                                                                           ),
                                             void (*c__g2l_state_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
                                                                       int*,
                                                                       double*
                                                                      ),
                                             void (*c__l2g_state_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
                                                                       int*,
                                                                       double*
                                                                      ),
                                             int* flag
                                            );
