cdef extern void c__pdaf_assimilate_3dvar (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                           void (*c__prodrinva_pdaf)(int*,
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
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_en3dvar_estkf (void (*c__collect_state_pdaf)(int*,
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
                                                   void (*c__init_obs_pdaf)(int*,
                                                                            int*,
                                                                            double*
                                                                           ),
                                                   void (*c__prodrinva_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
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
                                                   void (*c__init_obsvar_pdaf)(int*,
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
                                                  ) noexcept;
cdef extern void c__pdaf_assimilate_en3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                                    void (*c__init_obs_pdaf)(int*,
                                                                             int*,
                                                                             double*
                                                                            ),
                                                    void (*c__prodrinva_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
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
                                                    void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                   int*
                                                                                  ),
                                                    void (*c__obs_op_f_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                    void (*c__init_obs_f_pdaf)(int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                                    void (*c__init_obs_l_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                                    void (*c__prodrinva_l_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
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
                                                    void (*c__g2l_obs_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*
                                                                           ),
                                                    void (*c__init_obsvar_pdaf)(int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                    void (*c__init_obsvar_l_pdaf)(int*,
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
                                                   ) noexcept;
cdef extern void c__pdaf_assimilate_enkf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__add_obs_err_pdaf)(int*,
                                                                      int*,
                                                                      double*
                                                                     ),
                                          void (*c__init_obs_covar_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         bint*
                                                                        ),
                                          void (*c__next_observation_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_assimilate_estkf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
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
                                           void (*c__prodrinva_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                           void (*c__init_obsvar_pdaf)(int*,
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
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_etkf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__prodrinva_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                          void (*c__init_obsvar_pdaf)(int*,
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
                                         ) noexcept;
cdef extern void c__pdaf_assimilate_hyb3dvar_estkf (void (*c__collect_state_pdaf)(int*,
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
                                                    void (*c__init_obs_pdaf)(int*,
                                                                             int*,
                                                                             double*
                                                                            ),
                                                    void (*c__prodrinva_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
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
                                                    void (*c__init_obsvar_pdaf)(int*,
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
                                                   ) noexcept;
cdef extern void c__pdaf_assimilate_hyb3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                                     void (*c__init_obs_pdaf)(int*,
                                                                              int*,
                                                                              double*
                                                                             ),
                                                     void (*c__prodrinva_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
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
                                                     void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                    int*
                                                                                   ),
                                                     void (*c__obs_op_f_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
                                                                              double*
                                                                             ),
                                                     void (*c__init_obs_f_pdaf)(int*,
                                                                                int*,
                                                                                double*
                                                                               ),
                                                     void (*c__init_obs_l_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                double*
                                                                               ),
                                                     void (*c__prodrinva_l_pdaf)(int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 int*,
                                                                                 double*,
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
                                                     void (*c__g2l_obs_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             int*
                                                                            ),
                                                     void (*c__init_obsvar_pdaf)(int*,
                                                                                 int*,
                                                                                 double*,
                                                                                 double*
                                                                                ),
                                                     void (*c__init_obsvar_l_pdaf)(int*,
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
                                                    ) noexcept;
cdef extern void c__pdaf_assimilate_lenkf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
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
                                           void (*c__localize_covar_pdaf)(int*,
                                                                          int*,
                                                                          double*,
                                                                          double*
                                                                         ),
                                           void (*c__add_obs_err_pdaf)(int*,
                                                                       int*,
                                                                       double*
                                                                      ),
                                           void (*c__init_obs_covar_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          double*,
                                                                          bint*
                                                                         ),
                                           void (*c__next_observation_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*
                                                                           ),
                                           int* flag
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                            void (*c__init_obs_pdaf)(int*,
                                                                     int*,
                                                                     double*
                                                                    ),
                                            void (*c__init_obs_l_pdaf)(int*,
                                                                       int*,
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
                                            void (*c__prodrinva_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
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
                                            void (*c__g2l_obs_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*
                                                                   ),
                                            void (*c__init_obsvar_pdaf)(int*,
                                                                        int*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                            void (*c__init_obsvar_l_pdaf)(int*,
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
                                           ) noexcept;
cdef extern void c__pdaf_assimilate_letkf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                           void (*c__init_obs_l_pdaf)(int*,
                                                                      int*,
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
                                           void (*c__prodrinva_l_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
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
                                           void (*c__g2l_obs_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*
                                                                  ),
                                           void (*c__init_obsvar_pdaf)(int*,
                                                                       int*,
                                                                       double*,
                                                                       double*
                                                                      ),
                                           void (*c__init_obsvar_l_pdaf)(int*,
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
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_lnetf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_l_pdaf)(int*,
                                                                      int*,
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
                                           void (*c__likelihood_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
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
                                           void (*c__g2l_obs_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*
                                                                  ),
                                           void (*c__next_observation_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*
                                                                           ),
                                           int* flag
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_lknetf (void (*c__collect_state_pdaf)(int*,
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
                                            void (*c__init_obs_pdaf)(int*,
                                                                     int*,
                                                                     double*
                                                                    ),
                                            void (*c__init_obs_l_pdaf)(int*,
                                                                       int*,
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
                                            void (*c__prodrinva_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                            void (*c__prodrinva_hyb_l_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
                                                                            double*,
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
                                            void (*c__g2l_obs_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*,
                                                                    int*
                                                                   ),
                                            void (*c__init_obsvar_pdaf)(int*,
                                                                        int*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                            void (*c__init_obsvar_l_pdaf)(int*,
                                                                          int*,
                                                                          int*,
                                                                          double*,
                                                                          int*,
                                                                          double*
                                                                         ),
                                            void (*c__likelihood_l_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         double*
                                                                        ),
                                            void (*c__likelihood_hyb_l_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                            void (*c__next_observation_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*
                                                                            ),
                                            int* flag
                                           ) noexcept;
cdef extern void c__pdaf_assimilate_lseik (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                           void (*c__init_obs_l_pdaf)(int*,
                                                                      int*,
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
                                           void (*c__prodrinva_l_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
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
                                           void (*c__g2l_obs_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*
                                                                  ),
                                           void (*c__init_obsvar_pdaf)(int*,
                                                                       int*,
                                                                       double*,
                                                                       double*
                                                                      ),
                                           void (*c__init_obsvar_l_pdaf)(int*,
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
                                          ) noexcept;
cdef extern void c__pdaf_assimilate_netf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__likelihood_pdaf)(int*,
                                                                     int*,
                                                                     double*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                          void (*c__next_observation_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_assimilate_pf (void (*c__collect_state_pdaf)(int*,
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
                                        void (*c__init_obs_pdaf)(int*,
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
                                        void (*c__likelihood_pdaf)(int*,
                                                                   int*,
                                                                   double*,
                                                                   double*,
                                                                   double*
                                                                  ),
                                        void (*c__next_observation_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*
                                                                        ),
                                        int* flag
                                       ) noexcept;
cdef extern void c__pdaf_assimilate_seek (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__prodrinva_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                          void (*c__next_observation_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_assimilate_seik (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__prodrinva_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                          void (*c__next_observation_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_assimilate_prepost (void (*c__collect_state_pdaf)(int*,
                                                                           double*
                                                                          ),
                                             void (*c__distribute_state_pdaf)(int*,
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
                                            ) noexcept;
cdef extern void c__pdaf_deallocate () noexcept;
cdef extern void c__pdaf_diag_effsample (int* dim_sample,
                                         double* weights,
                                         double* effsample
                                        ) noexcept;
cdef extern void c__pdaf_diag_ensstats (int* dim,
                                        int* dim_ens,
                                        int* element,
                                        double* state,
                                        double* ens,
                                        double* skewness,
                                        double* kurtosis,
                                        int* status
                                       ) noexcept;
cdef extern void c__pdaf_diag_histogram (int* ncall,
                                         int* dim,
                                         int* dim_ens,
                                         int* element,
                                         double* state,
                                         double* ens,
                                         int* hist,
                                         double* delta,
                                         int* status
                                        ) noexcept;
cdef extern void c__pdaf_eofcovar (int* dim_state,
                                   int* nstates,
                                   int* nfields,
                                   int* dim_fields,
                                   int* offsets,
                                   int* remove_mstate,
                                   int* do_mv,
                                   double* states,
                                   double* stddev,
                                   double* svals,
                                   double* svec,
                                   double* meanstate,
                                   int* verbose,
                                   int* status
                                  ) noexcept;
cdef extern void c__pdaf_gather_dim_obs_f (int* dim_obs_p,
                                           int* dim_obs_f
                                          ) noexcept;
cdef extern void c__pdaf_gather_obs_f (double* obs_p,
                                       int* dimobs_p,
                                       double* obs_f,
                                       int* dimobs_f,
                                       int* status
                                      ) noexcept;
cdef extern void c__pdaf_gather_obs_f2 (double* coords_p,
                                        int* dimobs_p,
                                        double* coords_f,
                                        int* dimobs_f,
                                        int* nrows,
                                        int* status
                                       ) noexcept;
cdef extern void c__pdaf_generate_obs (void (*c__collect_state_pdaf)(int*,
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
                                       void (*c__init_obserr_f_pdaf)(int*,
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
                                      ) noexcept;
cdef extern void c__pdaf_get_assim_flag (int* did_assim
                                        ) noexcept;
cdef extern void c__pdaf_get_ensstats (int* dims,
                                       double** c_skew_ptr,
                                       double** c_kurt_ptr,
                                       int* status
                                      ) noexcept;
cdef extern void c__pdaf_get_localfilter (int* lfilter
                                         ) noexcept;
cdef extern void c__pdaf_get_memberid (int* memberid
                                      ) noexcept;
cdef extern void c__pdaf_get_obsmemberid (int* memberid
                                         ) noexcept;
cdef extern void c__pdaf_get_smootherens (double** c_sens_point,
                                          int* maxlag,
                                          int* dims,
                                          int* status
                                         ) noexcept;
cdef extern void c__pdaf_get_state (int* steps,
                                    double* time,
                                    int* doexit,
                                    void (*c__next_observation_pdaf)(int*,
                                                                     int*,
                                                                     int*,
                                                                     double*
                                                                    ),
                                    void (*c__distribute_state_pdaf)(int*,
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
                                   ) noexcept;
cdef extern void c__pdaf_init (int* filtertype,
                               int* subtype,
                               int* stepnull,
                               int* param_int,
                               int* dim_pint,
                               double* param_real,
                               int* dim_preal,
                               int* comm_model,
                               int* comm_filter,
                               int* comm_couple,
                               int* task_id,
                               int* n_modeltasks,
                               bint* in_filterpe,
                               void (*c__init_ens_pdaf)(int*,
                                                        int*,
                                                        int*,
                                                        double*,
                                                        double*,
                                                        double*,
                                                        int*
                                                       ),
                               int* in_screen,
                               int* flag
                              ) noexcept;
cdef extern void c__pdaf_local_weight (int* wtype,
                                       int* rtype,
                                       double* cradius,
                                       double* sradius,
                                       double* distance,
                                       int* nrows,
                                       int* ncols,
                                       double* a,
                                       double* var_obs,
                                       double* weight,
                                       int* verbose
                                      ) noexcept;
cdef extern void c__pdaf_print_info (int* printtype
                                    ) noexcept;
cdef extern void c__pdaf_put_state_3dvar (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
                                                                   int*,
                                                                   double*
                                                                  ),
                                          void (*c__prodrinva_pdaf)(int*,
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
                                         ) noexcept;
cdef extern void c__pdaf_put_state_en3dvar_estkf (void (*c__collect_state_pdaf)(int*,
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
                                                  void (*c__init_obs_pdaf)(int*,
                                                                           int*,
                                                                           double*
                                                                          ),
                                                  void (*c__prodrinva_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
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
                                                  void (*c__init_obsvar_pdaf)(int*,
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
                                                 ) noexcept;
cdef extern void c__pdaf_put_state_en3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                                   void (*c__init_obs_pdaf)(int*,
                                                                            int*,
                                                                            double*
                                                                           ),
                                                   void (*c__prodrinva_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
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
                                                   void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                  int*
                                                                                 ),
                                                   void (*c__obs_op_f_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
                                                                            double*
                                                                           ),
                                                   void (*c__init_obs_f_pdaf)(int*,
                                                                              int*,
                                                                              double*
                                                                             ),
                                                   void (*c__init_obs_l_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*
                                                                             ),
                                                   void (*c__prodrinva_l_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               int*,
                                                                               double*,
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
                                                   void (*c__g2l_obs_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           int*,
                                                                           int*
                                                                          ),
                                                   void (*c__init_obsvar_pdaf)(int*,
                                                                               int*,
                                                                               double*,
                                                                               double*
                                                                              ),
                                                   void (*c__init_obsvar_l_pdaf)(int*,
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
                                                  ) noexcept;
cdef extern void c__pdaf_put_state_enkf (void (*c__collect_state_pdaf)(int*,
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
                                         void (*c__init_obs_pdaf)(int*,
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
                                         void (*c__add_obs_err_pdaf)(int*,
                                                                     int*,
                                                                     double*
                                                                    ),
                                         void (*c__init_obs_covar_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        double*,
                                                                        bint*
                                                                       ),
                                         int* flag
                                        ) noexcept;
cdef extern void c__pdaf_put_state_estkf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__prodrinva_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                          void (*c__init_obsvar_pdaf)(int*,
                                                                      int*,
                                                                      double*,
                                                                      double*
                                                                     ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_put_state_etkf (void (*c__collect_state_pdaf)(int*,
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
                                         void (*c__init_obs_pdaf)(int*,
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
                                         void (*c__prodrinva_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   double*,
                                                                   double*,
                                                                   double*
                                                                  ),
                                         void (*c__init_obsvar_pdaf)(int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                         int* flag
                                        ) noexcept;
cdef extern void c__pdaf_put_state_generate_obs (void (*c__collect_state_pdaf)(int*,
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
                                                 void (*c__init_obserr_f_pdaf)(int*,
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
                                                ) noexcept;
cdef extern void c__pdaf_put_state_hyb3dvar_estkf (void (*c__collect_state_pdaf)(int*,
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
                                                   void (*c__init_obs_pdaf)(int*,
                                                                            int*,
                                                                            double*
                                                                           ),
                                                   void (*c__prodrinva_pdaf)(int*,
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
                                                   void (*c__init_obsvar_pdaf)(int*,
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
                                                  ) noexcept;
cdef extern void c__pdaf_put_state_hyb3dvar_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                                    void (*c__init_obs_pdaf)(int*,
                                                                             int*,
                                                                             double*
                                                                            ),
                                                    void (*c__prodrinva_pdaf)(int*,
                                                                              int*,
                                                                              int*,
                                                                              double*,
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
                                                    void (*c__init_dim_obs_f_pdaf)(int*,
                                                                                   int*
                                                                                  ),
                                                    void (*c__obs_op_f_pdaf)(int*,
                                                                             int*,
                                                                             int*,
                                                                             double*,
                                                                             double*
                                                                            ),
                                                    void (*c__init_obs_f_pdaf)(int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                                    void (*c__init_obs_l_pdaf)(int*,
                                                                               int*,
                                                                               int*,
                                                                               double*
                                                                              ),
                                                    void (*c__prodrinva_l_pdaf)(int*,
                                                                                int*,
                                                                                int*,
                                                                                int*,
                                                                                double*,
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
                                                    void (*c__g2l_obs_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            int*
                                                                           ),
                                                    void (*c__init_obsvar_pdaf)(int*,
                                                                                int*,
                                                                                double*,
                                                                                double*
                                                                               ),
                                                    void (*c__init_obsvar_l_pdaf)(int*,
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
                                                   ) noexcept;
cdef extern void c__pdaf_put_state_lenkf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
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
                                          void (*c__localize_covar_pdaf)(int*,
                                                                         int*,
                                                                         double*,
                                                                         double*
                                                                        ),
                                          void (*c__add_obs_err_pdaf)(int*,
                                                                      int*,
                                                                      double*
                                                                     ),
                                          void (*c__init_obs_covar_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         double*,
                                                                         bint*
                                                                        ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_put_state_lestkf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                           void (*c__init_obs_l_pdaf)(int*,
                                                                      int*,
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
                                           void (*c__prodrinva_l_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
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
                                           void (*c__g2l_obs_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*
                                                                  ),
                                           void (*c__init_obsvar_pdaf)(int*,
                                                                       int*,
                                                                       double*,
                                                                       double*
                                                                      ),
                                           void (*c__init_obsvar_l_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         int*,
                                                                         double*
                                                                        ),
                                           int* flag
                                          ) noexcept;
cdef extern void c__pdaf_put_state_letkf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
                                                                   int*,
                                                                   double*
                                                                  ),
                                          void (*c__init_obs_l_pdaf)(int*,
                                                                     int*,
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
                                          void (*c__prodrinva_l_pdaf)(int*,
                                                                      int*,
                                                                      int*,
                                                                      int*,
                                                                      double*,
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
                                          void (*c__g2l_obs_pdaf)(int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*
                                                                 ),
                                          void (*c__init_obsvar_pdaf)(int*,
                                                                      int*,
                                                                      double*,
                                                                      double*
                                                                     ),
                                          void (*c__init_obsvar_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        int*,
                                                                        double*
                                                                       ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_put_state_lnetf (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_l_pdaf)(int*,
                                                                     int*,
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
                                          void (*c__likelihood_l_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
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
                                          void (*c__g2l_obs_pdaf)(int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*
                                                                 ),
                                          int* outflag
                                         ) noexcept;
cdef extern void c__pdaf_put_state_lknetf (void (*c__collect_state_pdaf)(int*,
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
                                           void (*c__init_obs_pdaf)(int*,
                                                                    int*,
                                                                    double*
                                                                   ),
                                           void (*c__init_obs_l_pdaf)(int*,
                                                                      int*,
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
                                           void (*c__prodrinva_l_pdaf)(int*,
                                                                       int*,
                                                                       int*,
                                                                       int*,
                                                                       double*,
                                                                       double*,
                                                                       double*
                                                                      ),
                                           void (*c__prodrinva_hyb_l_pdaf)(int*,
                                                                           int*,
                                                                           int*,
                                                                           double*,
                                                                           double*,
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
                                           void (*c__g2l_obs_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*,
                                                                   int*
                                                                  ),
                                           void (*c__init_obsvar_pdaf)(int*,
                                                                       int*,
                                                                       double*,
                                                                       double*
                                                                      ),
                                           void (*c__init_obsvar_l_pdaf)(int*,
                                                                         int*,
                                                                         int*,
                                                                         double*,
                                                                         int*,
                                                                         double*
                                                                        ),
                                           void (*c__likelihood_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        double*,
                                                                        double*
                                                                       ),
                                           void (*c__likelihood_hyb_l_pdaf)(int*,
                                                                            int*,
                                                                            int*,
                                                                            int*,
                                                                            double*,
                                                                            double*,
                                                                            double*,
                                                                            double*
                                                                           ),
                                           int* outflag
                                          ) noexcept;
cdef extern void c__pdaf_put_state_lseik (void (*c__collect_state_pdaf)(int*,
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
                                          void (*c__init_obs_pdaf)(int*,
                                                                   int*,
                                                                   double*
                                                                  ),
                                          void (*c__init_obs_l_pdaf)(int*,
                                                                     int*,
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
                                          void (*c__prodrinva_l_pdaf)(int*,
                                                                      int*,
                                                                      int*,
                                                                      int*,
                                                                      double*,
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
                                          void (*c__g2l_obs_pdaf)(int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*,
                                                                  int*
                                                                 ),
                                          void (*c__init_obsvar_pdaf)(int*,
                                                                      int*,
                                                                      double*,
                                                                      double*
                                                                     ),
                                          void (*c__init_obsvar_l_pdaf)(int*,
                                                                        int*,
                                                                        int*,
                                                                        double*,
                                                                        int*,
                                                                        double*
                                                                       ),
                                          int* flag
                                         ) noexcept;
cdef extern void c__pdaf_put_state_netf (void (*c__collect_state_pdaf)(int*,
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
                                         void (*c__init_obs_pdaf)(int*,
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
                                         void (*c__likelihood_pdaf)(int*,
                                                                    int*,
                                                                    double*,
                                                                    double*,
                                                                    double*
                                                                   ),
                                         int* flag
                                        ) noexcept;
cdef extern void c__pdaf_put_state_pf (void (*c__collect_state_pdaf)(int*,
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
                                       void (*c__init_obs_pdaf)(int*,
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
                                       void (*c__likelihood_pdaf)(int*,
                                                                  int*,
                                                                  double*,
                                                                  double*,
                                                                  double*
                                                                 ),
                                       int* flag
                                      ) noexcept;
cdef extern void c__pdaf_put_state_prepost (void (*c__collect_state_pdaf)(int*,
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
                                           ) noexcept;
cdef extern void c__pdaf_put_state_seek (void (*c__collect_state_pdaf)(int*,
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
                                         void (*c__init_obs_pdaf)(int*,
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
                                         void (*c__prodrinva_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   double*,
                                                                   double*,
                                                                   double*
                                                                  ),
                                         int* flag
                                        ) noexcept;
cdef extern void c__pdaf_put_state_seik (void (*c__collect_state_pdaf)(int*,
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
                                         void (*c__init_obs_pdaf)(int*,
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
                                         void (*c__prodrinva_pdaf)(int*,
                                                                   int*,
                                                                   int*,
                                                                   double*,
                                                                   double*,
                                                                   double*
                                                                  ),
                                         void (*c__init_obsvar_pdaf)(int*,
                                                                     int*,
                                                                     double*,
                                                                     double*
                                                                    ),
                                         int* flag
                                        ) noexcept;
cdef extern void c__pdaf_reset_forget (double* forget_in
                                      ) noexcept;
cdef extern void c__pdaf_sampleens (int* dim,
                                    int* dim_ens,
                                    double* modes,
                                    double* svals,
                                    double* state,
                                    double* ens,
                                    int* verbose,
                                    int* flag
                                   ) noexcept;
cdef extern void c__pdaf_set_debug_flag (int* debugval
                                        ) noexcept;
cdef extern void c__pdaf_set_ens_pointer (double** c_ens_point,
                                          int* dims,
                                          int* status
                                         ) noexcept;
cdef extern void c__pdaf_set_smootherens (double** c_sens_point,
                                          int* maxlag,
                                          int* dims,
                                          int* status
                                         ) noexcept;
cdef extern void c__pdaf_seik_ttimesa (int* rank,
                                       int* dim_col,
                                       double* a,
                                       double* b
                                      ) noexcept;
cdef extern void c__pdaf_etkf_tleft (int* dim_ens,
                                     int* dim,
                                     double* a
                                    ) noexcept;
cdef extern void c__pdaf_estkf_omegaa (int* rank,
                                       int* dim_col,
                                       double* a,
                                       double* b
                                      ) noexcept;
cdef extern void c__pdaf_enkf_omega (int* seed,
                                     int* r,
                                     int* dim_ens,
                                     double* omega,
                                     double* norm,
                                     int* otype,
                                     int* screen
                                    ) noexcept;
cdef extern void c__pdaf_seik_omega (int* rank,
                                     double* omega,
                                     int* omegatype,
                                     int* screen
                                    ) noexcept;
cdef extern void c__pdaf_incremental (int* steps,
                                      void (*c__dist_stateinc_pdaf)(int*,
                                                                    double*,
                                                                    int*,
                                                                    int*
                                                                   )
                                     ) noexcept;
cdef extern void c__pdaf_add_increment (int* dim_p,
                                        double* state_p
                                       ) noexcept;
cdef extern void c__pdaf_local_weights (int* wtype,
                                        double* cradius,
                                        double* sradius,
                                        int* dim,
                                        double* distance,
                                        double* weight,
                                        int* verbose
                                       ) noexcept;
cdef extern void c__pdaf_diag_crps (int* dim,
                                    int* dim_ens,
                                    int* element,
                                    double* oens,
                                    double* obs,
                                    double* crps,
                                    double* reli,
                                    double* resol,
                                    double* uncert,
                                    int* status
                                   ) noexcept;
cdef extern void c__pdaf_force_analysis () noexcept;
cdef extern void c__pdaf_gather_obs_f2_flex (int* dim_obs_p,
                                             int* dim_obs_f,
                                             double* coords_p,
                                             double* coords_f,
                                             int* nrows,
                                             int* status
                                            ) noexcept;
cdef extern void c__pdaf_gather_obs_f_flex (int* dim_obs_p,
                                            int* dim_obs_f,
                                            double* obs_p,
                                            double* obs_f,
                                            int* status
                                           ) noexcept;
cdef extern void c__pdaf_prepost (void (*c__collect_state_pdaf)(int*,
                                                                double*
                                                               ),
                                  void (*c__distribute_state_pdaf)(int*,
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
                                 ) noexcept;
cdef extern void c__pdaf_set_memberid (int* memberid
                                      ) noexcept;
cdef extern void c__pdaf_set_comm_pdaf (int* in_comm_pdaf
                                       ) noexcept;
cdef extern void c__pdaf_set_offline_mode (int* screen
                                          ) noexcept;
cdef extern void c__pdaf_print_domain_stats (int* n_domains_p
                                            ) noexcept;
cdef extern void c__pdaf_init_local_obsstats () noexcept;
cdef extern void c__pdaf_incr_local_obsstats (int* dim_obs_l
                                             ) noexcept;
cdef extern void c__pdaf_print_local_obsstats (int* screen,
                                               int* n_domains_with_obs
                                              ) noexcept;
cdef extern void c__pdaf_omit_obs_omi (int* dim_p,
                                       int* dim_obs_p,
                                       int* dim_ens,
                                       double* state_p,
                                       double* ens_p,
                                       double* obs_p,
                                       void (*c__init_obs_pdaf)(int*,
                                                                int*,
                                                                double*
                                                               ),
                                       void (*c__obs_op_pdaf)(int*,
                                                              int*,
                                                              int*,
                                                              double*,
                                                              double*
                                                             ),
                                       int* compute_mean,
                                       int* screen
                                      ) noexcept;
cdef extern void c__pdafomi_init (int* n_obs
                                 ) noexcept;
cdef extern void c__pdafomi_set_doassim (int* i_obs,
                                         int* doassim
                                        ) noexcept;
cdef extern void c__pdafomi_set_disttype (int* i_obs,
                                          int* disttype
                                         ) noexcept;
cdef extern void c__pdafomi_set_ncoord (int* i_obs,
                                        int* ncoord
                                       ) noexcept;
cdef extern void c__pdafomi_set_id_obs_p (int* i_obs,
                                          int* nrows,
                                          int* dim_obs_p,
                                          int* id_obs_p
                                         ) noexcept;
cdef extern void c__pdafomi_set_icoeff_p (int* i_obs,
                                          int* nrows,
                                          int* dim_obs_p,
                                          double* icoeff_p
                                         ) noexcept;
cdef extern void c__pdafomi_set_domainsize (int* i_obs,
                                            int* ncoord,
                                            double* domainsize
                                           ) noexcept;
cdef extern void c__pdafomi_set_obs_err_type (int* i_obs,
                                              int* obs_err_type
                                             ) noexcept;
cdef extern void c__pdafomi_set_use_global_obs (int* i_obs,
                                                int* use_global_obs
                                               ) noexcept;
cdef extern void c__pdafomi_set_inno_omit (int* i_obs,
                                           double* inno_omit
                                          ) noexcept;
cdef extern void c__pdafomi_set_inno_omit_ivar (int* i_obs,
                                                double* inno_omit_ivar
                                               ) noexcept;
cdef extern void c__pdafomi_gather_obs (int* i_obs,
                                        int* dim_obs_p,
                                        double* obs_p,
                                        double* ivar_obs_p,
                                        double* ocoord_p,
                                        double* cradius,
                                        int* dim_obs
                                       ) noexcept;
cdef extern void c__pdafomi_gather_obsstate (int* i_obs,
                                             double* obsstate_p,
                                             double* obsstate_f,
                                             int* nobs_f_all
                                            ) noexcept;
cdef extern void c__pdafomi_set_domain_limits (double* lim_coords
                                              ) noexcept;
cdef extern void c__pdafomi_set_debug_flag (int* debugval
                                           ) noexcept;
cdef extern void c__pdafomi_deallocate_obs (int* i_obs
                                           ) noexcept;
cdef extern void c__pdafomi_obs_op_gridpoint (int* i_obs,
                                              double* state_p,
                                              int* dim_p,
                                              double* obs_f_all,
                                              int* nobs_f_all
                                             ) noexcept;
cdef extern void c__pdafomi_obs_op_gridavg (int* i_obs,
                                            int* nrows,
                                            double* state_p,
                                            int* dim_p,
                                            double* obs_f_all,
                                            int* nobs_f_all
                                           ) noexcept;
cdef extern void c__pdafomi_obs_op_interp_lin (int* i_obs,
                                               int* nrows,
                                               double* state_p,
                                               int* dim_p,
                                               double* obs_f_all,
                                               int* nobs_f_all
                                              ) noexcept;
cdef extern void c__pdafomi_obs_op_adj_gridavg (int* i_obs,
                                                int* nrows,
                                                double* state_p,
                                                int* dim_p,
                                                double* obs_f_all,
                                                int* nobs_f_all
                                               ) noexcept;
cdef extern void c__pdafomi_obs_op_adj_gridpoint (int* i_obs,
                                                  double* state_p,
                                                  int* dim_p,
                                                  double* obs_f_all,
                                                  int* nobs_f_all
                                                 ) noexcept;
cdef extern void c__pdafomi_obs_op_adj_interp_lin (int* i_obs,
                                                   int* nrows,
                                                   double* state_p,
                                                   int* dim_p,
                                                   double* obs_f_all,
                                                   int* nobs_f_all
                                                  ) noexcept;
cdef extern void c__pdafomi_get_interp_coeff_tri (double* gpc,
                                                  double* oc,
                                                  double* icoeff
                                                 ) noexcept;
cdef extern void c__pdafomi_get_interp_coeff_lin1d (double* gpc,
                                                    double* oc,
                                                    double* icoeff
                                                   ) noexcept;
cdef extern void c__pdafomi_get_interp_coeff_lin (int* num_gp,
                                                  int* n_dim,
                                                  double* gpc,
                                                  double* oc,
                                                  double* icoeff
                                                 ) noexcept;
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
                                             ) noexcept;
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
                                                     ) noexcept;
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
                                                      ) noexcept;
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
                                              ) noexcept;
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
                                                      ) noexcept;
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
                                                       ) noexcept;
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
                                             ) noexcept;
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
                                             ) noexcept;
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
                                         ) noexcept;
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
                                            ) noexcept;
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
                                                    ) noexcept;
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
                                                     ) noexcept;
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
                                                   ) noexcept;
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
                                             ) noexcept;
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
                                                     ) noexcept;
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
                                                      ) noexcept;
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
                                            ) noexcept;
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
                                            ) noexcept;
cdef extern void c__pdafomi_init_obs_f_cb (int* step,
                                           int* dim_obs_f,
                                           double* observation_f
                                          ) noexcept;
cdef extern void c__pdafomi_init_obsvar_cb (int* step,
                                            int* dim_obs_p,
                                            double* obs_p,
                                            double* meanvar
                                           ) noexcept;
cdef extern void c__pdafomi_g2l_obs_cb (int* domain_p,
                                        int* step,
                                        int* dim_obs_f,
                                        int* dim_obs_l,
                                        double* ostate_f,
                                        double* ostate_l
                                       ) noexcept;
cdef extern void c__pdafomi_init_obs_l_cb (int* domain_p,
                                           int* step,
                                           int* dim_obs_l,
                                           double* observation_l
                                          ) noexcept;
cdef extern void c__pdafomi_init_obsvar_l_cb (int* domain_p,
                                              int* step,
                                              int* dim_obs_l,
                                              double* obs_l,
                                              double* meanvar_l
                                             ) noexcept;
cdef extern void c__pdafomi_prodrinva_l_cb (int* domain_p,
                                            int* step,
                                            int* dim_obs_l,
                                            int* rank,
                                            double* obs_l,
                                            double* a_l,
                                            double* c_l
                                           ) noexcept;
cdef extern void c__pdafomi_likelihood_l_cb (int* domain_p,
                                             int* step,
                                             int* dim_obs_l,
                                             double* obs_l,
                                             double* resid_l,
                                             double* lhood_l
                                            ) noexcept;
cdef extern void c__pdafomi_prodrinva_cb (int* step,
                                          int* dim_obs_p,
                                          int* ncol,
                                          double* obs_p,
                                          double* a_p,
                                          double* c_p
                                         ) noexcept;
cdef extern void c__pdafomi_likelihood_cb (int* step,
                                           int* dim_obs,
                                           double* obs,
                                           double* resid,
                                           double* lhood
                                          ) noexcept;
cdef extern void c__pdafomi_add_obs_error_cb (int* step,
                                              int* dim_obs_p,
                                              double* c_p
                                             ) noexcept;
cdef extern void c__pdafomi_init_obscovar_cb (int* step,
                                              int* dim_obs,
                                              int* dim_obs_p,
                                              double* covar,
                                              double* m_state_p,
                                              bint* isdiag
                                             ) noexcept;
cdef extern void c__pdafomi_init_obserr_f_cb (int* step,
                                              int* dim_obs_f,
                                              double* obs_f,
                                              double* obserr_f
                                             ) noexcept;
cdef extern void c__pdafomi_prodrinva_hyb_l_cb (int* domain_p,
                                                int* step,
                                                int* dim_obs_l,
                                                int* rank,
                                                double* obs_l,
                                                double* alpha,
                                                double* a_l,
                                                double* c_l
                                               ) noexcept;
cdef extern void c__pdafomi_likelihood_hyb_l_cb (int* domain_p,
                                                 int* step,
                                                 int* dim_obs_l,
                                                 double* obs_l,
                                                 double* resid_l,
                                                 double* alpha,
                                                 double* lhood_l
                                                ) noexcept;
cdef extern void c__pdafomi_obsstats_l (int* screen
                                       ) noexcept;
cdef extern void c__pdafomi_weights_l (int* verbose,
                                       int* nobs_l,
                                       int* ncols,
                                       int* locweight,
                                       double* cradius,
                                       double* sradius,
                                       double* mata,
                                       double* ivar_obs_l,
                                       double* dist_l,
                                       double* weight_l
                                      ) noexcept;
cdef extern void c__pdafomi_weights_l_sgnl (int* verbose,
                                            int* nobs_l,
                                            int* ncols,
                                            int* locweight,
                                            double* cradius,
                                            double* sradius,
                                            double* mata,
                                            double* ivar_obs_l,
                                            double* dist_l,
                                            double* weight_l
                                           ) noexcept;
cdef extern void c__pdafomi_check_error (int* flag
                                        ) noexcept;
cdef extern void c__pdafomi_gather_obsdims () noexcept;
cdef extern void c__pdafomi_obsstats (int* screen
                                     ) noexcept;
cdef extern void c__pdafomi_init_dim_obs_l_iso (int* i_obs,
                                                int* ncoord,
                                                double* coords_l,
                                                int* locweight,
                                                double* cradius,
                                                double* sradius,
                                                int* cnt_obs_l
                                               ) noexcept;
cdef extern void c__pdafomi_init_dim_obs_l_noniso (int* i_obs,
                                                   int* ncoord,
                                                   double* coords_l,
                                                   int* locweight,
                                                   double* cradius,
                                                   double* sradius,
                                                   int* cnt_obs_l
                                                  ) noexcept;
cdef extern void c__pdafomi_init_dim_obs_l_noniso_locweights (int* i_obs,
                                                              int* ncoord,
                                                              double* coords_l,
                                                              int* locweights,
                                                              double* cradius,
                                                              double* sradius,
                                                              int* cnt_obs_l
                                                             ) noexcept;
cdef extern void c__pdafomi_localize_covar_iso (int* i_obs,
                                                int* dim_p,
                                                int* dim_obs,
                                                int* ncoord,
                                                int* locweight,
                                                double* cradius,
                                                double* sradius,
                                                double* coords,
                                                double* hp,
                                                double* hph
                                               ) noexcept;
cdef extern void c__pdafomi_localize_covar_noniso (int* i_obs,
                                                   int* dim_p,
                                                   int* dim_obs,
                                                   int* ncoord,
                                                   int* locweight,
                                                   double* cradius,
                                                   double* sradius,
                                                   double* coords,
                                                   double* hp,
                                                   double* hph
                                                  ) noexcept;
cdef extern void c__pdafomi_localize_covar_noniso_locweights (int* i_obs,
                                                              int* dim_p,
                                                              int* dim_obs,
                                                              int* ncoord,
                                                              int* locweights,
                                                              double* cradius,
                                                              double* sradius,
                                                              double* coords,
                                                              double* hp,
                                                              double* hph
                                                             ) noexcept;
cdef extern void c__pdafomi_omit_by_inno_l_cb (int* domain_p,
                                               int* dim_obs_l,
                                               double* resid_l,
                                               double* obs_l
                                              ) noexcept;
cdef extern void c__pdafomi_omit_by_inno_cb (int* dim_obs_f,
                                             double* resid_f,
                                             double* obs_f
                                            ) noexcept;
