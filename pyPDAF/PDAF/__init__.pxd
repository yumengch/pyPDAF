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
                                          );
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
                                                  );
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
                                                   );
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
                                         );
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
                                          );
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
                                         );
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
                                                   );
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
                                                    );
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
                                          );
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
                                           );
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
                                          );
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
                                          );
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
                                          );
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
                                         );
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
                                       );
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
                                         );
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
                                         );
cdef extern void c__pdaf_deallocate ();
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
                                   );
cdef extern void c__pdaf_diag_effsample (int* dim_sample,
                                         double* weights,
                                         double* effsample
                                        );
cdef extern void c__pdaf_diag_ensstats (int* dim,
                                        int* dim_ens,
                                        int* element,
                                        double* state,
                                        double* ens,
                                        double* skewness,
                                        double* kurtosis,
                                        int* status
                                       );
cdef extern void c__pdaf_diag_histogram (int* ncall,
                                         int* dim,
                                         int* dim_ens,
                                         int* element,
                                         double* state,
                                         double* ens,
                                         int* hist,
                                         double* delta,
                                         int* status
                                        );
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
                                  );
cdef extern void c__pdaf_force_analysis ();
cdef extern void c__pdaf_gather_dim_obs_f (int* dim_obs_p,
                                           int* dim_obs_f
                                          );
cdef extern void c__pdaf_gather_obs_f (double* obs_p,
                                       int* dimobs_p,
                                       double* obs_f,
                                       int* dimobs_f,
                                       int* status
                                      );
cdef extern void c__pdaf_gather_obs_f2 (double* coords_p,
                                        int* dimobs_p,
                                        double* coords_f,
                                        int* dimobs_f,
                                        int* nrows,
                                        int* status
                                       );
cdef extern void c__pdaf_gather_obs_f2_flex (int* dim_obs_p,
                                             int* dim_obs_f,
                                             double* coords_p,
                                             double* coords_f,
                                             int* nrows,
                                             int* status
                                            );
cdef extern void c__pdaf_gather_obs_f_flex (int* dim_obs_p,
                                            int* dim_obs_f,
                                            double* obs_p,
                                            double* obs_f,
                                            int* status
                                           );
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
                                      );
cdef extern void c__pdaf_get_assim_flag (int* did_assim
                                        );
cdef extern void c__pdaf_get_localfilter (int* lfilter
                                         );
cdef extern void c__pdaf_get_memberid (int* memberid
                                      );
cdef extern void c__pdaf_get_obsmemberid (int* memberid
                                         );
cdef extern void c__pdaf_get_smootherens (double** c_sens_point,
                                          int* maxlag,
                                          int* dims,
                                          int* status
                                         );
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
                                   );
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
                              );
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
                                      );
cdef extern void c__pdaf_local_weights (int* wtype,
                                        double* cradius,
                                        double* sradius,
                                        int* dim,
                                        double* distance,
                                        double* weight,
                                        int* verbose
                                       );
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
                                 );
cdef extern void c__pdaf_print_info (int* printtype
                                    );
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
                                         );
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
                                                 );
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
                                                  );
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
                                        );
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
                                         );
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
                                        );
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
                                                );
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
                                                  );
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
                                                   );
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
                                         );
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
                                          );
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
                                         );
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
                                         );
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
                                         );
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
                                        );
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
                                      );
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
                                           );
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
                                        );
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
                                        );
cdef extern void c__pdaf_reset_forget (double* forget_in
                                      );
cdef extern void c__pdaf_sampleens (int* dim,
                                    int* dim_ens,
                                    double* modes,
                                    double* svals,
                                    double* state,
                                    double* ens,
                                    int* verbose,
                                    int* flag
                                   );
cdef extern void c__pdaf_set_ens_pointer (double** c_ens_point,
                                          int* dims,
                                          int* status
                                         );
cdef extern void c__pdaf_set_memberid (int* memberid
                                      );
cdef extern void c__pdaf_set_smootherens (double** c_sens_point,
                                          int* maxlag,
                                          int* dims,
                                          int* status
                                         );
cdef extern void c__pdaf_seik_ttimesa (int* rank,
                                       int* dim_col,
                                       double* a,
                                       double* b
                                      );
cdef extern void c__pdaf_etkf_tleft (int* dim_ens,
                                     int* dim,
                                     double* a
                                    );
cdef extern void c__pdaf_estkf_omegaa (int* rank,
                                       int* dim_col,
                                       double* a,
                                       double* b
                                      );
cdef extern void c__pdaf_enkf_omega (int* seed,
                                     int* r,
                                     int* dim_ens,
                                     double* omega,
                                     double* norm,
                                     int* otype,
                                     int* screen
                                    );
cdef extern void c__pdaf_seik_omega (int* rank,
                                     double* omega,
                                     int* omegatype,
                                     int* screen
                                    );
cdef extern void c__pdaf_incremental (int* steps,
                                      void (*c__dist_stateinc_pdaf)(int*,
                                                                    double*,
                                                                    int*,
                                                                    int*
                                                                   )
                                     );
cdef extern void c__pdaf_add_increment (int* dim_p,
                                        double* state_p
                                       );
