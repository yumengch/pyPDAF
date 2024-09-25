docstrings = dict()
docstrings['assimilate_estkf'] = "Using ESTKF (error space transform " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This function should be called at each model time step. " \
                                 "The ESTKF is a more efficient equivalent to the ETKF. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_estkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for ensemble mean)" \
                                 "5. py__init_obs_pdaf " \
                                 "6. py__obs_op_pdaf (for each ensemble member)" \
                                 "7. py__init_obsvar_pdaf " \
                                 "(only relevant for adaptive forgetting factor schemes)" \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. core DA algorithm " \
                                 "10. py__prepoststep_state_pdaf " \
                                 "11. py__distribute_state_pdaf " \
                                 "12. py__next_observation_pdaf "
docstrings['assimilate_lestkf'] = "Using Local ESTKF (error space transform " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step. " \
                                 "The LESTKF is a more efficient equivalent to the LETKF. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lestkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "14. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "15. py__prodRinvA_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__prepoststep_state_pdaf " \
                                 "19. py__distribute_state_pdaf " \
                                 "20. py__next_observation_pdaf "
docstrings['assimilate_etkf'] = "Using ETKF (ensemble transform " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_etkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for ensemble mean)" \
                                 "5. py__init_obs_pdaf " \
                                 "6. py__obs_op_pdaf (for each ensemble member)" \
                                 "7. py__init_obsvar_pdaf " \
                                 "(only relevant for adaptive forgetting factor schemes)" \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. core DA algorithm " \
                                 "10. py__prepoststep_state_pdaf " \
                                 "11. py__distribute_state_pdaf " \
                                 "12. py__next_observation_pdaf "
docstrings['assimilate_letkf'] = "Using local ensemble transform Kalman filter for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_letkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "14. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "15. py__prodRinvA_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__prepoststep_state_pdaf " \
                                 "19. py__distribute_state_pdaf " \
                                 "20. py__next_observation_pdaf "
docstrings['assimilate_enkf'] = "Using stochastic EnKF (ensemble " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_enkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for ensemble mean)" \
                                 "5. py__add_obs_err_pdaf " \
                                 "6. py__init_obs_pdaf " \
                                 "7. py__init_obscovar_pdaf " \
                                 "8. py__obs_op_pdaf (for each ensemble member)" \
                                 "9. core DA algorithm " \
                                 "10. py__prepoststep_state_pdaf " \
                                 "11. py__distribute_state_pdaf " \
                                 "12. py__next_observation_pdaf "
docstrings['assimilate_lenkf'] = "Using stochastic EnKF (ensemble " \
                                 "Kalman filter) with covariance localisation " \
                                 "for DA without OMI. " \
                                 "This is the only scheme for covariance localisation in PDAF. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lenkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for each ensemble member)" \
                                 "5. py__localize_pdaf " \
                                 "6. py__add_obs_err_pdaf " \
                                 "7. py__init_obs_pdaf " \
                                 "8. py__init_obscovar_pdaf " \
                                 "9. py__obs_op_pdaf (repeated to reduce storage)" \
                                 "10. core DA algorithm" \
                                 "11. py__prepoststep_state_pdaf " \
                                 "12. py__distribute_state_pdaf " \
                                 "13. py__next_observation_pdaf "
docstrings['assimilate_3dvar'] = "Using 3DVar for DA without OMI. " \
								 "This is a deterministic filtering scheme. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_3dvar` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:" \
                                 "6. py__cvt_pdaf " \
                                 "7. py__obs_op_lin_pdaf " \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. py__obs_op_adj_pdaf " \
                                 "10. py__cvt_adj_pdaf " \
                                 "11. core DA algorithm " \
                                 "After the iterations: " \
                                 "12. py__cvt_pdaf " \
                                 "13. py__prepoststep_state_pdaf " \
                                 "14. py__distribute_state_pdaf " \
                                 "15. py__next_observation_pdaf "
docstrings['assimilate_en3dvar_estkf'] = "Using 3DEnVar for DA without OMI. " \
								 "The background error covariance matrix is estimated by ensemble. " \
								 "The 3DEnVar only calculates the analysis of the ensemble mean. " \
								 "An ESTKF is used to generate ensemble perturbations. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_en3dvar_estkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:" \
                                 "6. py__cvt_ens_pdaf " \
                                 "7. py__obs_op_lin_pdaf " \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. py__obs_op_adj_pdaf " \
                                 "10. py__cvt_adj_ens_pdaf " \
                                 "11. core 3DEnVar algorithm " \
                                 "After the iterations: " \
                                 "12. py__cvt_ens_pdaf " \
                                 "Perform ESTKF: " \
                                 "13. py__init_dim_obs_pdaf " \
                                 "14. py__obs_op_pdaf (for ensemble mean)" \
                                 "15. py__init_obs_pdaf " \
                                 "16. py__obs_op_pdaf (for each ensemble member)" \
                                 "17. py__init_obsvar_pdaf " \
                                 "(only relevant for adaptive forgetting factor schemes)" \
                                 "18. py__prodRinvA_pdaf " \
                                 "19. core ESTKF algorithm " \
                                 "20. py__prepoststep_state_pdaf " \
                                 "21. py__distribute_state_pdaf " \
                                 "22. py__next_observation_pdaf "
docstrings['assimilate_en3dvar_lestkf'] = "Using 3DEnVar for DA without OMI. " \
								 "The background error covariance matrix is estimated by ensemble. " \
								 "The 3DEnVar only calculates the analysis of the ensemble mean. " \
								 "An LESTKF is used to generate ensemble perturbations. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_en3dvar_lestkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:" \
                                 "6. py__cvt_ens_pdaf " \
                                 "7. py__obs_op_lin_pdaf " \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. py__obs_op_adj_pdaf " \
                                 "10. py__cvt_adj_ens_pdaf " \
                                 "11. core DA algorithm " \
                                 "After the iterations: " \
                                 "12. py__cvt_ens_pdaf " \
                                 "Perform LESTKF: " \
                                 "13. py__init_n_domains_p_pdaf " \
                                 "14. py__init_dim_obs_pdaf " \
                                 "15. py__obs_op_pdaf (for each ensemble member)" \
                                 "16. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "17. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "18. py__init_dim_l_pdaf " \
                                 "19. py__init_dim_obs_l_pdaf " \
                                 "20. py__g2l_state_pdaf " \
                                 "21. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "22. py__init_obs_l_pdaf "\
                                 "23. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "24. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "25. py__prodRinvA_pdaf " \
                                 "26. core DA algorithm " \
                                 "27. py__l2g_state_pdaf " \
                                 "28. py__prepoststep_state_pdaf " \
                                 "29. py__distribute_state_pdaf " \
                                 "30. py__next_observation_pdaf "
docstrings['assimilate_hyb3dvar_estkf'] = "Using Hybrid 3DEnVar for DA without OMI. " \
										  "Here, the background error covariance is hybridised by a static background error covariance, " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "ESTKF in this implementation. " \
                                          "This function should be called at each model time step. \n\n" \
                                          "The function is a combination of `pyPDAF.PDAF.put_state_hyb3dvar_estkf` " \
                                          "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                          "in the following sequence: \n" \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:" \
			                              "6. py__cvt_pdaf " \
			                              "7. py__cvt_ens_pdaf " \
			                              "8. py__obs_op_lin_pdaf " \
			                              "9. py__prodRinvA_pdaf " \
			                              "10. py__obs_op_adj_pdaf " \
			                              "11. py__cvt_adj_pdaf " \
			                              "12. py__cvt_adj_ens_pdaf " \
			                              "13. core 3DEnVar algorithm " \
			                              "After the iterations: " \
			                              "14. py__cvt_pdaf " \
			                              "15. py__cvt_ens_pdaf " \
			                              "Perform ESTKF: " \
			                              "16. py__init_dim_obs_pdaf " \
			                              "17. py__obs_op_pdaf (for ensemble mean)" \
			                              "18. py__init_obs_pdaf " \
			                              "19. py__obs_op_pdaf (for each ensemble member)" \
			                              "20. py__init_obsvar_pdaf " \
			                              "(only relevant for adaptive forgetting factor schemes)" \
			                              "21. py__prodRinvA_pdaf " \
			                              "22. core ESTKF algorithm " \
			                              "23. py__prepoststep_state_pdaf " \
			                              "24. py__distribute_state_pdaf " \
			                              "25. py__next_observation_pdaf "
docstrings['assimilate_hyb3dvar_lestkf'] = "Using Hybrid 3DEnVar for DA without OMI. " \
										  "Here, the background error covariance is hybridised by a static background error covariance, " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "LESTKF in this implementation. " \
                                          "This function should be called at each model time step. \n\n" \
                                          "The function is a combination of `pyPDAF.PDAF.put_state_hyb3dvar_lestkf` " \
                                          "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                          "in the following sequence: \n" \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:" \
			                              "6. py__cvt_pdaf " \
			                              "7. py__cvt_ens_pdaf " \
			                              "8. py__obs_op_lin_pdaf " \
			                              "9. py__prodRinvA_pdaf " \
			                              "10. py__obs_op_adj_pdaf " \
			                              "11. py__cvt_adj_pdaf " \
			                              "12. py__cvt_adj_ens_pdaf " \
			                              "13. core DA algorithm " \
			                              "After the iterations: " \
			                              "14. py__cvt_pdaf " \
			                              "15. py__cvt_ens_pdaf " \
			                              "Perform LESTKF: " \
			                              "16. py__init_n_domains_p_pdaf " \
			                              "17. py__init_dim_obs_pdaf " \
			                              "18. py__obs_op_pdaf (for each ensemble member)" \
			                              "19. py__init_obs_pdaf " \
			                              "(if global adaptive forgetting factor is used " \
			                              "(type_forget=1 in pyPDAF.PDAF.init)) " \
			                              "20. py__init_obsvar_pdaf " \
			                              "(if global adaptive forgetting factor is used) " \
			                              "loop over each local domain:" \
			                              "21. py__init_dim_l_pdaf " \
			                              "22. py__init_dim_obs_l_pdaf " \
			                              "23. py__g2l_state_pdaf " \
			                              "24. py__g2l_obs_pdaf " \
			                              "(localise mean ensemble in observation space)" \
			                              "25. py__init_obs_l_pdaf "\
			                              "26. py__g2l_obs_pdaf " \
			                              "(localise each ensemble member in observation space) " \
			                              "27. py__init_obsvar_l_pdaf " \
			                              "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
			                              "28. py__prodRinvA_pdaf " \
			                              "29. core DA algorithm " \
			                              "30. py__l2g_state_pdaf " \
			                              "31. py__prepoststep_state_pdaf " \
			                              "32. py__distribute_state_pdaf " \
			                              "33. py__next_observation_pdaf "
docstrings['assimilate_lnetf'] = "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lnetf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "loop over each local domain:" \
                                 "6. py__init_dim_l_pdaf " \
                                 "7. py__init_dim_obs_l_pdaf " \
                                 "8. py__g2l_state_pdaf " \
                                 "9. py__init_obs_l_pdaf "\
                                 "10. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)" \
                                 "11. py__likelihood_l_pdaf " \
                                 "12. core DA algorithm " \
                                 "13. py__l2g_state_pdaf " \
                                 "14. py__prepoststep_state_pdaf " \
                                 "15. py__distribute_state_pdaf " \
                                 "16. py__next_observation_pdaf "
docstrings['assimilate_lknetf'] = "This function will is a hybridised LETKF and LNETF " \
                                 "for DA without OMI. The LNETF computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "The hybridisation with LETKF is expected to lead to improved performance for " \
                                 "quasi-Gaussian problems. " \
                                 "The function should be called at each model step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lknetf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "14. py__prodRinvA_pdaf " \
                                 "15. py__likelihood_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member)" \
                                 "19. py__likelihood_hyb_l_pdaf" \
                                 "20. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "21. py__prodRinvA_hyb_l_pdaf " \
                                 "22. py__prepoststep_state_pdaf " \
                                 "23. py__distribute_state_pdaf " \
                                 "24. py__next_observation_pdaf "
docstrings['assimilate_lseik'] = "Using local singular evolutive interpolated Kalman filter for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step.\n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lseik` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "14. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "15. py__prodRinvA_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__prepoststep_state_pdaf " \
                                 "19. py__distribute_state_pdaf " \
                                 "20. py__next_observation_pdaf "
docstrings['assimilate_netf'] = "This function will use Nonlinear Ensemble Transform Filter (NETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "The function should be called at each model step. \n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_netf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__init_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__likelihood_pdaf " \
                                 "7. core DA algorithm " \
                                 "8. py__prepoststep_state_pdaf " \
                                 "9. py__distribute_state_pdaf " \
                                 "10. py__next_observation_pdaf "
docstrings['assimilate_pf'] = "This function will use particle filter for DA without OMI. " \
                              "This is a fully nonlinear filter. " \
                              "The function should be called at each model step. \n\n" \
                              "The function is a combination of `pyPDAF.PDAF.put_state_pf` " \
                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                              "in the following sequence: \n" \
                              "1. py__collect_state_pdaf " \
                              "2. py__prepoststep_state_pdaf " \
                              "3. py__init_dim_obs_pdaf " \
                              "4. py__init_obs_pdaf " \
                              "5. py__obs_op_pdaf (for each ensemble member)" \
                              "6. py__likelihood_pdaf " \
                              "7. core DA algorithm " \
                              "8. py__prepoststep_state_pdaf " \
                              "9. py__distribute_state_pdaf " \
                              "10. py__next_observation_pdaf "
docstrings['assimilate_seek'] = "This function will use singular evolutive extended Kalman filter for DA without OMI. " \
                                 "This is a deterministic Kalman filter. " \
                                 "The function should be called at each model step.\n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_seek` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for ensemble mean)" \
                                 "5. py__init_obs_pdaf " \
                                 "6. py__obs_op_pdaf (for each ensemble member)" \
                                 "7. py__prodRinvA_pdaf " \
                                 "8. core DA algorithm " \
                                 "9. py__prepoststep_state_pdaf " \
                                 "10. py__distribute_state_pdaf " \
                                 "11. py__next_observation_pdaf "
docstrings['assimilate_seik'] = "This function will use singular evolutive interpolated Kalman filter for DA without OMI. " \
                                 "The function should be called at each model step.\n\n" \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_seik` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for ensemble mean)" \
                                 "5. py__init_obs_pdaf " \
                                 "6. py__obs_op_pdaf (for each ensemble member)" \
                                 "7. py__init_obsvar_pdaf " \
                                 "(only relevant for adaptive forgetting factor schemes)" \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. core DA algorithm " \
                                 "10. py__prepoststep_state_pdaf " \
                                 "11. py__distribute_state_pdaf " \
                                 "12. py__next_observation_pdaf "
docstrings['assimilate_prepost'] = "This function does not perform any DA. " \
                                   "It is used to perform a preprocess and postprocess of the ensemble. \n" \
                                   "The function is a combination of `pyPDAF.PDAF.put_state_prepost` " \
                                   "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                   "in the following sequence: \n" \
                                   "1. py__collect_state_pdaf " \
                                   "2. py__prepoststep_state_pdaf " \
                                   "3. py__prepoststep_state_pdaf " \
                                   "4. py__distribute_state_pdaf " \
                                   "5. py__next_observation_pdaf "
docstrings['generate_obs'] = "This function generates synthetic observations based on each member of model forecast. " \
                             "This is based on the usual implementation strategy for PDAF.\n\n" \
                             "The function is a combination of `pyPDAF.PDAF.put_state_generate_obs` " \
                             "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                             "in the following sequence: \n" \
                             "1. py__collect_state_pdaf " \
                             "2. py__prepoststep_state_pdaf " \
                             "3. py__init_dim_obs_pdaf " \
                             "4. py__obs_op_pdaf" \
                             "7. py__init_obserr_f_pdaf " \
                             "5. py__get_obs_f_pdaf " \
                             "6. py__prepoststep_state_pdaf " \
                             "7. py__distribute_state_pdaf " \
                             "8. py__next_observation_pdaf "

docstrings['put_state_3dvar'] = "Using 3DVar for DA without OMI. " \
								"This is a deterministic filtering scheme. " \
                                "This function is usually used in 'flexible' parallelisation, " \
                                "but 3dvar is deterministic and does not require ensemble. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process the state vector and " \
                                "distribute the state vector " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n\n" \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n" \
                                "1. py__collect_state_pdaf " \
                                "2. py__prepoststep_state_pdaf " \
                                "3. py__init_dim_obs_pdaf " \
                                "4. py__obs_op_pdaf " \
                                "5. py__init_obs_pdaf " \
                                "Starting the iterative optimisation:" \
                                "6. py__cvt_pdaf " \
                                "7. py__obs_op_lin_pdaf " \
                                "8. py__prodRinvA_pdaf " \
                                "9. py__obs_op_adj_pdaf " \
                                "10. py__cvt_adj_pdaf " \
                                "11. core DA algorithm " \
                                "After the iterations: " \
                                "12. py__cvt_pdaf "
docstrings['put_state_en3dvar_estkf'] = "Using 3DEnVar for DA without OMI. " \
								 "The background error covariance matrix is estimated by ensemble. " \
								 "The 3DEnVar only calculates the analysis of the ensemble mean. " \
								 "An ESTKF is used to generate ensemble perturbations. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:" \
                                 "6. py__cvt_ens_pdaf " \
                                 "7. py__obs_op_lin_pdaf " \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. py__obs_op_adj_pdaf " \
                                 "10. py__cvt_adj_ens_pdaf " \
                                 "11. core 3DEnVar algorithm " \
                                 "After the iterations: " \
                                 "12. py__cvt_ens_pdaf " \
                                 "Perform ESTKF: " \
                                 "13. py__init_dim_obs_pdaf " \
                                 "14. py__obs_op_pdaf (for ensemble mean)" \
                                 "15. py__init_obs_pdaf " \
                                 "16. py__obs_op_pdaf (for each ensemble member)" \
                                 "17. py__init_obsvar_pdaf " \
                                 "(only relevant for adaptive forgetting factor schemes)" \
                                 "18. py__prodRinvA_pdaf " \
                                 "19. core ESTKF algorithm "
docstrings['put_state_en3dvar_lestkf'] = "Using 3DEnVar for DA without OMI. " \
								 "The background error covariance matrix is estimated by ensemble. " \
								 "The 3DEnVar only calculates the analysis of the ensemble mean. " \
								 "An LESTKF is used to generate ensemble perturbations. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:" \
                                 "6. py__cvt_ens_pdaf " \
                                 "7. py__obs_op_lin_pdaf " \
                                 "8. py__prodRinvA_pdaf " \
                                 "9. py__obs_op_adj_pdaf " \
                                 "10. py__cvt_adj_ens_pdaf " \
                                 "11. core DA algorithm " \
                                 "After the iterations: " \
                                 "12. py__cvt_ens_pdaf " \
                                 "Perform LESTKF: " \
                                 "13. py__init_n_domains_p_pdaf " \
                                 "14. py__init_dim_obs_pdaf " \
                                 "15. py__obs_op_pdaf (for each ensemble member)" \
                                 "16. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "17. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "18. py__init_dim_l_pdaf " \
                                 "19. py__init_dim_obs_l_pdaf " \
                                 "20. py__g2l_state_pdaf " \
                                 "21. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "22. py__init_obs_l_pdaf "\
                                 "23. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "24. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "25. py__prodRinvA_pdaf " \
                                 "26. core DA algorithm " \
                                 "27. py__l2g_state_pdaf "
docstrings['put_state_enkf'] = "Using stochastic EnKF (ensemble " \
                               "Kalman filter) for DA without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n\n" \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n" \
                               "1. py__collect_state_pdaf " \
                               "2. py__prepoststep_state_pdaf " \
                               "3. py__init_dim_obs_pdaf " \
                               "4. py__obs_op_pdaf (for ensemble mean)" \
                               "5. py__add_obs_err_pdaf " \
                               "6. py__init_obs_pdaf " \
                               "7. py__init_obscovar_pdaf " \
                               "8. py__obs_op_pdaf (for each ensemble member)" \
                               "9. core DA algorithm "
docstrings['put_state_estkf'] = "Using ESTKF (error space transform " \
                                "Kalman filter) for DA without OMI. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. " \
                                "The ESTKF is a more efficient equivalent to the ETKF. \n\n" \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n" \
                                "1. py__collect_state_pdaf " \
                                "2. py__prepoststep_state_pdaf " \
                                "3. py__init_dim_obs_pdaf " \
                                "4. py__obs_op_pdaf (for ensemble mean)" \
                                "5. py__init_obs_pdaf " \
                                "6. py__obs_op_pdaf (for each ensemble member)" \
                                "7. py__init_obsvar_pdaf " \
                                "(only relevant for adaptive forgetting factor schemes)" \
                                "8. py__prodRinvA_pdaf " \
                                "9. core DA algorithm "
docstrings['put_state_etkf'] = "Using ETKF (ensemble transform " \
                               "Kalman filter) for DA without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n\n" \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n" \
                               "1. py__collect_state_pdaf " \
                               "2. py__prepoststep_state_pdaf " \
                               "3. py__init_dim_obs_pdaf " \
                               "4. py__obs_op_pdaf (for ensemble mean)" \
                               "5. py__init_obs_pdaf " \
                               "6. py__obs_op_pdaf (for each ensemble member)" \
                               "7. py__init_obsvar_pdaf " \
                               "(only relevant for adaptive forgetting factor schemes)" \
                               "8. py__prodRinvA_pdaf " \
                               "9. core DA algorithm "
docstrings['put_state_generate_obs'] = "This function generates synthetic observations based on each member of model forecast. " \
                                       "This function is for the case where the ensemble size is larger than the number of processors. \n\n" \
                                       "The function executes the user-supplied function " \
                                       "in the following sequence: \n" \
                                       "1. py__collect_state_pdaf " \
                                       "2. py__prepoststep_state_pdaf " \
                                       "3. py__init_dim_obs_pdaf " \
                                       "4. py__obs_op_pdaf" \
                                       "7. py__init_obserr_f_pdaf " \
                                       "5. py__get_obs_f_pdaf " \
                                       "6. py__prepoststep_state_pdaf " \
                                       "7. py__distribute_state_pdaf " \
                                       "8. py__next_observation_pdaf "
docstrings['put_state_hyb3dvar_estkf'] = "Using Hybrid 3DEnVar for DA without OMI. " \
										  "Here, the background error covariance is hybridised by a static background error covariance, " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "ESTKF in this implementation. \n\n" \
                                          "This function is usually used in 'flexible' parallelisation. " \
                                          "i.e., the ensemble size is larger than the available number of processes. " \
                                          "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                          "to the model after this function. " \
                                          "This function should be called at each model time step. \n\n" \
                                          "The function executes the user-supplied function " \
                                          "in the following sequence: \n" \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:" \
			                              "6. py__cvt_pdaf " \
			                              "7. py__cvt_ens_pdaf " \
			                              "8. py__obs_op_lin_pdaf " \
			                              "9. py__prodRinvA_pdaf " \
			                              "10. py__obs_op_adj_pdaf " \
			                              "11. py__cvt_adj_pdaf " \
			                              "12. py__cvt_adj_ens_pdaf " \
			                              "13. core 3DEnVar algorithm " \
			                              "After the iterations: " \
			                              "14. py__cvt_pdaf " \
			                              "15. py__cvt_ens_pdaf " \
			                              "Perform ESTKF: " \
			                              "16. py__init_dim_obs_pdaf " \
			                              "17. py__obs_op_pdaf (for ensemble mean)" \
			                              "18. py__init_obs_pdaf " \
			                              "19. py__obs_op_pdaf (for each ensemble member)" \
			                              "20. py__init_obsvar_pdaf " \
			                              "(only relevant for adaptive forgetting factor schemes)" \
			                              "21. py__prodRinvA_pdaf " \
			                              "22. core ESTKF algorithm "
docstrings['put_state_hyb3dvar_lestkf'] = "Using Hybrid 3DEnVar for DA without OMI. " \
										  "Here, the background error covariance is hybridised by a static background error covariance, " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "LESTKF in this implementation. \n\n" \
                                          "This function is usually used in 'flexible' parallelisation. " \
                                          "i.e., the ensemble size is larger than the available number of processes. " \
                                          "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                          "to the model after this function. " \
                                          "This function should be called at each model time step. \n\n" \
                                          "The function executes the user-supplied function " \
                                          "in the following sequence: \n" \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:" \
			                              "6. py__cvt_pdaf " \
			                              "7. py__cvt_ens_pdaf " \
			                              "8. py__obs_op_lin_pdaf " \
			                              "9. py__prodRinvA_pdaf " \
			                              "10. py__obs_op_adj_pdaf " \
			                              "11. py__cvt_adj_pdaf " \
			                              "12. py__cvt_adj_ens_pdaf " \
			                              "13. core DA algorithm " \
			                              "After the iterations: " \
			                              "14. py__cvt_pdaf " \
			                              "15. py__cvt_ens_pdaf " \
			                              "Perform LESTKF: " \
			                              "16. py__init_n_domains_p_pdaf " \
			                              "17. py__init_dim_obs_pdaf " \
			                              "18. py__obs_op_pdaf (for each ensemble member)" \
			                              "19. py__init_obs_pdaf " \
			                              "(if global adaptive forgetting factor is used " \
			                              "(type_forget=1 in pyPDAF.PDAF.init)) " \
			                              "20. py__init_obsvar_pdaf " \
			                              "(if global adaptive forgetting factor is used) " \
			                              "loop over each local domain:" \
			                              "21. py__init_dim_l_pdaf " \
			                              "22. py__init_dim_obs_l_pdaf " \
			                              "23. py__g2l_state_pdaf " \
			                              "24. py__g2l_obs_pdaf " \
			                              "(localise mean ensemble in observation space)" \
			                              "25. py__init_obs_l_pdaf "\
			                              "26. py__g2l_obs_pdaf " \
			                              "(localise each ensemble member in observation space) " \
			                              "27. py__init_obsvar_l_pdaf " \
			                              "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
			                              "28. py__prodRinvA_pdaf " \
			                              "29. core DA algorithm " \
			                              "30. py__l2g_state_pdaf "
docstrings['put_state_lenkf'] = "Using stochastic EnKF (ensemble " \
                                "Kalman filter) with covariance localisation " \
                                "for DA without OMI. " \
                                "This is the only scheme for covariance localisation in PDAF. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n\n" \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n" \
                                "1. py__collect_state_pdaf " \
                                "2. py__prepoststep_state_pdaf " \
                                "3. py__init_dim_obs_pdaf " \
                                "4. py__obs_op_pdaf (for each ensemble member)" \
                                "5. py__localize_pdaf " \
                                "6. py__add_obs_err_pdaf " \
                                "7. py__init_obs_pdaf " \
                                "8. py__init_obscovar_pdaf " \
                                "9. py__obs_op_pdaf (repeated to reduce storage)" \
                                "10. core DA algorithm"
docstrings['put_state_lestkf'] = "Using Local ESTKF (error space transform " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. " \
                                 "The LESTKF is a more efficient equivalent to the LETKF. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "14. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "15. py__prodRinvA_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf "
docstrings['put_state_letkf'] = "Using local ensemble transform Kalman filter for DA without OMI. " \
                                "This is a domain localisation method. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n\n" \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n" \
                                "1. py__collect_state_pdaf " \
                                "2. py__prepoststep_state_pdaf " \
                                "3. py__init_n_domains_p_pdaf " \
                                "4. py__init_dim_obs_pdaf " \
                                "5. py__obs_op_pdaf (for each ensemble member)" \
                                "6. py__init_obs_pdaf " \
                                "(if global adaptive forgetting factor is used " \
                                "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                "7. py__init_obsvar_pdaf " \
                                "(if global adaptive forgetting factor is used) " \
                                "loop over each local domain:" \
                                "8. py__init_dim_l_pdaf " \
                                "9. py__init_dim_obs_l_pdaf " \
                                "10. py__g2l_state_pdaf " \
                                "11. py__g2l_obs_pdaf " \
                                "(localise mean ensemble in observation space)" \
                                "12. py__init_obs_l_pdaf "\
                                "13. py__g2l_obs_pdaf " \
                                "(localise each ensemble member in observation space) " \
                                "14. py__init_obsvar_l_pdaf " \
                                "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                "15. py__prodRinvA_l_pdaf " \
                                "16. core DA algorithm " \
                                "17. py__l2g_state_pdaf "
docstrings['put_state_lnetf'] = "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. \n\n" \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "loop over each local domain:" \
                                 "6. py__init_dim_l_pdaf " \
                                 "7. py__init_dim_obs_l_pdaf " \
                                 "8. py__g2l_state_pdaf " \
                                 "9. py__init_obs_l_pdaf "\
                                 "10. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)" \
                                 "11. py__likelihood_l_pdaf " \
                                 "12. core DA algorithm " \
                                 "13. py__l2g_state_pdaf "
docstrings['put_state_lknetf'] = "This function will is a hybridised LETKF and LNETF " \
                                 "for DA without OMI. The LNETF computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "The hybridisation with LETKF is expected to lead to improved performance for " \
                                 "quasi-Gaussian problems.  \n\n" \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "14. py__prodRinvA_pdaf " \
                                 "15. py__likelihood_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member)" \
                                 "19. py__likelihood_hyb_l_pdaf" \
                                 "20. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "21. py__prodRinvA_hyb_l_pdaf "
docstrings['put_state_lseik'] = "Using local singular evolutive interpolated Kalman filter for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__init_obs_pdaf " \
                                 "(if global adaptive forgetting factor is used " \
                                 "(type_forget=1 in pyPDAF.PDAF.init)) " \
                                 "7. py__init_obsvar_pdaf " \
                                 "(if global adaptive forgetting factor is used) " \
                                 "loop over each local domain:" \
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)" \
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space) " \
                                 "14. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "15. py__prodRinvA_l_pdaf " \
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf "
docstrings['put_state_netf'] = "This function will use Nonlinear Ensemble Transform Filter (NETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n\n" \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n" \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__init_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "6. py__likelihood_pdaf " \
                                 "7. core DA algorithm "
docstrings['put_state_pf'] = "This function will use particle filter for DA without OMI. " \
                             "This is a fully nonlinear filter. " \
                             "This function is usually used in 'flexible' parallelisation. " \
                             "i.e., the ensemble size is larger than the available number of processes. " \
                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                             "to the model after this function. " \
                             "This function should be called at each model time step. \n\n" \
                             "The function executes the user-supplied function " \
                             "in the following sequence: \n" \
                             "1. py__collect_state_pdaf " \
                             "2. py__prepoststep_state_pdaf " \
                             "3. py__init_dim_obs_pdaf " \
                             "4. py__init_obs_pdaf " \
                             "5. py__obs_op_pdaf (for each ensemble member)" \
                             "6. py__likelihood_pdaf " \
                             "7. core DA algorithm "
docstrings['put_state_prepost'] = "This function does not perform any DA. " \
                                   "It is used to preprocess the ensemble. " \
                                   "To distribute the ensemble back to the model with post-processing, `pyPDAF.PDAF.get_state` function should be used afterwards.\n" \
                                   "The sequence of the user-supplied functions is: \n" \
                                   "1. py__collect_state_pdaf " \
                                   "2. py__prepoststep_state_pdaf "
docstrings['put_state_seek'] = "This function will use singular evolutive extended Kalman filter for DA without OMI. " \
                               "This function is usually used in 'flexible' parallelisation, " \
                               "but SEEK is deterministic and does not require ensemble. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process the state vector and " \
                               "distribute the state vector " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n\n" \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n" \
                               "1. py__collect_state_pdaf " \
                               "2. py__prepoststep_state_pdaf " \
                               "3. py__init_dim_obs_pdaf " \
                               "4. py__obs_op_pdaf (for ensemble mean)" \
                               "5. py__init_obs_pdaf " \
                               "6. py__obs_op_pdaf (for each ensemble member)" \
                               "7. py__prodRinvA_pdaf " \
                               "8. core DA algorithm "
docstrings['put_state_seik'] = "This function will use singular evolutive interpolated Kalman filter for DA without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n\n" \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n" \
                               "1. py__collect_state_pdaf " \
                               "2. py__prepoststep_state_pdaf " \
                               "3. py__init_dim_obs_pdaf " \
                               "4. py__obs_op_pdaf (for ensemble mean)" \
                               "5. py__init_obs_pdaf " \
                               "6. py__obs_op_pdaf (for each ensemble member)" \
                               "7. py__init_obsvar_pdaf " \
                               "(only relevant for adaptive forgetting factor schemes)" \
                               "8. py__prodRinvA_pdaf " \
                               "9. core DA algorithm "

docstrings['deallocate'] = "This function finalise the PDAF systems including deaclloating all arrays in PDAF."

docstrings['diag_effsample'] = "This function calculates the effective sample size of a particle filter " \
                               "as defined in Doucet et al. 2001 p. 333. " \
                               "It is defined as the inverse of the sum of the squared particle filter weights. " \
                               "If the `n_eff=dim_sample`, all weights are identical, the filter has no influence. "  \
                               "If `n_eff=0`, the filter is collapsed. " \
                               "This is typically called during the analysis step of a particle filter, "\
                               "e.g. in the analysis step of NETF and LNETF."
docstrings['diag_ensstats']  = "This function returns the skewness and kurtosis of the ensemble of a given element of the state vector. " \
                               "The definition used for kurtosis follows that used by Lawson and Hansen, Mon. Wea. Rev. 132 (2004) 1966."
docstrings['diag_histogram'] = "This function returns a rank histogram of the ensemble. " \
                               "A rank histogram is used to diagnose the reliability of the ensemble. " \
                               "A perfectly reliable ensemble should have a uniform rank histogram. " \
                               "The function can be called in the pre/poststep routine of PDAF "\
                               "both before and after the analysis step to collect the histogram information."

docstrings['eofcovar'] = "This function performs an EOF analysis of a state vectors at multiple time steps " \
                         "by singular value decomposition. " \
                         "A multivariate scaling can be performed to ensure that " \
                         "all fields in the state vectors have unit variance. " \
                         "The function returns a the singular vectors and square root of the singular values of the covariance matrix. " \
                         "It also returns the time mean of the state vectors and the temporal anomaly used for the EOF analysis.\n\n" \
                         "These information can be used to initialise an ensemble. \n\n" \
                         "One can store the singular vectors and corresponding values in a file. " \
                         "When one wants to initialise an ensemble from these EOFs for " \
                         "a data assimilation application, one can use the function " \
                         "`pyPDAF.PDAF.sampleens` to generate an ensemble of a chosen size " \
                         "(up to the number of EOFs plus one) " \
                         "by second-order exact sampling. " \
                         "Thus, it can be useful to store more EOFs than one finally " \
                         "might want to use to have the flexibility to cary the ensemble size."

docstrings['gather_dim_obs_f'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF) " \
                                 "this function returns the total observation dimension " \
                                 "from process-local observation dimensions. " \
                                 "The function has to be called once before using any of the functions " \
                                 "`pyPDAF.PDAF.gather_obs_f` or `pyPDAF.PDAF.gather_obs_f2`. " \
                                 "This is because it stores the information on the process-local observation dimensions " \
                                 "to allocate actual observation vectors. " \
                                 "The routine is typically used in the routine `py__init_dim_obs_f_pdaf` " \
                                 "if the analysis step of the local filters is parallelized."

docstrings['gather_obs_f'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF) " \
                             "this function returns the total observation vector " \
                             "from process-local observations. " \
                             "The function depends on " \
                             "`pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. " \
                             "Further, the related routine `pyPDAF.PDAF.gather_obs_f2` is used to " \
                             "gather the associated 2D observation coordinates."

docstrings['gather_obs_f2'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF) " \
                             "this function returns the full observation coordinates " \
                             "from process-local observation coordinates. " \
                             "The function depends on " \
                             "`pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. " \
                             "Further, the related routine `pyPDAF.PDAF.gather_obs_f` is used to " \
                             "gather the associated observation vectors. \n\n" \
                             "The routine is typically used in the routines `py__init_dim_obs_f_pdaf` " \
                             "if the analysis step of the local filters is parallelized."

docstrings['get_assim_flag'] = "This function returns the flag that indicates if the DA is performed in the last time step. " \
                               "It only works for online DA systems. "

docstrings['get_ensstats'] = "This is a diagnotics function for LKNETF which returns the skewness and kutosis used there. "

docstrings['get_localfilter'] = "This function returns whether a local filter is used. "

docstrings['get_memberid'] = "This function returns the ensemble member id on the current process. \n" \
                             "For example, it can be called during the ensemble integration if ensemble-specific forcing is applied. " \
                             "It can also be used in the user-supplied functions such as `py__collect_state_pdaf` and `py__distribute_state_pdaf`."

docstrings['get_obsmemberid'] = "This function returns the ensemble member id when observation operator is being applied. \n" \
                                "This function is used specifically for user-supplied function `py__obs_op_pdaf`."

docstrings['get_smootherens'] = "This function returns the smoothed ensemble in earlier time steps. " \
                                "It is only used when the smoother options is used ."

docstrings['get_state'] = "This function distribute the analysis state vector back to the model. " \
                          "It also post-processes the ensemble and sets the next assimilation time. \n\n" \
                          "The function executes the user-supplied function in the following sequence: \n" \
                          "1. py__prepoststep_state_pdaf " \
                          "2. py__distribute_state_pdaf " \
                          "3. py__next_observation_pdaf "

docstrings['init'] = "This function initialises the PDAF system. " \
                     "It is called once at the beginning of the assimilation. " \
                     "The function specifies the type of DA methods, parameters of the filters, the MPI communicators, and other parallel options." \
                     "The function also provides an initial ensemble to the PDAF system by the user-supplied function which can be distribute to the model by `pyPDAF.PDAF.get_state`. \n" \
                     "For the options and parameters of DA methods, check the https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF. \n" \
                     "The parallisation module in the repository example can be used directly for most cases. " \
                     "Explanation of the parallelisation strategy in PDAF can be found in https://pdaf.awi.de/trac/wiki/ImplementationConceptOnline#Parallelizationofthedataassimilationprogram " \
                     "and https://pdaf.awi.de/trac/wiki/AdaptParallelization"

docstrings['local_weight'] = "The function is used for localisation in the analysis step of a filter " \
                             "and computes a weight according to the specified distance " \
                             "and the settings for the localising function. " \
                             "Typically the function is called in `py__prodRinvA_l_pdaf` " \
                             "in the domain-localised filters. " \
                             "Also, the function is typically called for the LEnKF " \
                             "in the `py__localize_covar_pdaf`. \n" \
                             "This function is usually only used in user-codes that do not use PDAF-OMI."

docstrings['print_info'] = None
docstrings['reset_forget'] = None
docstrings['SampleEns'] = None
docstrings['set_debug_flag'] = None
docstrings['set_ens_pointer'] = None
docstrings['set_smootherens'] = None
docstrings['seik_TtimesA'] = None
docstrings['etkf_Tleft'] = None
docstrings['estkf_OmegaA'] = None
docstrings['enkf_omega'] = None
docstrings['seik_omega'] = None
docstrings['incremental'] = None
docstrings['add_increment'] = None
docstrings['local_weights'] = None
docstrings['diag_CRPS'] = None
docstrings['force_analysis'] = None
docstrings['gather_obs_f2_flex'] = None
docstrings['gather_obs_f_flex'] = None
docstrings['prepost'] = None
docstrings['set_memberid'] = None
docstrings['set_comm_pdaf'] = None
docstrings['set_offline_mode'] = None
docstrings['print_domain_stats'] = None
docstrings['init_local_obsstats'] = None
docstrings['incr_local_obsstats'] = None
docstrings['print_local_obsstats'] = None


docstrings['omit_obs_omi'] = None
