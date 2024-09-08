docstrings = dict()
docstrings['assimilate_estkf'] = "Using ESTKF (error space transform " \
                                 "Kalman filter) for DA without OMI." \
                                 "This function should be called at each model time step." \
                                 "The ESTKF is a more efficient equivalent to the ETKF"
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "Kalman filter) for DA without OMI." \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step." \
                                 "The LESTKF is a more efficient equivalent to the LETKF"
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "loop over each local domain:"
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)"
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
                                 "Kalman filter) for DA without OMI." \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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

docstrings['assimilate_letkf'] = "Using local ensemble transform Kalman filter for DA without OMI." \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "loop over each local domain:"
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)"
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
                                 "Kalman filter) for DA without OMI." \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "for DA without OMI." \
                                 "This is the only scheme for covariance localisation in PDAF". \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf (for each ensemble member)" \
                                 "5. py__localize_pdaf "
                                 "6. py__add_obs_err_pdaf " \
                                 "7. py__init_obs_pdaf " \
                                 "8. py__init_obscovar_pdaf " \
                                 "9. py__obs_op_pdaf (repeated to reduce storage)" \
                                 "10. core DA algorithm"
                                 "11. py__prepoststep_state_pdaf " \
                                 "12. py__distribute_state_pdaf " \
                                 "13. py__next_observation_pdaf "

docstrings['assimilate_3dvar'] = "Using 3DVar for DA without OMI." \
								 "This is a deterministic filtering scheme " \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:"
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

docstrings['assimilate_en3dvar_estkf'] = "Using 3DEnVar for DA without OMI." \
								 "The background error covariance matrix is estimated by ensemble " \
								 "The 3DEnVar only calculates the ensemble mean " \
								 "An ESTKF is used to generate ensemble perturbations. " \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:"
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

docstrings['assimilate_en3dvar_lestkf'] = "Using 3DEnVar for DA without OMI." \
								 "The background error covariance matrix is estimated by ensemble " \
								 "The 3DEnVar only calculates the ensemble mean " \
								 "An ESTKF is used to generate ensemble perturbations. " \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_dim_obs_pdaf " \
                                 "4. py__obs_op_pdaf " \
                                 "5. py__init_obs_pdaf " \
                                 "Starting the iterative optimisation:"
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
                                 "loop over each local domain:"
                                 "18. py__init_dim_l_pdaf " \
                                 "19. py__init_dim_obs_l_pdaf " \
                                 "20. py__g2l_state_pdaf " \
                                 "21. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)"
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
										  "Here, the background error covariance is hybridised by a static background error covariance " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "ESTKF in this implementation."
			                              "The function executes the user-supplied function " \
			                              "in the following sequence: " \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:"
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
										  "Here, the background error covariance is hybridised by a static background error covariance " \
										  "and a flow-dependent background error covariance estimated from ensemble. " \
										  "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
										  "LESTKF in this implementation."
			                              "The function executes the user-supplied function " \
			                              "in the following sequence: " \
			                              "1. py__collect_state_pdaf " \
			                              "2. py__prepoststep_state_pdaf " \
			                              "3. py__init_dim_obs_pdaf " \
			                              "4. py__obs_op_pdaf " \
			                              "5. py__init_obs_pdaf " \
			                              "Starting the iterative optimisation:"
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
			                              "loop over each local domain:"
			                              "21. py__init_dim_l_pdaf " \
			                              "22. py__init_dim_obs_l_pdaf " \
			                              "23. py__g2l_state_pdaf " \
			                              "24. py__g2l_obs_pdaf " \
			                              "(localise mean ensemble in observation space)"
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
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
                                 "1. py__collect_state_pdaf " \
                                 "2. py__prepoststep_state_pdaf " \
                                 "3. py__init_n_domains_p_pdaf " \
                                 "4. py__init_dim_obs_pdaf " \
                                 "5. py__obs_op_pdaf (for each ensemble member)" \
                                 "loop over each local domain:"
                                 "6. py__init_dim_l_pdaf " \
                                 "7. py__init_dim_obs_l_pdaf " \
                                 "8. py__g2l_state_pdaf " \
                                 "9. py__init_obs_l_pdaf "\
                                 "10. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)"
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
                                 "quasi-Gaussian problems."
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "loop over each local domain:"
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise each ensemble member in observation space)"
                                 "12. py__init_obs_l_pdaf "\
                                 "13. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "14. py__prodRinvA_pdaf " \
                                 "15. py__likelihood_l_pdaf " \                                
                                 "16. core DA algorithm " \
                                 "17. py__l2g_state_pdaf " \
                                 "18. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member)" \
                                 "19. py__likelihood_hyb_l_pdaf"
                                 "20. py__init_obsvar_l_pdaf " \
                                 "(only called if local adaptive forgetting factor (type_forget=2) is used) "\
                                 "21. py__prodRinvA_hyb_l_pdaf " \
                                 "22. py__prepoststep_state_pdaf " \
                                 "23. py__distribute_state_pdaf " \
                                 "24. py__next_observation_pdaf "

docstrings['assimilate_lseik'] = "Using local singular evolutive interpolated Kalman filter for DA without OMI." \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "loop over each local domain:"
                                 "8. py__init_dim_l_pdaf " \
                                 "9. py__init_dim_obs_l_pdaf " \
                                 "10. py__g2l_state_pdaf " \
                                 "11. py__g2l_obs_pdaf " \
                                 "(localise mean ensemble in observation space)"
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
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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
                                 "The function should be called at each model step." \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: " \
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

docstrings['assimilate_prepost'] = "This function does not "