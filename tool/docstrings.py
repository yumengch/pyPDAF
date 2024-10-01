docstrings = dict()
docstrings['assimilate_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_global_nondiagR`. \n    " \
                                 "Using ESTKF (error space transform " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This function should be called at each model time step. " \
                                 "The ESTKF is a more efficient equivalent to the ETKF. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_estkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf \n    " \
                                 "2. py__prepoststep_state_pdaf \n    " \
                                 "3. py__init_dim_obs_pdaf \n    " \
                                 "4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                 "5. py__init_obs_pdaf \n    " \
                                 "6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes) \n    " \
                                 "8. py__prodRinvA_pdaf \n    " \
                                 "9. core DA algorithm \n    " \
                                 "10. py__prepoststep_state_pdaf \n    " \
                                 "11. py__distribute_state_pdaf \n    " \
                                 "12. py__next_observation_pdaf "
docstrings['assimilate_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                  "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                  "Using Local ESTKF (error space transform " \
                                  "Kalman filter) for DA without OMI. " \
                                  "This is a domain localisation method. " \
                                  "This function should be called at each model time step. " \
                                  "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                  "The function is a combination of `pyPDAF.PDAF.put_state_lestkf` " \
                                  "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                  "in the following sequence: \n    " \
                                  "1. py__collect_state_pdaf \n    "\
                                  "2. py__prepoststep_state_pdaf \n    "\
                                  "3. py__init_n_domains_p_pdaf \n    "\
                                  "4. py__init_dim_obs_pdaf \n    "\
                                  "5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                  "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init)) \n    "\
                                  "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used) \n    "\
                                  "loop over each local domain:\n    " \
                                  "8. py__init_dim_l_pdaf \n    "\
                                  "9. py__init_dim_obs_l_pdaf \n    "\
                                  "10. py__g2l_state_pdaf \n    "\
                                  "11. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                  "12. py__init_obs_l_pdaf \n    " \
                                  "13. py__g2l_obs_pdaf (localise each ensemble member in observation space) \n    "\
                                  "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used) \n    " \
                                  "15. py__prodRinvA_l_pdaf \n    "\
                                  "16. core DA algorithm \n    " \
                                  "17. py__l2g_state_pdaf \n    "\
                                  "18. py__prepoststep_state_pdaf \n    "\
                                  "19. py__distribute_state_pdaf \n    "\
                                  "20. py__next_observation_pdaf "
docstrings['assimilate_etkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e. `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_global_nondiagR`. \n    " \
                                "Using ETKF (ensemble transform " \
                                "Kalman filter) for DA without OMI. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function is a combination of `pyPDAF.PDAF.put_state_etkf` " \
                                "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                "5. py__init_obs_pdaf\n    " \
                                "6. py__obs_op_pdaf (for each ensemble member\n    " \
                                "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes\n    " \
                                "8. py__prodRinvA_pdaf\n    " \
                                "9. core DA algorithm\n    " \
                                "10. py__prepoststep_state_pdaf\n    " \
                                "11. py__distribute_state_pdaf\n    " \
                                "12. py__next_observation_pdaf\n    "
docstrings['assimilate_letkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                 "Using local ensemble transform Kalman filter for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_letkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_n_domains_p_pdaf\n    " \
                                 "4. py__init_dim_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n" \
                                 "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init))\n    " \
                                 "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                 "loop over each local domain:\n    " \
                                 "8. py__init_dim_l_pdaf\n    " \
                                 "9. py__init_dim_obs_l_pdaf\n    " \
                                 "10. py__g2l_state_pdaf\n    " \
                                 "11. py__g2l_obs_pdaf (localise mean ensemble in observation space\n    " \
                                 "12. py__init_obs_l_pdaf\n    "\
                                 "13. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                 "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)\n    "\
                                 "15. py__prodRinvA_l_pdaf\n    " \
                                 "16. core DA algorithm\n    " \
                                 "17. py__l2g_state_pdaf\n    " \
                                 "18. py__prepoststep_state_pdaf\n    " \
                                 "19. py__distribute_state_pdaf\n    " \
                                 "20. py__next_observation_pdaf\n    "
docstrings['assimilate_enkf'] =  "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`. \n    " \
                                 "Using stochastic EnKF (ensemble " \
                                 "Kalman filter) for DA without OMI. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_enkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_dim_obs_pdaf\n    " \
                                 "4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                 "5. py__add_obs_err_pdaf\n    " \
                                 "6. py__init_obs_pdaf\n    " \
                                 "7. py__init_obscovar_pdaf\n    " \
                                 "8. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "9. core DA algorithm\n    " \
                                 "10. py__prepoststep_state_pdaf\n    " \
                                 "11. py__distribute_state_pdaf\n    " \
                                 "12. py__next_observation_pdaf\n    "
docstrings['assimilate_lenkf'] = "It is recommended to OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.omi_assimilate_lenkf` or `pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`. \n    " \
                                 "Using stochastic EnKF (ensemble " \
                                 "Kalman filter) with covariance localisation " \
                                 "for DA without OMI. " \
                                 "This is the only scheme for covariance localisation in PDAF. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lenkf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_dim_obs_pdaf\n    " \
                                 "4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "5. py__localize_pdaf\n    " \
                                 "6. py__add_obs_err_pdaf\n    " \
                                 "7. py__init_obs_pdaf\n    " \
                                 "8. py__init_obscovar_pdaf\n    " \
                                 "9. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                 "10. core DA algorith\n    " \
                                 "11. py__prepoststep_state_pdaf\n    " \
                                 "12. py__distribute_state_pdaf\n    " \
                                 "13. py__next_observation_pdaf\n    "
docstrings['assimilate_3dvar'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency." \
                                 "I.e., `pyPDAF.PDAF.omi_assimilate_3dvar` or `pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`. \n    " \
                                 "Using 3DVar for DA without OMI. " \
                                 "This is a deterministic filtering scheme. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_3dvar` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_dim_obs_pdaf\n    " \
                                 "4. py__obs_op_pdaf\n    " \
                                 "5. py__init_obs_pdaf\n    " \
                                 "Starting the iterative optimisation:\n    " \
                                 "6. py__cvt_pdaf\n    " \
                                 "7. py__obs_op_lin_pdaf\n    " \
                                 "8. py__prodRinvA_pdaf\n    " \
                                 "9. py__obs_op_adj_pdaf\n    " \
                                 "10. py__cvt_adj_pdaf\n    " \
                                 "11. core DA algorithm\n    " \
                                 "After the iterations: \n    " \
                                 "12. py__cvt_pdaf\n    " \
                                 "13. py__prepoststep_state_pdaf\n    " \
                                 "14. py__distribute_state_pdaf\n    " \
                                 "15. py__next_observation_pdaf\n    "
docstrings['assimilate_en3dvar_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                         "I.e., `pyPDAF.PDAF.omi_assimilate_en3dvar_estkf` or `pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`. \n    " \
                                         "Using 3DEnVar for DA without OMI. " \
                                         "The background error covariance matrix is estimated by ensemble. " \
                                         "The 3DEnVar only calculates the analysis of the ensemble mean. " \
                                         "An ESTKF is used to generate ensemble perturbations. " \
                                         "This function should be called at each model time step. \n    \n    " \
                                         "The function is a combination of `pyPDAF.PDAF.put_state_en3dvar_estkf` " \
                                         "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                         "in the following sequence: \n    " \
                                         "1. py__collect_state_pdaf\n    " \
                                         "2. py__prepoststep_state_pdaf\n    " \
                                         "3. py__init_dim_obs_pdaf\n    " \
                                         "4. py__obs_op_pdaf\n    " \
                                         "5. py__init_obs_pdaf\n    " \
                                         "Starting the iterative optimisation:\n    " \
                                         "6. py__cvt_ens_pdaf\n    " \
                                         "7. py__obs_op_lin_pdaf\n    " \
                                         "8. py__prodRinvA_pdaf\n    " \
                                         "9. py__obs_op_adj_pdaf\n    " \
                                         "10. py__cvt_adj_ens_pdaf\n    " \
                                         "11. core 3DEnVar algorithm\n    " \
                                         "After the iterations: \n    " \
                                         "12. py__cvt_ens_pdaf\n    " \
                                         "Perform ESTKF: " \
                                         "13. py__init_dim_obs_pdaf\n    " \
                                         "14. py__obs_op_pdaf (for ensemble mean)\n    " \
                                         "15. py__init_obs_pdaf\n    " \
                                         "16. py__obs_op_pdaf (for each ensemble member)\n    " \
                                         "17. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                         "18. py__prodRinvA_pdaf\n    " \
                                         "19. core ESTKF algorithm\n    " \
                                         "20. py__prepoststep_state_pdaf\n    " \
                                         "21. py__distribute_state_pdaf\n    " \
                                         "22. py__next_observation_pdaf\n    "
docstrings['assimilate_en3dvar_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                          "I.e., `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf` or `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`. \n    " \
                                          "Using 3DEnVar for DA without OMI. " \
                                          "The background error covariance matrix is estimated by ensemble. " \
                                          "The 3DEnVar only calculates the analysis of the ensemble mean. " \
                                          "An LESTKF is used to generate ensemble perturbations. " \
                                          "This function should be called at each model time step. \n    \n    " \
                                          "The function is a combination of `pyPDAF.PDAF.put_state_en3dvar_lestkf` " \
                                          "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                          "in the following sequence: \n    " \
                                          "1. py__collect_state_pdaf\n    " \
                                          "2. py__prepoststep_state_pdaf\n    " \
                                          "3. py__init_dim_obs_pdaf\n    " \
                                          "4. py__obs_op_pdaf\n    " \
                                          "5. py__init_obs_pdaf\n    " \
                                          "Starting the iterative optimisation:\n    " \
                                          "6. py__cvt_ens_pdaf\n    " \
                                          "7. py__obs_op_lin_pdaf\n    " \
                                          "8. py__prodRinvA_pdaf\n    " \
                                          "9. py__obs_op_adj_pdaf\n    " \
                                          "10. py__cvt_adj_ens_pdaf\n    " \
                                          "11. core DA algorithm\n    " \
                                          "After the iterations: \n    " \
                                          "12. py__cvt_ens_pdaf\n    " \
                                          "Perform LESTKF: \n    " \
                                          "13. py__init_n_domains_p_pdaf\n    " \
                                          "14. py__init_dim_obs_pdaf\n    " \
                                          "15. py__obs_op_pdaf (for each ensemble member\n    " \
                                          "16. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init))(if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init)(if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init))(if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init))(if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init)(if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                          "17. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                          "loop over each local domain:\n    " \
                                          "18. py__init_dim_l_pdaf\n    " \
                                          "19. py__init_dim_obs_l_pdaf\n    " \
                                          "20. py__g2l_state_pdaf\n    " \
                                          "21. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                          "22. py__init_obs_l_pdaf\n    "\
                                          "23. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                          "24. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)(only called if local adaptive forgetting factor (type_forget=2\n    "\
                                          "25. py__prodRinvA_l_pdaf\n    " \
                                          "26. core DA algorithm\n    " \
                                          "27. py__l2g_state_pdaf\n    " \
                                          "28. py__prepoststep_state_pdaf\n    " \
                                          "29. py__distribute_state_pdaf\n    " \
                                          "30. py__next_observation_pdaf\n    "
docstrings['assimilate_hyb3dvar_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                          "I.e., `pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf` or `pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`. \n    " \
                                          "Using Hybrid 3DEnVar for DA without OMI. " \
                                          "Here, the background error covariance is hybridised by a static background error covariance, " \
                                          "and a flow-dependent background error covariance estimated from ensemble. " \
                                          "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
                                          "ESTKF in this implementation. " \
                                          "This function should be called at each model time step. \n    \n    " \
                                          "The function is a combination of `pyPDAF.PDAF.put_state_hyb3dvar_estkf` " \
                                          "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                          "in the following sequence: \n    " \
                                          "1. py__collect_state_pdaf\n    " \
                                          "2. py__prepoststep_state_pdaf\n    " \
                                          "3. py__init_dim_obs_pdaf\n    " \
                                          "4. py__obs_op_pdaf\n    " \
                                          "5. py__init_obs_pdaf\n    " \
                                          "Starting the iterative optimisation:\n    " \
                                          "6. py__cvt_pdaf\n    " \
                                          "7. py__cvt_ens_pdaf\n    " \
                                          "8. py__obs_op_lin_pdaf\n    " \
                                          "9. py__prodRinvA_pdaf\n    " \
                                          "10. py__obs_op_adj_pdaf\n    " \
                                          "11. py__cvt_adj_pdaf\n    " \
                                          "12. py__cvt_adj_ens_pdaf\n    " \
                                          "13. core 3DEnVar algorithm\n    " \
                                          "After the iterations: \n    " \
                                          "14. py__cvt_pdaf\n    " \
                                          "15. py__cvt_ens_pdaf\n    " \
                                          "Perform ESTKF: " \
                                          "16. py__init_dim_obs_pdaf\n    " \
                                          "17. py__obs_op_pdaf (for ensemble mean\n    " \
                                          "18. py__init_obs_pdaf\n    " \
                                          "19. py__obs_op_pdaf (for each ensemble member\n    " \
                                          "20. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                          "21. py__prodRinvA_pdaf\n    " \
                                          "22. core ESTKF algorithm\n    " \
                                          "23. py__prepoststep_state_pdaf\n    " \
                                          "24. py__distribute_state_pdaf\n    " \
                                          "25. py__next_observation_pdaf\n    "
docstrings['assimilate_hyb3dvar_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                           "I.e., `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf` or `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`. \n    " \
                                           "Using Hybrid 3DEnVar for DA without OMI. " \
                                           "Here, the background error covariance is hybridised by a static background error covariance, " \
                                           "and a flow-dependent background error covariance estimated from ensemble. " \
                                           "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by " \
                                           "LESTKF in this implementation. " \
                                           "This function should be called at each model time step. \n    \n    " \
                                           "The function is a combination of `pyPDAF.PDAF.put_state_hyb3dvar_lestkf` " \
                                           "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                           "in the following sequence: \n    " \
                                           "1. py__collect_state_pdaf\n    " \
                                           "2. py__prepoststep_state_pdaf\n    " \
                                           "3. py__init_dim_obs_pdaf\n    " \
                                           "4. py__obs_op_pdaf\n    " \
                                           "5. py__init_obs_pdaf\n    " \
                                           "Starting the iterative optimisation:\n    " \
                                           "6. py__cvt_pdaf\n    " \
                                           "7. py__cvt_ens_pdaf\n    " \
                                           "8. py__obs_op_lin_pdaf\n    " \
                                           "9. py__prodRinvA_pdaf\n    " \
                                           "10. py__obs_op_adj_pdaf\n    " \
                                           "11. py__cvt_adj_pdaf\n    " \
                                           "12. py__cvt_adj_ens_pdaf\n    " \
                                           "13. core DA algorithm\n    " \
                                           "After the iterations: \n    " \
                                           "14. py__cvt_pdaf\n    " \
                                           "15. py__cvt_ens_pdaf\n    " \
                                           "Perform LESTKF: \n    " \
                                           "16. py__init_n_domains_p_pdaf\n    " \
                                           "17. py__init_dim_obs_pdaf\n    " \
                                           "18. py__obs_op_pdaf (for each ensemble member\n    " \
                                           "19. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in `pyPDAF.PDAF.init`))\n    " \
                                           "20. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                           "loop over each local domain:\n    " \
                                           "21. py__init_dim_l_pdaf\n    " \
                                           "22. py__init_dim_obs_l_pdaf\n    " \
                                           "23. py__g2l_state_pdaf\n    " \
                                           "24. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                           "25. py__init_obs_l_pdaf\n    "\
                                           "26. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                           "27. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)\n    "\
                                           "28. py__prodRinvA_l_pdaf\n    " \
                                           "29. core DA algorithm\n    " \
                                           "30. py__l2g_state_pdaf\n    " \
                                           "31. py__prepoststep_state_pdaf\n    " \
                                           "32. py__distribute_state_pdaf\n    " \
                                           "33. py__next_observation_pdaf\n    "
docstrings['assimilate_lnetf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`. \n    " \
                                 "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lnetf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_n_domains_p_pdaf\n    " \
                                 "4. py__init_dim_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                 "loop over each local domain:\n    " \
                                 "6. py__init_dim_l_pdaf\n    " \
                                 "7. py__init_dim_obs_l_pdaf\n    " \
                                 "8. py__g2l_state_pdaf\n    " \
                                 "9. py__init_obs_l_pdaf\n    "\
                                 "10. py__g2l_obs_pdaf(localise each ensemble member in observation space)\n    " \
                                 "11. py__likelihood_l_pdaf\n    " \
                                 "12. core DA algorithm\n    " \
                                 "13. py__l2g_state_pdaf\n    " \
                                 "14. py__prepoststep_state_pdaf\n    " \
                                 "15. py__distribute_state_pdaf\n    " \
                                 "16. py__next_observation_pdaf\n    "
docstrings['assimilate_lknetf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                  "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`. \n    " \
                                  "This function will is a hybridised LETKF and LNETF " \
                                  "for DA without OMI. The LNETF computes the distribution up to " \
                                  "the second moment similar to KF but using a nonlinear weighting similar to " \
                                  "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                  "The hybridisation with LETKF is expected to lead to improved performance for " \
                                  "quasi-Gaussian problems. " \
                                  "The function should be called at each model step. \n    \n    " \
                                  "The function is a combination of `pyPDAF.PDAF.put_state_lknetf` " \
                                  "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                  "in the following sequence: \n    " \
                                  "1. py__collect_state_pdaf\n    " \
                                  "2. py__prepoststep_state_pdaf\n    " \
                                  "3. py__init_n_domains_p_pdaf\n    " \
                                  "4. py__init_dim_obs_pdaf\n    " \
                                  "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                  "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                  "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                  "loop over each local domain:\n    " \
                                  "8. py__init_dim_l_pdaf\n    " \
                                  "9. py__init_dim_obs_l_pdaf\n    " \
                                  "10. py__g2l_state_pdaf\n    " \
                                  "11. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                  "12. py__init_obs_l_pdaf\n    "\
                                  "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                  "14. py__prodRinvA_pdaf\n    " \
                                  "15. py__likelihood_l_pdaf\n    " \
                                  "16. core DA algorithm\n    " \
                                  "17. py__l2g_state_pdaf\n    " \
                                  "18. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                  "19. py__likelihood_hyb_l_pda\n    " \
                                  "20. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                  "21. py__prodRinvA_hyb_l_pdaf\n    " \
                                  "22. py__prepoststep_state_pdaf\n    " \
                                  "23. py__distribute_state_pdaf\n    " \
                                  "24. py__next_observation_pdaf\n    "
docstrings['assimilate_lseik'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                 "Using local singular evolutive interpolated Kalman filter for DA without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function should be called at each model time step.\n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_lseik` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_n_domains_p_pdaf\n    " \
                                 "4. py__init_dim_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                 "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                 "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                 "loop over each local domain:\n    " \
                                 "8. py__init_dim_l_pdaf\n    " \
                                 "9. py__init_dim_obs_l_pdaf\n    " \
                                 "10. py__g2l_state_pdaf\n    " \
                                 "11. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                 "12. py__init_obs_l_pdaf\n    "\
                                 "13. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                 "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                 "15. py__prodRinvA_l_pdaf\n    " \
                                 "16. core DA algorithm\n    " \
                                 "17. py__l2g_state_pdaf\n    " \
                                 "18. py__prepoststep_state_pdaf\n    " \
                                 "19. py__distribute_state_pdaf\n    " \
                                 "20. py__next_observation_pdaf\n    "
docstrings['assimilate_netf'] =  "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`. \n    " \
                                 "This function will use Nonlinear Ensemble Transform Filter (NETF) " \
                                 "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "The function should be called at each model step. \n    \n    " \
                                 "The function is a combination of `pyPDAF.PDAF.put_state_netf` " \
                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_dim_obs_pdaf\n    " \
                                 "4. py__init_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                 "6. py__likelihood_pdaf\n    " \
                                 "7. core DA algorithm\n    " \
                                 "8. py__prepoststep_state_pdaf\n    " \
                                 "9. py__distribute_state_pdaf\n    " \
                                 "10. py__next_observation_pdaf\n    "
docstrings['assimilate_pf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                              "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`. \n    " \
                              "This function will use particle filter for DA without OMI. " \
                              "This is a fully nonlinear filter. " \
                              "The function should be called at each model step. \n    \n    " \
                              "The function is a combination of `pyPDAF.PDAF.put_state_pf` " \
                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                              "in the following sequence: \n    " \
                              "1. py__collect_state_pdaf\n    " \
                              "2. py__prepoststep_state_pdaf\n    " \
                              "3. py__init_dim_obs_pdaf\n    " \
                              "4. py__init_obs_pdaf\n    " \
                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                              "6. py__likelihood_pdaf\n    " \
                              "7. core DA algorithm\n    " \
                              "8. py__prepoststep_state_pdaf\n    " \
                              "9. py__distribute_state_pdaf\n    " \
                              "10. py__next_observation_pdaf\n    "
docstrings['assimilate_seek'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_global_nondiagR`. \n    " \
                                "This function will use singular evolutive extended Kalman filter for DA without OMI. " \
                                "This is a deterministic Kalman filter. " \
                                "The function should be called at each model step.\n    \n    " \
                                "The function is a combination of `pyPDAF.PDAF.put_state_seek` " \
                                "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                "5. py__init_obs_pdaf\n    " \
                                "6. py__obs_op_pdaf (for each ensemble member\n    " \
                                "7. py__prodRinvA_pdaf\n    " \
                                "8. core DA algorithm\n    " \
                                "9. py__prepoststep_state_pdaf\n    " \
                                "10. py__distribute_state_pdaf\n    " \
                                "11. py__next_observation_pdaf\n    "
docstrings['assimilate_seik'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.omi_assimilate_global` or `pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`. \n    " \
                                "This function will use singular evolutive interpolated Kalman filter for DA without OMI. " \
                                "The function should be called at each model step.\n    \n    " \
                                "The function is a combination of `pyPDAF.PDAF.put_state_seik` " \
                                "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                "5. py__init_obs_pdaf\n    " \
                                "6. py__obs_op_pdaf (for each ensemble member\n    " \
                                "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                "8. py__prodRinvA_pdaf\n    " \
                                "9. core DA algorithm\n    " \
                                "10. py__prepoststep_state_pdaf\n    " \
                                "11. py__distribute_state_pdaf\n    " \
                                "12. py__next_observation_pdaf\n    "
docstrings['assimilate_prepost'] = "This function does not perform any DA. " \
                                   "It is used to perform a preprocess and postprocess of the ensemble. " \
                                   "Compared to `pyPDAF.PDAF.prepost`, this function sets assimilation flag.\n    " \
                                   "The function is a combination of `pyPDAF.PDAF.put_state_prepost` " \
                                   "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                   "in the following sequence: \n    " \
                                   "1. py__collect_state_pdaf\n    " \
                                   "2. py__prepoststep_state_pdaf (preprocess, step < 0)\n    " \
                                   "3. py__prepoststep_state_pdaf (postprocess, step > 0\n    " \
                                   "4. py__distribute_state_pdaf\n    " \
                                   "5. py__next_observation_pdaf\n    "
docstrings['generate_obs'] = "When diagonal observation error covariance matrix is used, " \
                             "it is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                             "I.e., `pyPDAF.PDAF.omi_generate_obs`. \n    " \
                             "This function generates synthetic observations based on each member of model forecast. " \
                             "This is based on the usual implementation strategy for PDAF.\n    \n    " \
                             "The function is a combination of `pyPDAF.PDAF.put_state_generate_obs` " \
                             "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                             "in the following sequence: \n    " \
                             "1. py__collect_state_pdaf\n    " \
                             "2. py__prepoststep_state_pdaf\n    " \
                             "3. py__init_dim_obs_pdaf\n    " \
                             "4. py__obs_op_pda\n    " \
                             "5. py__init_obserr_f_pdaf\n    " \
                             "6. py__get_obs_f_pdaf\n    " \
                             "7. py__prepoststep_state_pdaf\n    " \
                             "8. py__distribute_state_pdaf\n    " \
                             "9. py__next_observation_pdaf\n    "

docstrings['put_state_3dvar'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                "I.e., `pyPDAF.PDAF.omi_put_state_global` or `pyPDAF.PDAF.omi_put_state_global_nondiagR`. \n    " \
                                "Using 3DVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                "This is a deterministic filtering scheme. " \
                                "This function is usually used in 'flexible' parallelisation, " \
                                "but 3dvar is deterministic and does not require ensemble.\n    " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process the state vector and " \
                                "distribute the state vector " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf\n    " \
                                "5. py__init_obs_pdaf\n    " \
                                "Starting the iterative optimisation:\n    " \
                                "6. py__cvt_pdaf\n    " \
                                "7. py__obs_op_lin_pdaf\n    " \
                                "8. py__prodRinvA_pdaf\n    " \
                                "9. py__obs_op_adj_pdaf\n    " \
                                "10. py__cvt_adj_pdaf\n    " \
                                "11. core DA algorithm\n    " \
                                "After the iterations: \n    " \
                                "12. py__cvt_pdaf\n    "
docstrings['put_state_en3dvar_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                        "I.e., `pyPDAF.PDAF.omi_put_state_en3dvar_estkf` or `pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`. \n    \n    \n    " \
                                        "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                        "The background error covariance matrix is estimated by ensemble. " \
                                        "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                        "An ESTKF is used to generate ensemble perturbations. " \
                                        "This function is usually used in 'flexible' parallelisation. " \
                                        "i.e., the ensemble size is larger than the available number of processes. " \
                                        "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                        "to the model after this function. " \
                                        "This function should be called at each model time step. \n    \n    " \
                                        "The function executes the user-supplied function " \
                                        "in the following sequence: \n    " \
                                        "1. py__collect_state_pdaf\n    " \
                                        "2. py__prepoststep_state_pdaf\n    " \
                                        "3. py__init_dim_obs_pdaf\n    " \
                                        "4. py__obs_op_pdaf\n    " \
                                        "5. py__init_obs_pdaf\n    " \
                                        "Starting the iterative optimisation:\n    " \
                                        "6. py__cvt_ens_pdaf\n    " \
                                        "7. py__obs_op_lin_pdaf\n    " \
                                        "8. py__prodRinvA_pdaf\n    " \
                                        "9. py__obs_op_adj_pdaf\n    " \
                                        "10. py__cvt_adj_ens_pdaf\n    " \
                                        "11. core 3DEnVar algorithm\n    " \
                                        "After the iterations: \n    " \
                                        "12. py__cvt_ens_pdaf\n    " \
                                        "Perform ESTKF: " \
                                        "13. py__init_dim_obs_pdaf\n    " \
                                        "14. py__obs_op_pdaf (for ensemble mean\n    " \
                                        "15. py__init_obs_pdaf\n    " \
                                        "16. py__obs_op_pdaf (for each ensemble member\n    " \
                                        "17. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                        "18. py__prodRinvA_pdaf\n    " \
                                        "19. core ESTKF algorithm\n    "
docstrings['put_state_en3dvar_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                         "I.e., `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf` or `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                         "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                         "The background error covariance matrix is estimated by ensemble. " \
                                         "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                         "An LESTKF is used to generate ensemble perturbations. " \
                                         "This function is usually used in 'flexible' parallelisation. " \
                                         "i.e., the ensemble size is larger than the available number of processes. " \
                                         "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                         "to the model after this function. " \
                                         "This function should be called at each model time step. \n    \n    " \
                                         "The function executes the user-supplied function " \
                                         "in the following sequence: \n    " \
                                         "1. py__collect_state_pdaf\n    " \
                                         "2. py__prepoststep_state_pdaf\n    " \
                                         "3. py__init_dim_obs_pdaf\n    " \
                                         "4. py__obs_op_pdaf\n    " \
                                         "5. py__init_obs_pdaf\n    " \
                                         "Starting the iterative optimisation:\n    " \
                                         "6. py__cvt_ens_pdaf\n    " \
                                         "7. py__obs_op_lin_pdaf\n    " \
                                         "8. py__prodRinvA_pdaf\n    " \
                                         "9. py__obs_op_adj_pdaf\n    " \
                                         "10. py__cvt_adj_ens_pdaf\n    " \
                                         "11. core DA algorithm\n    " \
                                         "After the iterations: \n    " \
                                         "12. py__cvt_ens_pdaf\n    " \
                                         "Perform LESTKF: \n    " \
                                         "13. py__init_n_domains_p_pdaf\n    " \
                                         "14. py__init_dim_obs_pdaf\n    " \
                                         "15. py__obs_op_pdaf (for each ensemble member\n    " \
                                         "16. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                         "17. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                         "loop over each local domain:\n    " \
                                         "18. py__init_dim_l_pdaf\n    " \
                                         "19. py__init_dim_obs_l_pdaf\n    " \
                                         "20. py__g2l_state_pdaf\n    " \
                                         "21. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                         "22. py__init_obs_l_pdaf\n    "\
                                         "23. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                         "24. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                         "25. py__prodRinvA_pdaf\n    " \
                                         "26. core DA algorithm\n    " \
                                         "27. py__l2g_state_pdaf\n    "
docstrings['put_state_enkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                               "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.put_state_enkf_nondiagR`. \n    " \
                               "Using stochastic EnKF (ensemble " \
                               "Kalman filter) for DA without post-processing and analysis distribution to forecsat without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n    \n    " \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n    " \
                               "1. py__collect_state_pdaf\n    " \
                               "2. py__prepoststep_state_pdaf\n    " \
                               "3. py__init_dim_obs_pdaf\n    " \
                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                               "5. py__add_obs_err_pdaf\n    " \
                               "6. py__init_obs_pdaf\n    " \
                               "7. py__init_obscovar_pdaf\n    " \
                               "8. py__obs_op_pdaf (for each ensemble member\n    " \
                               "9. core DA algorithm\n    "
docstrings['put_state_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.put_state_global_nondiagR`. \n    " \
                                "Using ESTKF (error space transform " \
                                "Kalman filter) for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. " \
                                "The ESTKF is a more efficient equivalent to the ETKF. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                "5. py__init_obs_pdaf\n    " \
                                "6. py__obs_op_pdaf (for each ensemble member\n    " \
                                "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                "8. py__prodRinvA_pdaf\n    " \
                                "9. core DA algorithm\n    "
docstrings['put_state_etkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                               "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.put_state_global_nondiagR`. \n    " \
                               "Using ETKF (ensemble transform " \
                               "Kalman filter) for DA without post-processing and analysis distribution to forecsat without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n    \n    " \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n    " \
                               "1. py__collect_state_pdaf\n    " \
                               "2. py__prepoststep_state_pdaf\n    " \
                               "3. py__init_dim_obs_pdaf\n    " \
                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                               "5. py__init_obs_pdaf\n    " \
                               "6. py__obs_op_pdaf (for each ensemble member\n    " \
                               "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                               "8. py__prodRinvA_pdaf\n    " \
                               "9. core DA algorithm\n    "
docstrings['put_state_generate_obs'] = "When diagonal observation error covariance matrix is used, " \
                                       "it is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.omi_put_state_generate_obs`. \n    " \
                                       "This function generates synthetic observations based on each member of model forecast. " \
                                       "This function is for the case where the ensemble size is larger than the number of processors. \n    \n    " \
                                       "The function executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_dim_obs_pdaf\n    " \
                                       "4. py__obs_op_pda\n    " \
                                       "5. py__init_obserr_f_pdaf\n    " \
                                       "6. py__get_obs_f_pdaf\n    " \
                                       "7. py__prepoststep_state_pdaf\n    " \
                                       "8. py__distribute_state_pdaf\n    " \
                                       "9. py__next_observation_pdaf\n    "
docstrings['put_state_hyb3dvar_estkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                         "I.e., `pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf` or `pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`. \n    \n    \n    " \
                                         "Using Hybrid 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                         "Here, the background error covariance is hybridised by a static background error covariance, " \
                                         "and a flow-dependent background error covariance estimated from ensemble. " \
                                         "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                         "ESTKF in this implementation. \n    \n    " \
                                         "This function is usually used in 'flexible' parallelisation. " \
                                         "i.e., the ensemble size is larger than the available number of processes. " \
                                         "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                         "to the model after this function. " \
                                         "This function should be called at each model time step. \n    \n    " \
                                         "The function executes the user-supplied function " \
                                         "in the following sequence: \n    " \
                                         "1. py__collect_state_pdaf\n    " \
                                         "2. py__prepoststep_state_pdaf\n    " \
                                         "3. py__init_dim_obs_pdaf\n    " \
                                         "4. py__obs_op_pdaf\n    " \
                                         "5. py__init_obs_pdaf\n    " \
                                         "Starting the iterative optimisation:\n    " \
                                         "6. py__cvt_pdaf\n    " \
                                         "7. py__cvt_ens_pdaf\n    " \
                                         "8. py__obs_op_lin_pdaf\n    " \
                                         "9. py__prodRinvA_pdaf\n    " \
                                         "10. py__obs_op_adj_pdaf\n    " \
                                         "11. py__cvt_adj_pdaf\n    " \
                                         "12. py__cvt_adj_ens_pdaf\n    " \
                                         "13. core 3DEnVar algorithm\n    " \
                                         "After the iterations: \n    " \
                                         "14. py__cvt_pdaf\n    " \
                                         "15. py__cvt_ens_pdaf\n    " \
                                         "Perform ESTKF: " \
                                         "16. py__init_dim_obs_pdaf\n    " \
                                         "17. py__obs_op_pdaf (for ensemble mean\n    " \
                                         "18. py__init_obs_pdaf\n    " \
                                         "19. py__obs_op_pdaf (for each ensemble member\n    " \
                                         "20. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                         "21. py__prodRinvA_pdaf\n    " \
                                         "22. core ESTKF algorithm\n    "
docstrings['put_state_hyb3dvar_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                          "I.e., `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf` or `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                          "Using Hybrid 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                          "Here, the background error covariance is hybridised by a static background error covariance, " \
                                          "and a flow-dependent background error covariance estimated from ensemble. " \
                                          "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                          "LESTKF in this implementation. \n    \n    " \
                                          "This function is usually used in 'flexible' parallelisation. " \
                                          "i.e., the ensemble size is larger than the available number of processes. " \
                                          "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                          "to the model after this function. " \
                                          "This function should be called at each model time step. \n    \n    " \
                                          "The function executes the user-supplied function " \
                                          "in the following sequence: \n    " \
                                          "1. py__collect_state_pdaf\n    " \
                                          "2. py__prepoststep_state_pdaf\n    " \
                                          "3. py__init_dim_obs_pdaf\n    " \
                                          "4. py__obs_op_pdaf\n    " \
                                          "5. py__init_obs_pdaf\n    " \
                                          "Starting the iterative optimisation:\n    " \
                                          "6. py__cvt_pdaf\n    " \
                                          "7. py__cvt_ens_pdaf\n    " \
                                          "8. py__obs_op_lin_pdaf\n    " \
                                          "9. py__prodRinvA_pdaf\n    " \
                                          "10. py__obs_op_adj_pdaf\n    " \
                                          "11. py__cvt_adj_pdaf\n    " \
                                          "12. py__cvt_adj_ens_pdaf\n    " \
                                          "13. core DA algorithm\n    " \
                                          "After the iterations: \n    " \
                                          "14. py__cvt_pdaf\n    " \
                                          "15. py__cvt_ens_pdaf\n    " \
                                          "Perform LESTKF: \n    " \
                                          "16. py__init_n_domains_p_pdaf\n    " \
                                          "17. py__init_dim_obs_pdaf\n    " \
                                          "18. py__obs_op_pdaf (for each ensemble member\n    " \
                                          "19. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in `pyPDAF.PDAF.init`))\n    " \
                                          "20. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                          "loop over each local domain:\n    " \
                                          "21. py__init_dim_l_pdaf\n    " \
                                          "22. py__init_dim_obs_l_pdaf\n    " \
                                          "23. py__g2l_state_pdaf\n    " \
                                          "24. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                          "25. py__init_obs_l_pdaf\n    "\
                                          "26. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                          "27. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)\n    "\
                                          "28. py__prodRinvA_pdaf\n    " \
                                          "29. core DA algorithm\n    " \
                                          "30. py__l2g_state_pdaf\n    "
docstrings['put_state_lenkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.omi_put_state_lenkf` or `pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`. \n    " \
                                "Using stochastic EnKF (ensemble " \
                                "Kalman filter) with covariance localisation " \
                                "for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                "This is the only scheme for covariance localisation in PDAF. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_dim_obs_pdaf\n    " \
                                "4. py__obs_op_pdaf (for each ensemble member\n    " \
                                "5. py__localize_pdaf\n    " \
                                "6. py__add_obs_err_pdaf\n    " \
                                "7. py__init_obs_pdaf\n    " \
                                "8. py__init_obscovar_pdaf\n    " \
                                "9. py__obs_op_pdaf (repeated to reduce storage\n    " \
                                "10. core DA algorith\n    "
docstrings['put_state_lestkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                 "Using Local ESTKF (error space transform " \
                                 "Kalman filter) for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                 "This is a domain localisation method. " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. " \
                                 "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_n_domains_p_pdaf\n    " \
                                 "4. py__init_dim_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                 "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                 "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                 "loop over each local domain:\n    " \
                                 "8. py__init_dim_l_pdaf\n    " \
                                 "9. py__init_dim_obs_l_pdaf\n    " \
                                 "10. py__g2l_state_pdaf\n    " \
                                 "11. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                 "12. py__init_obs_l_pdaf\n    "\
                                 "13. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                 "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                 "15. py__prodRinvA_l_pdaf\n    " \
                                 "16. core DA algorithm\n    " \
                                 "17. py__l2g_state_pdaf\n    "
docstrings['put_state_letkf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                "Using local ensemble transform Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                "This is a domain localisation method. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_n_domains_p_pdaf\n    " \
                                "4. py__init_dim_obs_pdaf\n    " \
                                "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                "loop over each local domain:\n    " \
                                "8. py__init_dim_l_pdaf\n    " \
                                "9. py__init_dim_obs_l_pdaf\n    " \
                                "10. py__g2l_state_pdaf\n    " \
                                "11. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                "12. py__init_obs_l_pdaf\n    "\
                                "13. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                "15. py__prodRinvA_l_pdaf\n    " \
                                "16. core DA algorithm\n    " \
                                "17. py__l2g_state_pdaf\n    "
docstrings['put_state_lnetf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`. \n    " \
                                "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                "for DA without post-processing and analysis distribution to forecsat without OMI. The nonlinear filter computes the distribution up to " \
                                "the second moment similar to KF but using a nonlinear weighting similar to " \
                                "particle filter. This leads to an equal weights assumption for prior ensemble. \n    \n    " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_n_domains_p_pdaf\n    " \
                                "4. py__init_dim_obs_pdaf\n    " \
                                "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                "loop over each local domain:\n    " \
                                "6. py__init_dim_l_pdaf\n    " \
                                "7. py__init_dim_obs_l_pdaf\n    " \
                                "8. py__g2l_state_pdaf\n    " \
                                "9. py__init_obs_l_pdaf\n    "\
                                "10. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                "11. py__likelihood_l_pdaf\n    " \
                                "12. core DA algorithm\n    " \
                                "13. py__l2g_state_pdaf\n    "
docstrings['put_state_lknetf'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                 "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`. \n    " \
                                 "This function will is a hybridised LETKF and LNETF " \
                                 "for DA without post-processing and analysis distribution to forecsat without OMI. The LNETF computes the distribution up to " \
                                 "the second moment similar to KF but using a nonlinear weighting similar to " \
                                 "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                 "The hybridisation with LETKF is expected to lead to improved performance for " \
                                 "quasi-Gaussian problems.  \n    \n    " \
                                 "This function is usually used in 'flexible' parallelisation. " \
                                 "i.e., the ensemble size is larger than the available number of processes. " \
                                 "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                 "to the model after this function. " \
                                 "This function should be called at each model time step. \n    \n    " \
                                 "The function executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_n_domains_p_pdaf\n    " \
                                 "4. py__init_dim_obs_pdaf\n    " \
                                 "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                 "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                 "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                 "loop over each local domain:\n    " \
                                 "8. py__init_dim_l_pdaf\n    " \
                                 "9. py__init_dim_obs_l_pdaf\n    " \
                                 "10. py__g2l_state_pdaf\n    " \
                                 "11. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                 "12. py__init_obs_l_pdaf\n    "\
                                 "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                 "14. py__prodRinvA_pdaf\n    " \
                                 "15. py__likelihood_l_pdaf\n    " \
                                 "16. core DA algorithm\n    " \
                                 "17. py__l2g_state_pdaf\n    " \
                                 "18. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                 "19. py__likelihood_hyb_l_pda\n    " \
                                 "20. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                 "21. py__prodRinvA_hyb_l_pdaf\n    "
docstrings['put_state_lseik'] = "It is recommended to use local module with OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                "Using local singular evolutive interpolated Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                "This is a domain localisation method. " \
                                "This function is usually used in 'flexible' parallelisation. " \
                                "i.e., the ensemble size is larger than the available number of processes. " \
                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                "to the model after this function. " \
                                "This function should be called at each model time step. \n    \n    " \
                                "The function executes the user-supplied function " \
                                "in the following sequence: \n    " \
                                "1. py__collect_state_pdaf\n    " \
                                "2. py__prepoststep_state_pdaf\n    " \
                                "3. py__init_n_domains_p_pdaf\n    " \
                                "4. py__init_dim_obs_pdaf\n    " \
                                "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                "loop over each local domain:\n    " \
                                "8. py__init_dim_l_pdaf\n    " \
                                "9. py__init_dim_obs_l_pdaf\n    " \
                                "10. py__g2l_state_pdaf\n    " \
                                "11. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                "12. py__init_obs_l_pdaf\n    "\
                                "13. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                "14. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                "15. py__prodRinvA_l_pdaf\n    " \
                                "16. core DA algorithm\n    " \
                                "17. py__l2g_state_pdaf\n    "
docstrings['put_state_netf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                               "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`. \n    " \
                               "This function will use Nonlinear Ensemble Transform Filter (NETF) " \
                               "for DA without post-processing and analysis distribution to forecsat without OMI. The nonlinear filter computes the distribution up to " \
                               "the second moment similar to KF but using a nonlinear weighting similar to " \
                               "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n    \n    " \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n    " \
                               "1. py__collect_state_pdaf\n    " \
                               "2. py__prepoststep_state_pdaf\n    " \
                               "3. py__init_dim_obs_pdaf\n    " \
                               "4. py__init_obs_pdaf\n    " \
                               "5. py__obs_op_pdaf (for each ensemble member\n    " \
                               "6. py__likelihood_pdaf\n    " \
                               "7. core DA algorithm\n    "
docstrings['put_state_pf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                             "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`. \n    " \
                             "This function will use particle filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                             "This is a fully nonlinear filter. " \
                             "This function is usually used in 'flexible' parallelisation. " \
                             "i.e., the ensemble size is larger than the available number of processes. " \
                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                             "to the model after this function. " \
                             "This function should be called at each model time step. \n    \n    " \
                             "The function executes the user-supplied function " \
                             "in the following sequence: \n    " \
                             "1. py__collect_state_pdaf\n    " \
                             "2. py__prepoststep_state_pdaf\n    " \
                             "3. py__init_dim_obs_pdaf\n    " \
                             "4. py__init_obs_pdaf\n    " \
                             "5. py__obs_op_pdaf (for each ensemble member\n    " \
                             "6. py__likelihood_pdaf\n    " \
                             "7. core DA algorithm\n    "
docstrings['put_state_prepost'] = "This function does not perform any DA. " \
                                   "It is used to preprocess the ensemble. " \
                                   "To distribute the ensemble back to the model with post-processing, `pyPDAF.PDAF.get_state` function should be used afterwards.\n    " \
                                   "The sequence of the user-supplied functions is: \n    " \
                                   "1. py__collect_state_pdaf\n    " \
                                   "2. py__prepoststep_state_pdaf\n    "
docstrings['put_state_seek'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                               "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.put_state_global_nondiagR`. \n    " \
                               "This function will use singular evolutive extended Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                               "This function is usually used in 'flexible' parallelisation, " \
                               "but SEEK is deterministic and does not require ensemble. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process the state vector and " \
                               "distribute the state vector " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n    \n    " \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n    " \
                               "1. py__collect_state_pdaf\n    " \
                               "2. py__prepoststep_state_pdaf\n    " \
                               "3. py__init_dim_obs_pdaf\n    " \
                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                               "5. py__init_obs_pdaf\n    " \
                               "6. py__obs_op_pdaf (for each ensemble member\n    " \
                               "7. py__prodRinvA_pdaf\n    " \
                               "8. core DA algorithm\n    "
docstrings['put_state_seik'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                               "I.e., `pyPDAF.PDAF.put_state_global` or `pyPDAF.PDAF.put_state_global_nondiagR`. \n    " \
                               "This function will use singular evolutive interpolated Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                               "This function is usually used in 'flexible' parallelisation. " \
                               "i.e., the ensemble size is larger than the available number of processes. " \
                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                               "to the model after this function. " \
                               "This function should be called at each model time step. \n    \n    " \
                               "The function executes the user-supplied function " \
                               "in the following sequence: \n    " \
                               "1. py__collect_state_pdaf\n    " \
                               "2. py__prepoststep_state_pdaf\n    " \
                               "3. py__init_dim_obs_pdaf\n    " \
                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                               "5. py__init_obs_pdaf\n    " \
                               "6. py__obs_op_pdaf (for each ensemble member\n    " \
                               "7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                               "8. py__prodRinvA_pdaf\n    " \
                               "9. core DA algorithm\n    "

docstrings['deallocate'] = "This function finalise the PDAF systems including deaclloating all arrays in PDAF."

docstrings['diag_effsample'] = "This function calculates the effective sample size of a particle filter " \
                               "as defined in Doucet et al. 2001 p. 333.\n    " \
                               "It is defined as the inverse of the sum of the squared particle filter weights. " \
                               "If the `n_eff=dim_sample`, all weights are identical, the filter has no influence. "  \
                               "If `n_eff=0`, the filter is collapsed. " \
                               "This is typically called during the analysis step of a particle filter, "\
                               "e.g. in the analysis step of NETF and LNETF."
docstrings['diag_ensstats']  = "This function returns the skewness and kurtosis of the ensemble of a given element of the state vector. " \
                               "The definition used for kurtosis follows that used by Lawson and Hansen, Mon. Wea. Rev. 132 (2004) 1966\n    "
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
                         "It also returns the time mean of the state vectors and the temporal anomaly used for the EOF analysis.\n    \n    " \
                         "These information can be used to initialise an ensemble. \n    \n    " \
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
                                 "`pyPDAF.PDAF.gather_obs_f` or `pyPDAF.PDAF.gather_obs_f2`.\n    " \
                                 "This is because it stores the information on the process-local observation dimensions " \
                                 "to allocate actual observation vectors. " \
                                 "The routine is typically used in the routine `py__init_dim_obs_f_pdaf` " \
                                 "if the analysis step of the local filters is parallelized."

docstrings['gather_obs_f'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF) " \
                             "this function returns the total observation vector " \
                             "from process-local observations. " \
                             "The function depends on " \
                             "`pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. " \
                             "Further, the related routine `pyPDAF.PDAF.gather_obs_f2` is used to\n    " \
                             "gather the associated 2D observation coordinates\n    "

docstrings['gather_obs_f2'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF)\n    " \
                             "this function returns the full observation coordinates " \
                             "from process-local observation coordinates. " \
                             "The function depends on " \
                             "`pyPDAF.PDAF.gather_dim_obs_f` which defines the process-local observation dimensions. " \
                             "Further, the related routine `pyPDAF.PDAF.gather_obs_f` is used to " \
                             "gather the associated observation vectors. \n    \n    " \
                             "The routine is typically used in the routines `py__init_dim_obs_f_pdaf` " \
                             "if the analysis step of the local filters is parallelized."

docstrings['get_assim_flag'] = "This function returns the flag that indicates if the DA is performed in the last time step. " \
                               "It only works for online DA systems. "

docstrings['get_ensstats'] = "This is a diagnotics function for LKNETF which returns the skewness and kutosis used there. "

docstrings['get_localfilter'] = "This function returns whether a local filter is used. "

docstrings['get_memberid'] = "This function returns the ensemble member id on the current process. \n    " \
                             "For example, it can be called during the ensemble integration if ensemble-specific forcing is applied. " \
                             "It can also be used in the user-supplied functions such as `py__collect_state_pdaf` and `py__distribute_state_pdaf`."

docstrings['get_obsmemberid'] = "This function returns the ensemble member id when observation operator is being applied. \n    " \
                                "This function is used specifically for user-supplied function `py__obs_op_pdaf`."

docstrings['get_smootherens'] = "This function returns the smoothed ensemble in earlier time steps. " \
                                "It is only used when the smoother options is used ."

docstrings['get_state'] = "This function distribute the analysis state vector back to the model. " \
                          "It also post-processes the ensemble and sets the next assimilation time. \n    \n    " \
                          "The function executes the user-supplied function in the following sequence: \n    " \
                          "1. py__prepoststep_state_pdaf\n    " \
                          "2. py__distribute_state_pdaf\n    " \
                          "3. py__next_observation_pdaf\n    "

docstrings['init'] = "This function initialises the PDAF system. " \
                     "It is called once at the beginning of the assimilation. " \
                     "The function specifies the type of DA methods, parameters of the filters, the MPI communicators, and other parallel options." \
                     "The function also provides an initial ensemble to the PDAF system by the user-supplied function which can be distribute to the model by `pyPDAF.PDAF.get_state`. \n    " \
                     "For the options and parameters of DA methods, check the https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF. \n    " \
                     "The parallisation module in the repository example can be used directly for most cases. " \
                     "Explanation of the parallelisation strategy in PDAF can be found in https://pdaf.awi.de/trac/wiki/ImplementationConceptOnline#Parallelizationofthedataassimilationprogram " \
                     "and https://pdaf.awi.de/trac/wiki/AdaptParallelization"

docstrings['local_weight'] = "The function is used for localisation in the analysis step of a filter " \
                             "and computes a weight according to the specified distance " \
                             "and the settings for the localising function. " \
                             "Typically the function is called in `py__prodRinvA_l_pdaf` " \
                             "in the domain-localised filters. " \
                             "Also, the function is typically called for the LEnKF " \
                             "in the `py__localize_covar_pdaf`. \n    " \
                             "This function is usually only used in user-codes that do not use PDAF-OMI."

docstrings['print_info'] = "This function prints the wallclock time and memory measured by PDAF. This is called at the end of the DA program. \n    " \
                           "The function displays the following information: \n    " \
                           "- Memory required for the ensemble array, state vector, and transform matrix" \
                           "- Memory required by the analysis step" \
                           "- Memory required to perform the ensemble transformation"

docstrings['reset_forget'] = "This function allows a user to reset the forgetting factor manually during the assimilation process. " \
                             "For the local ensemble Kalman filters the forgetting factor can be set either globally of differently " \
                             "for each local analysis domain. " \
                             "For the LNETF and the global filters only a global setting of the forgeting factor is possible. " \
                             "In addition, the implementation of adaptive choices for the forgetting factor (beyond what is implemented in PDAF) are possible."

docstrings['SampleEns'] = "This function generates an ensemble from singular values and their vectors (EOF modes) centred on given mean state. " \
                          "The singular values and vectors are derived from the ensemble anomalies which can be obtained from a long model trajectory using " \
                          "`pyPDAF.PDAF.eofcovar`."

docstrings['set_debug_flag'] = "This function activates the debug output of the PDAF. Starting from the use of this function, the debug infomation is " \
                               "sent to screen output.  The screen output end when the debug flag is set to 0. " \
                               "We recommend using debugging output for single local domain, e.g. `if domain_p = 1: pyPDAF.PDAF.set_debug_flag(1)`.\n    "

docstrings['set_ens_pointer'] = "This function returns the ensemble array in a numpy array where the internal array data has the same memoery address as PDAF ensemble array."

docstrings['set_smootherens'] = "This function can be used in the offline implementation when a smoother is used. " \
                                "It is typically called in `py__init_ens_pdaf` in the call to `pyPDAF.PDAF.PDAF_init`. " \
                                "The function `pyPDAF.PDAF.set_smootherens` is used when the smoother extension of a filter is used. " \
                                "In this case, the smoothed ensemble states at earlier times are stored in an internal array of PDAF. " \
                                "To be able to smooth post times, the smoother algorithm must have access to the past ensembles. " \
                                "In the offline mode the user has to manually fill the smoother ensemble array from ensembles read in from files. " \
                                "In the online mode, the smoother array is filled automatically during the cycles of forecast phases and analysis steps. "

docstrings['seik_TtimesA'] = "This is an internal function in PDAF where it perform matrix calculation of B = TA. This allows for two types of T matrix. " \
                             "The resulting matrix B is the transformation matrix act on the full forecast ensemble. " \
                             "Mathematical description of the function is the second term of Eq. (23) and the T matrix is defined in Eq. (13) in\n    " \
                             "Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). A unification of ensemble square root Kalman filters. Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    "
docstrings['etkf_Tleft'] = "This is an internal function in PDAF where it perform matrix calculation of B = TA. " \
                           "This function performs the second term of Eq. (34) i\n    " \
                           "Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). A unification of ensemble square root Kalman filters. Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    "
docstrings['estkf_OmegaA'] = "This function is an internal function in PDAF. This function performs the second term of Eq. (29) i\n    " \
                           "Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). A unification of ensemble square root Kalman filters. Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    "
docstrings['enkf_omega'] = "Generation of a random matrix with orthogonal basis following SEEK approach for EnKF with given properties."
docstrings['seik_omega'] = "Generation of a random matrix with orthogonal basis following SEIK approach."
docstrings['incremental'] = "This is a helper function to apply analysis increment to model state in model forecast phase. It simply calls the user-supplied function. "
docstrings['add_increment'] = "This function directly adds analysis increment to given state vector without the need for user-supplied functions."
docstrings['local_weights'] = "This function returns a vector of the localisation weights based on distance and localisation functions and radii. " \
                              "This function is particularly useful for mannually apply covariance localisations for state or observation errors."
docstrings['diag_CRPS'] = "Obtain a continuous rank probability score for an ensemble. The implementation is based on " \
                          "This function follows follows Hersbach, H., 2000: Decomposition of the Continuous Ranked Probability Score for Ensemble Prediction Systems. Wea. Forecasting, 15, 559–570, https://doi.org/10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2\n    "
docstrings['force_analysis'] = "This function overwrite member index of the ensemble state by local_dim_ens (number of ensembles for current process, in full parallel setup, this is 1.) and the counter cnt_steps by nsteps-1.\n    " \
                               "This forces that the analysis step is executed at the next call to PDAF assimilation functions."
docstrings['gather_obs_f2_flex'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF)\n    " \
                             "this function returns the full observation coordinates " \
                             "from process-local observation coordinates. `pyPDAF.PDAF.gather_obs_f_flex` is used to get corresponding observations. " \
                             "Unlike `pyPDAF.PDAF.gather_obs_f2`, the function does not use depends on\n    " \
                             "`pyPDAF.PDAF.gather_dim_obs_f`"
docstrings['gather_obs_f_flex'] = "In the local filters (LESKTF, LETKF, LSEIK, LNETF) " \
                             "this function returns the total observation vector " \
                             "from process-local observations. `pyPDAF.PDAF.gather_obs_f2_flex` is used to get corresponding coordinates.\n    " \
                             "Unlike `pyPDAF.PDAF.gather_obs_f`, the function does not use depends on " \
                             "`pyPDAF.PDAF.gather_dim_obs_f`"
docstrings['prepost'] = "This function does not perform any DA. " \
                        "It is used to perform a preprocess and postprocess of the ensemble. " \
                        "Compared to `pyPDAF.PDAF.assimilate_prepost`, this function does not set assimilation flag.\n    " \
                        "The function is a combination of `pyPDAF.PDAF.put_state_prepost` " \
                        "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                        "in the following sequence: \n    " \
                        "1. py__collect_state_pdaf\n    " \
                        "2. py__prepoststep_state_pdaf\n    " \
                        "3. py__prepoststep_state_pdaf\n    " \
                        "4. py__distribute_state_pdaf\n    " \
                        "5. py__next_observation_pdaf\n    "
docstrings['set_memberid'] = "This function sets the ensemble member index to given value."
docstrings['set_comm_pdaf'] = "This function sets the MPI communicator of PDAF. This is by default `MPI_COMM_WORLD`. "
docstrings['set_offline_mode'] = "This function activates offline mode."
docstrings['print_domain_stats'] = "This function make screen output of statistics of the local domains on current process."
docstrings['init_local_obsstats'] = "This function initialise the observation statistics of local domain. " \
                                    "This statistics can be updated by `pyPDAF.PDAF.incr_local_obsstats`, " \
                                    "and can be viewed by `pyPDAF.PDAF.print_local_obsstats`."
docstrings['incr_local_obsstats'] = "This function update the observation statistics of local domain. " \
                                    "This statistics should be initialised by `pyPDAF.PDAF.init_local_obsstats`, " \
                                    "and can be viewed by `pyPDAF.PDAF.print_local_obsstats`."
docstrings['print_local_obsstats'] = "This function print the observation statistics of local domain on screen. " \
                                    "This statistics should be initialised by `pyPDAF.PDAF.init_local_obsstats`, " \
                                    "and can be updated by `pyPDAF.PDAF.incr_local_obsstats`."

docstrings['omit_obs_omi'] = "This function computes innovation and omit corresponding observations in assimilation if the innovation is too large. " \
                             "This function is used by some of the global filters, e.g. EnKF, LEnKF, PF, NETF, with OMI."
docstrings['diag_CRPS_nompi'] = "Obtain a continuous rank probability score for an ensemble without using MPI parallelisation. The implementation is based on " \
                          "This function follows follows Hersbach, H., 2000: Decomposition of the Continuous Ranked Probability Score for Ensemble Prediction Systems. Wea. Forecasting, 15, 559–570, https://doi.org/10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2\n    "


docstrings['omi_init'] = "This function initialise the number of observation types in OMI by allocating an array of `obs_f` derived types instances. This should be called before any other OMI functions."
docstrings['omi_set_doassim'] = "This function sets the `doassim` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                "If `doassim` is set to 0, the observation is not assimilated in the DA system. " \
                                "See https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdoassim"
docstrings['omi_set_disttype'] = "This function sets the `disttype` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                 "`disttype` determines the way the distance between observation and model grid is calculated in OMI. " \
                                 "See https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype"
docstrings['omi_set_ncoord'] = "This function sets the `ncoord` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                               "This is the dimension of coordinates of the observation. "
docstrings['omi_set_id_obs_p'] = "This function sets the `id_obs_p` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                 "`id_obs_p(nrows, dim_obs_p)` is a 2D array of integers.\n    " \
                                 "The value of `nrows` depends on the observation operator used for an observation. " \
                                 "Examples: \n     `nrows=1` if observations are located on model grid point. In this case,\n    " \
                                 "`id_obs_p` stores the index of the state vector (starting from 1) corresponds to the observations,\n    " \
                                 "e.g. `id_obs_p[0, j] = i means that the `j`-th observation is located on the same grid point and is the " \
                                 "same physical quantity as the `i`-th element of the state vector. \n    " \
                                 "If `nrows=4`, each observations corresponds to 4 indices of the state vector. In this case,\n    " \
                                 "the location of these state vector is used to perform bi-linear interpolation from model grid to observation location. " \
                                 "This information is used in the `pyPDAF.PDAF.omi_obs_op_gridavg` and `pyPDAF.PDAF.omi_obs_op_interp_lin` functions. " \
                                 "When interpolation is needed, the weighting of the interpolation is done in the `pyPDAF.PDAF.omi_get_interp_coeff_lin`,  " \
                                 "`pyPDAF.PDAF.omi_get_interp_coeff_lin1D`, and `pyPDAF.PDAF.omi_get_interp_coeff_tri` functions.\n    "
docstrings['omi_set_icoeff_p'] = "This function sets the `icoeff_p` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                 "`icoeff_p(nrows, dim_obs_p)` is a 2D array of real number used to implement\n    " \
                                 "interpolations. This is used in tandem with `id_obs_p`. " \
                                 "Checking the documentation of `pyPDAF.PDAF.omi_set_id_obs_p` for some details. " \
                                 "Also, see https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients for setting these values."
docstrings['omi_set_domainsize'] = "This function sets the `domainsize` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                   "`domainsize(ncoord)` is the size of the domain in each spatial dimension. " \
                                   "This information is used to compute the Cartesian disance with periodic boundary. " \
                                   "If the value of one dimension is `<=0`, no periodicity is assumed in that dimension. "
docstrings['omi_set_obs_err_type'] = "This function sets the `obs_err_type` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                     "`obs_err_type` is an integer that specifies the type of observation error. "
docstrings['omi_set_use_global_obs'] = "This function sets the `use_global_obs` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                       "In the domain-localized filters (LESTK, LETKF, LSEIK, LNETF) " \
                                       "observations are assimilated that are located within the localization around some grid point. " \
                                       "When a model uses parallelisation with domain-decomposition some of these observations might belong to a different process-domain. " \
                                       "In the default mode (use_global_obs=1) PDAF-OMI gathers all globally available observations so that each process has access to all observations.\n    " \
                                       "It can be more efficient to limit the observations on a process-domain to those observations that are located inside the domain or within the localization radius around it. " \
                                       "Then, in the local analyses less observations have to be checked for their distance. " \
                                       "Setting use_global_obs=0 activates this feature. However, it needs additional preparations to make PDAF-OMI aware of the limiting coordinates of a process sub-domain. " \
                                       "See https://pdaf.awi.de/trac/wiki/OMI_use_global_obs for the use of `pyPDAF.PDAF.omi_set_domain_limits`."
docstrings['omi_set_inno_omit'] = "This function sets the `inno_omit` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                  "Setting this variable to a value > 0.0 activates the functionality that observations are omitted (made irrelevant) from the analysis update " \
                                  "if the difference of their value and the ensemble mean to too large. " \
                                  "If inno_omit=2.0, an observation would be omitted if the squared difference between the observed ensemble mean state and the observation value is larger than 2 times the observation error variance\n    " \
                                  "See https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#Omittingobservationsthatarepotentialoutliers"
docstrings['omi_set_inno_omit_ivar'] = "This function sets the `inno_omit_ivar` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                                       "This is used to specify the inverse of the observations variance to omit the observation. " \
                                       "By default it is `1e-12` for a large observation error, but users can adjust this value to ensure that the observation is omitted based on applications\n    "
docstrings['omi_gather_obs'] = "This function is typically called in the user-supplied function `py__init_dim_obs_pdaf`. " \
                               "This function returns the full observation dimensioin from process-local observations. " \
                               "It also sets the observation vector, its coordinates, and the inverse of the observation variance. " \
                               "This function furtuer sets the localisation radius in OMI."
docstrings['omi_gather_obsstate'] = "This function is used to implement custom observation operators. " \
                                    "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Implementingyourownobservationoperator"
docstrings['omi_set_domain_limits'] = "This is used to set the domain limits for the use of `pyPDAF.PDAF.omi_set_use_global_obs`." \
                                      "Currently, it only supports 2D limitations. See https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_set_domain_limit\n    "
docstrings['omi_set_debug_flag'] = "This sets the debug flag for OMI. If set to 1, debug information is printed to the screen.\n    " \
                                   "The debug flag can be set to 0 to stop the debugging. See https://pdaf.awi.de/trac/wiki/OMI_debugging"
docstrings['omi_deallocate_obs'] = "It deallocates teh OMI-internal obsrevation arrays but this should not be called as it is called internally in PDAF."
docstrings['omi_obs_op_gridpoint'] = "The routine provides an observation operator for observations that are the same grid point values in the state vector. " \
                                     "The function is used in the user-supplied function `py__obs_op_pdaf`. " \
                                     "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAF-OMIObservationOperators"
docstrings['omi_obs_op_gridavg'] = "The routine provides an observation operator for observations that are the average of some grid point values in the state vector. " \
                                   "The function is used in the user-supplied function `py__obs_op_pdaf`. " \
                                   "This function is used with `id_obs_p` property of `obs_f`" \
                                   "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAF-OMIObservationOperators"
docstrings['omi_obs_op_interp_lin'] = "The routine provides an observation operator for observations that are the average of some grid point values in the state vector. " \
                                      "The function is used in the user-supplied function `py__obs_op_pdaf`. " \
                                      "This function is used with `id_obs_p` and `icoeff_p` property of `obs_f`" \
                                      "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAF-OMIObservationOperators"
docstrings['omi_obs_op_adj_gridavg'] = "This is the adjoint of `pyPDAF.PDAF.omi_obs_op_gridavg`." \
                                       "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Adjointobservationoperators"
docstrings['omi_obs_op_adj_gridpoint'] = "This is the adjoint of `pyPDAF.PDAF.omi_obs_op_gridpoint`." \
                                         "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Adjointobservationoperators"
docstrings['omi_obs_op_adj_interp_lin'] = "This is the adjoint of `pyPDAF.PDAF.omi_obs_op_interp_lin`." \
                                         "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Adjointobservationoperators"
docstrings['omi_get_interp_coeff_tri'] = "This function returns the coefficient for linear interpolation in 2D on unstructure triangular grid.\n    "
docstrings['omi_get_interp_coeff_lin1D'] = "This function returns the coefficient for linear interpolation in 1D.\n    "
docstrings['omi_get_interp_coeff_lin'] = "This function returns the coefficient for linear interpolation up to 3D.\n    " \
                                         "See https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin"

docstrings['omi_assimilate_3dvar'] = "Using 3DVar for DA with diagonal observation error covariance matrix.\n    " \
                                     "This is a deterministic filtering scheme. " \
                                     "This function should be called at each model time step. \n    \n    " \
                                     "The function is a combination of `pyPDAF.PDAF.omi_put_state_3dvar`\n    " \
                                     "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                     "in the following sequence: \n    " \
                                     "1. py__collect_state_pdaf\n    " \
                                     "2. py__prepoststep_state_pdaf\n    " \
                                     "3. py__init_dim_obs_pdaf\n    " \
                                     "4. py__obs_op_pdaf\n    " \
                                     "Starting the iterative optimisation:\n    " \
                                     "5. py__cvt_pdaf\n    " \
                                     "6. py__obs_op_lin_pdaf\n    " \
                                     "7. py__obs_op_adj_pdaf\n    " \
                                     "8. py__cvt_adj_pdaf\n    " \
                                     "9. core DA algorithm\n    " \
                                     "After the iterations: \n    " \
                                     "10. py__cvt_pdaf\n    " \
                                     "11. py__prepoststep_state_pdaf\n    " \
                                     "12. py__distribute_state_pdaf\n    " \
                                     "13. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_en3dvar_estkf'] =  "Using 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                              "The background error covariance matrix is estimated by ensemble. " \
                                              "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                              "An ESTKF is used to generate ensemble perturbations. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf\n    " \
                                              "Starting the iterative optimisation:\n    " \
                                              "5. py__cvt_ens_pdaf\n    " \
                                              "6. py__obs_op_lin_pdaf\n    " \
                                              "7. py__obs_op_adj_pdaf\n    " \
                                              "8. py__cvt_adj_ens_pdaf\n    " \
                                              "9. core 3DEnVar algorithm\n    " \
                                              "After the iterations: \n    " \
                                              "10. py__cvt_ens_pdaf\n    " \
                                              "Perform ESTKF: " \
                                              "11. py__init_dim_obs_pdaf\n    " \
                                              "12. py__obs_op_pdaf (for ensemble mean\n    " \
                                              "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "14. core ESTKF algorithm\n    " \
                                              "15. py__prepoststep_state_pdaf\n    " \
                                              "16. py__distribute_state_pdaf\n    " \
                                              "17. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_en3dvar_lestkf'] = "It is recommended to use `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf` for better efficiency. \n    \n   \n    " \
                                              "Using 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                              "The background error covariance matrix is estimated by ensemble. " \
                                              "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                              "An LESTKF is used to generate ensemble perturbations. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_en3dvar_lestkf`\n    " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf\n    " \
                                              "Starting the iterative optimisation:\n    " \
                                              "5. py__cvt_ens_pdaf\n    " \
                                              "6. py__obs_op_lin_pdaf\n    " \
                                              "7. py__obs_op_adj_pdaf\n    " \
                                              "8. py__cvt_adj_ens_pdaf\n    " \
                                              "9. core DA algorithm\n    " \
                                              "After the iterations: \n    " \
                                              "10. py__cvt_ens_pdaf\n    " \
                                              "Perform LESTKF: \n    " \
                                              "11. py__init_n_domains_p_pdaf\n    " \
                                              "12. py__init_dim_obs_pdaf\n    " \
                                              "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "loop over each local domain:\n    " \
                                              "14. py__init_dim_l_pdaf\n    " \
                                              "15. py__init_dim_obs_l_pdaf\n    " \
                                              "16. py__g2l_state_pdaf (localise mean ensemble in observation space)\n    " \
                                              "17. core DA algorithm\n    " \
                                              "18. py__l2g_state_pdaf\n    " \
                                              "19. py__prepoststep_state_pdaf\n    " \
                                              "20. py__distribute_state_pdaf\n    " \
                                              "21. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_global'] = "Using global filters for DA except for 3DVars with diagonal observation error covariance matrix.\n    " \
                                      "The function is a combination of `pyPDAF.PDAF.put_state_global` " \
                                      "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                      "This function should be called at each model time step. \n    \n    " \
                                      "in the following sequence: \n    " \
                                      "1. py__collect_state_pdaf\n    " \
                                      "2. py__prepoststep_state_pdaf\n    " \
                                      "3. py__init_dim_obs_pdaf\n    " \
                                      "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                      "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                      "6. core DA algorithm\n    " \
                                      "7. py__prepoststep_state_pdaf\n    " \
                                      "8. py__distribute_state_pdaf\n    " \
                                      "9. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_hyb3dvar_estkf'] = "Using Hybrid 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                              "Here, the background error covariance is hybridised by a static background error covariance, " \
                                              "and a flow-dependent background error covariance estimated from ensemble. " \
                                              "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                              "ESTKF in this implementation. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf\n    " \
                                              "Starting the iterative optimisation:\n    " \
                                              "5. py__cvt_pdaf\n    " \
                                              "6. py__cvt_ens_pdaf\n    " \
                                              "7. py__obs_op_lin_pdaf\n    " \
                                              "8. py__obs_op_adj_pdaf\n    " \
                                              "9. py__cvt_adj_pdaf\n    " \
                                              "10. py__cvt_adj_ens_pdaf\n    " \
                                              "11. core 3DEnVar algorithm\n    " \
                                              "After the iterations: \n    " \
                                              "12. py__cvt_pdaf\n    " \
                                              "13. py__cvt_ens_pdaf\n    " \
                                              "Perform ESTKF: " \
                                              "14. py__init_dim_obs_pdaf\n    " \
                                              "15. py__obs_op_pdaf (for ensemble mean\n    " \
                                              "16. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "17. core ESTKF algorithm\n    " \
                                              "18. py__prepoststep_state_pdaf\n    " \
                                              "19. py__distribute_state_pdaf\n    " \
                                              "20. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_hyb3dvar_lestkf'] = "It is recommended to use `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf` for better efficiency. \n    \n   \n    " \
                                               "Using Hybrid 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                               "Here, the background error covariance is hybridised by a static background error covariance, " \
                                               "and a flow-dependent background error covariance estimated from ensemble. " \
                                               "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                               "LESTKF in this implementation. " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "The function is a combination of `pyPDAF.PDAF.omi_put_state_hyb3dvar_lestkf`\n    " \
                                               "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_dim_obs_pdaf\n    " \
                                               "4. py__obs_op_pdaf\n    " \
                                               "Starting the iterative optimisation:\n    " \
                                               "5. py__cvt_pdaf\n    " \
                                               "6. py__cvt_ens_pdaf\n    " \
                                               "7. py__obs_op_lin_pdaf\n    " \
                                               "8. py__obs_op_adj_pdaf\n    " \
                                               "9. py__cvt_adj_pdaf\n    " \
                                               "10. py__cvt_adj_ens_pdaf\n    " \
                                               "11. core DA algorithm\n    " \
                                               "After the iterations: \n    " \
                                               "12. py__cvt_pdaf\n    " \
                                               "13. py__cvt_ens_pdaf\n    " \
                                               "Perform LESTKF: \n    " \
                                               "14. py__init_n_domains_p_pdaf\n    " \
                                               "15. py__init_dim_obs_pdaf\n    " \
                                               "16. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "loop over each local domain:\n    " \
                                               "17. py__init_dim_l_pdaf\n    " \
                                               "18. py__init_dim_obs_l_pdaf\n    " \
                                               "19. py__g2l_state_pdaf\n    " \
                                               "20. py__init_obs_l_pdaf\n    "\
                                               "21. core DA algorithm\n    " \
                                               "22. py__l2g_state_pdaf\n    " \
                                               "23. py__prepoststep_state_pdaf\n    " \
                                               "24. py__distribute_state_pdaf\n    " \
                                               "25. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_lenkf'] = "Using stochastic EnKF (ensemble " \
                                     "Kalman filter) with covariance localisation " \
                                     "for DA  with diagonal observation error covariance matrix. " \
                                     "This is the only scheme for covariance localisation in PDAF. " \
                                     "This function should be called at each model time step. \n    \n    " \
                                     "The function is a combination of `pyPDAF.PDAF.omi_put_state_lenkf` " \
                                     "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                     "in the following sequence: \n    " \
                                     "1. py__collect_state_pdaf\n    " \
                                     "2. py__prepoststep_state_pdaf\n    " \
                                     "3. py__init_dim_obs_pdaf\n    " \
                                     "4. py__obs_op_pdaf (for each ensemble member\n    " \
                                     "5. py__localize_pdaf\n    " \
                                     "6. py__obs_op_pdaf (repeated to reduce storage\n    " \
                                     "7. core DA algorith\n    " \
                                     "8. py__prepoststep_state_pdaf\n    " \
                                     "9. py__distribute_state_pdaf\n    " \
                                     "10. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_local'] = "It is recommended to use `pyPDAF.PDAF.localomi_assimilate` for better efficiency." \
                                     "Using domain localised filters for DA with diagonal observation error covariance matrix. " \
                                     "This function should be called at each model time step. \n    \n    " \
                                     "The function is a combination of `pyPDAF.PDAF.omi_put_state_local` " \
                                     "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                     "in the following sequence: \n    " \
                                     "1. py__collect_state_pdaf\n    " \
                                     "2. py__prepoststep_state_pdaf\n    " \
                                     "3. py__init_n_domains_p_pdaf\n    " \
                                     "4. py__init_dim_obs_pdaf\n    " \
                                     "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                     "loop over each local domain:\n    " \
                                     "6. py__init_dim_l_pdaf\n    " \
                                     "7. py__init_dim_obs_l_pdaf\n    " \
                                     "8. py__g2l_state_pdaf\n    " \
                                     "9. py__init_obs_l_pdaf\n    "\
                                     "10. core DA algorithm\n    " \
                                     "11. py__l2g_state_pdaf\n    " \
                                     "12. py__prepoststep_state_pdaf\n    " \
                                     "13. py__distribute_state_pdaf\n    " \
                                     "14. py__next_observation_pdaf\n    "
docstrings['omi_generate_obs'] = "This function generates synthetic observations based on each member of model forecast with diagonal observation error covariance matrix. " \
                                 "This is based on the usual implementation strategy for PDAF.\n    \n    " \
                                 "The function calls `pyPDAF.PDAF.generate_obs` " \
                                 " and executes the user-supplied function " \
                                 "in the following sequence: \n    " \
                                 "1. py__collect_state_pdaf\n    " \
                                 "2. py__prepoststep_state_pdaf\n    " \
                                 "3. py__init_dim_obs_pdaf\n    " \
                                 "4. py__obs_op_pda\n    " \
                                 "5. py__get_obs_f_pdaf\n    " \
                                 "6. py__prepoststep_state_pdaf\n    " \
                                 "7. py__distribute_state_pdaf\n    " \
                                 "8. py__next_observation_pdaf\n    "
docstrings['omi_put_state_3dvar'] = "Using 3DVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                    "This is a deterministic filtering scheme. " \
                                    "This function is usually used in 'flexible' parallelisation, " \
                                    "but 3dvar is deterministic and does not require ensemble.\n    " \
                                    "A `pyPDAF.PDAF.get_state` function should be used to post-process the state vector and " \
                                    "distribute the state vector " \
                                    "to the model after this function. " \
                                    "This function should be called at each model time step. \n    \n    " \
                                    "The function executes the user-supplied function " \
                                    "in the following sequence: \n    " \
                                    "1. py__collect_state_pdaf\n    " \
                                    "2. py__prepoststep_state_pdaf\n    " \
                                    "3. py__init_dim_obs_pdaf\n    " \
                                    "4. py__obs_op_pdaf\n    " \
                                    "Starting the iterative optimisation:\n    " \
                                    "5. py__cvt_pdaf\n    " \
                                    "6. py__obs_op_lin_pdaf\n    " \
                                    "7. py__obs_op_adj_pdaf\n    " \
                                    "8. py__cvt_adj_pdaf\n    " \
                                    "9. core DA algorithm\n    " \
                                    "After the iterations: \n    " \
                                    "10. py__cvt_pdaf\n    "
docstrings['omi_put_state_en3dvar_estkf'] = "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                            "The background error covariance matrix is estimated by ensemble. " \
                                            "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                            "An ESTKF is used to generate ensemble perturbations. " \
                                            "This function is usually used in 'flexible' parallelisation. " \
                                            "i.e., the ensemble size is larger than the available number of processes. " \
                                            "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                            "to the model after this function. " \
                                            "This function should be called at each model time step. \n    \n    " \
                                            "The function executes the user-supplied function " \
                                            "in the following sequence: \n    " \
                                            "1. py__collect_state_pdaf\n    " \
                                            "2. py__prepoststep_state_pdaf\n    " \
                                            "3. py__init_dim_obs_pdaf\n    " \
                                            "4. py__obs_op_pdaf\n    " \
                                            "Starting the iterative optimisation:\n    " \
                                            "5. py__cvt_ens_pdaf\n    " \
                                            "6. py__obs_op_lin_pdaf\n    " \
                                            "7. py__obs_op_adj_pdaf\n    " \
                                            "8. py__cvt_adj_ens_pdaf\n    " \
                                            "9. core 3DEnVar algorithm\n    " \
                                            "After the iterations: \n    " \
                                            "10. py__cvt_ens_pdaf\n    " \
                                            "Perform ESTKF: " \
                                            "11. py__init_dim_obs_pdaf\n    " \
                                            "12. py__obs_op_pdaf (for ensemble mean\n    " \
                                            "13. py__obs_op_pdaf (for each ensemble member (only relevant for adaptive forgetting factor schemes)\n    " \
                                            "14. core ESTKF algorithm\n    "
docstrings['omi_put_state_en3dvar_lestkf'] = "It is recommended to use `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf` for better efficiency. \n    \n   \n    " \
                                             "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                             "The background error covariance matrix is estimated by ensemble. " \
                                             "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                             "An LESTKF is used to generate ensemble perturbations. " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_dim_obs_pdaf\n    " \
                                             "4. py__obs_op_pdaf\n    " \
                                             "Starting the iterative optimisation:\n    " \
                                             "5. py__cvt_ens_pdaf\n    " \
                                             "6. py__obs_op_lin_pdaf\n    " \
                                             "7. py__obs_op_adj_pdaf\n    " \
                                             "8. py__cvt_adj_ens_pdaf\n    " \
                                             "9. core DA algorithm\n    " \
                                             "After the iterations: \n    " \
                                             "10. py__cvt_ens_pdaf\n    " \
                                             "Perform LESTKF: \n    " \
                                             "11. py__init_n_domains_p_pdaf\n    " \
                                             "12. py__init_dim_obs_pdaf\n    " \
                                             "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "loop over each local domain:\n    " \
                                             "14. py__init_dim_l_pdaf\n    " \
                                             "15. py__init_dim_obs_l_pdaf\n    " \
                                             "16. py__g2l_state_pdaf\n    " \
                                             "17. core DA algorithm\n    " \
                                             "18. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_generate_obs'] = "This function generates synthetic observations based on each member of model forecast with diagonal observation error covariance matrix. " \
                                           "This function is for the case where the ensemble size is larger than the number of processors. \n    \n    " \
                                           "The function executes the user-supplied function " \
                                           "in the following sequence: \n    " \
                                           "1. py__collect_state_pdaf\n    " \
                                           "2. py__prepoststep_state_pdaf\n    " \
                                           "3. py__init_dim_obs_pdaf\n    " \
                                           "4. py__obs_op_pda\n    " \
                                           "5. py__get_obs_f_pdaf\n    "
docstrings['omi_put_state_global'] = "Using global filters for DA except for 3DVars without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                     "This function should be called at each model time step. \n    \n    " \
                                     "in the following sequence: \n    " \
                                     "1. py__collect_state_pdaf\n    " \
                                     "2. py__prepoststep_state_pdaf\n    " \
                                     "3. py__init_dim_obs_pdaf\n    " \
                                     "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                     "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                     "6. core DA algorithm\n    "
docstrings['omi_put_state_hyb3dvar_estkf'] =  "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                              "Here, the background error covariance is hybridised by a static background error covariance, " \
                                              "and a flow-dependent background error covariance estimated from ensemble. " \
                                              "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                              "ESTKF in this implementation. \n    \n    " \
                                              "This function is usually used in 'flexible' parallelisation. " \
                                              "i.e., the ensemble size is larger than the available number of processes. " \
                                              "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                              "to the model after this function. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf\n    " \
                                              "Starting the iterative optimisation:\n    " \
                                              "5. py__cvt_ens_pdaf\n    " \
                                              "6. py__obs_op_lin_pdaf\n    " \
                                              "7. py__obs_op_adj_pdaf\n    " \
                                              "8. py__cvt_adj_ens_pdaf\n    " \
                                              "9. core 3DEnVar algorithm\n    " \
                                              "After the iterations: \n    " \
                                              "10. py__cvt_ens_pdaf\n    " \
                                              "Perform ESTKF: " \
                                              "11. py__init_dim_obs_pdaf\n    " \
                                              "12. py__obs_op_pdaf (for ensemble mean\n    " \
                                              "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "14. core ESTKF algorithm\n    "
docstrings['omi_put_state_hyb3dvar_lestkf'] =  "It is recommended to use `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf` for better efficiency. \n    \n   \n    " \
                                               "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                               "Here, the background error covariance is hybridised by a static background error covariance, " \
                                               "and a flow-dependent background error covariance estimated from ensemble. " \
                                               "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                               "LESTKF in this implementation. \n    \n    " \
                                               "This function is usually used in 'flexible' parallelisation. " \
                                               "i.e., the ensemble size is larger than the available number of processes. " \
                                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                               "to the model after this function. " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "The function executes the user-supplied function " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_dim_obs_pdaf\n    " \
                                               "4. py__obs_op_pdaf\n    " \
                                               "Starting the iterative optimisation:\n    " \
                                               "5. py__cvt_ens_pdaf\n    " \
                                               "6. py__obs_op_lin_pdaf\n    " \
                                               "7. py__obs_op_adj_pdaf\n    " \
                                               "8. py__cvt_adj_ens_pdaf\n    " \
                                               "9. core DA algorithm\n    " \
                                               "After the iterations: \n    " \
                                               "10. py__cvt_ens_pdaf\n    " \
                                               "Perform LESTKF: \n    " \
                                               "11. py__init_n_domains_p_pdaf\n    " \
                                               "12. py__init_dim_obs_pdaf\n    " \
                                               "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "loop over each local domain:\n    " \
                                               "14. py__init_dim_l_pdaf\n    " \
                                               "15. py__init_dim_obs_l_pdaf\n    " \
                                               "16. py__g2l_state_pdaf (localise mean ensemble in observation space)\n    " \
                                               "17. core DA algorithm\n    " \
                                               "18. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_lenkf'] = "Using stochastic EnKF (ensemble " \
                                    "Kalman filter) with covariance localisation " \
                                    "for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix. " \
                                    "This is the only scheme for covariance localisation in PDAF. " \
                                    "This is the only scheme for covariance localisation in PDAF. " \
                                    "This function is usually used in 'flexible' parallelisation. " \
                                    "i.e., the ensemble size is larger than the available number of processes. " \
                                    "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                    "to the model after this function. " \
                                    "This function should be called at each model time step. \n    \n    " \
                                    "The function executes the user-supplied function " \
                                    "in the following sequence: \n    " \
                                    "1. py__collect_state_pdaf\n    " \
                                    "2. py__prepoststep_state_pdaf\n    " \
                                    "3. py__init_dim_obs_pdaf\n    " \
                                    "4. py__obs_op_pdaf (for each ensemble member\n    " \
                                    "5. py__localize_pdaf\n    " \
                                    "6. py__obs_op_pdaf (repeated to reduce storage\n    " \
                                    "7. core DA algorith\n    "
docstrings['omi_put_state_local'] = "It is recommended to use `pyPDAF.PDAF.localomi_put_state` for better efficiency. \n    \n    " \
                                    "Using domain localised filters for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix. " \
                                    "This is a domain localisation method. " \
                                    "This function is usually used in 'flexible' parallelisation. " \
                                    "i.e., the ensemble size is larger than the available number of processes. " \
                                    "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                    "to the model after this function. " \
                                    "This function should be called at each model time step. " \
                                    "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                    "The function executes the user-supplied function " \
                                    "in the following sequence: \n    " \
                                    "1. py__collect_state_pdaf\n    " \
                                    "2. py__prepoststep_state_pdaf\n    " \
                                    "3. py__init_n_domains_p_pdaf\n    " \
                                    "4. py__init_dim_obs_pdaf\n    " \
                                    "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                    "loop over each local domain:\n    " \
                                    "6. py__init_dim_l_pdaf\n    " \
                                    "7. py__init_dim_obs_l_pdaf\n    " \
                                    "8. py__g2l_state_pdaf\n    " \
                                    "9. core DA algorithm\n    " \
                                    "10. py__l2g_state_pdaf\n    "

docstrings['omi_init_obs_f_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to initialise the observation vector. " \
                                  "This could be used to modify the observation vector when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_init_obsvar_cb'] = "This function is an internal PDAF function that is used as a call-back function to initialise the observation error variance. " \
                                   "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_g2l_obs_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to convert between global and local observation vectors in domain localisation.\n    "
docstrings['omi_init_obs_l_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to initialise local observation vector in domain localisation. " \
                                  "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_init_obsvar_l_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to initialise local observation vector in domain localisation. " \
                                     "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_prodRinvA_l_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to perform the matrix multiplication inverse of local observation error covariance and a matrix A in domain localisation. " \
                                   "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_likelihood_l_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to compute the likelihood of the observation for a given ensemble member according to the observations used for the local analysis in the localized LNETF. " \
                                 "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`. See https://pdaf.awi.de/trac/wiki/U_likelihood_l"
docstrings['omi_prodRinvA_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to perform the matrix multiplication inverse of observation errro covariance and a matrix A. " \
                                 "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_likelihood_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to compute the likelihood of the observation for a given ensemble member according to the observations used for the local analysis for NETF or particle filter. " \
                                 "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`. See https://pdaf.awi.de/trac/wiki/U_likelihood_l"
docstrings['omi_add_obs_error_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to add random observation error to stochastic EnKF. " \
                                     "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`. See https://pdaf.awi.de/trac/wiki/U_likelihood_l"
docstrings['omi_init_obscovar_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to construct a full observation error covariance matrix used only in stochastic EnKF. " \
                                     "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_init_obserr_f_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to construct a full observation error covariance matrix used only in stochastic EnKF. " \
                                     "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_prodRinvA_hyb_l_cb'] = "This function is an internal PDAF-OMI function that is used as a call-back function to perform the matrix multiplication inverse of local observation error covariance and a matrix A in LKNETF. " \
                                       "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."
docstrings['omi_likelihood_hyb_l_cb'] = "This is an internal PDAF-OMI function that is used as a call-back function to compute the likelihood of the observation for a given ensemble member according to the observations used for the local analysis in LKNETF. " \
                                        "This could be used to modify the observation variance when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`. See https://pdaf.awi.de/trac/wiki/U_likelihood_l"

docstrings['omi_obsstats_l'] = "This function is called in the update routine of local filters and write statistics on locally used and excluded observations."
docstrings['omi_weights_l'] = "This function computes a weight vector according to the distances of observations from the local analysis domain with a vector of localisation radius."
docstrings['omi_weights_l_sgnl'] = "This function computes a weight vector according to the distances of observations from the local analysis domain with given localisation radius."
docstrings['omi_check_error'] = "This function returns the value of the PDAF-OMI internal error flag."
docstrings['omi_gather_obsdims'] = "This function gathers the information about the full dimension of each observation type in each process-local subdomain."
docstrings['omi_obsstats'] = "The function is called in the update routine of global filters and writes statistics on used and excluded observations."
docstrings['omi_init_dim_obs_l_iso'] = "The function has to be called in `init_dim_obs_l_OBTYPE` in each observation module if a domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used. " \
                                       "It initialises the local observation information for PDAF-OMI for a single local analysis domain. This is used for isotropic localisation " \
                                       "where the localisation radius is the same in all directions."
docstrings['omi_init_dim_obs_l_noniso'] = "The function has to be called in `init_dim_obs_l_OBTYPE` in each observation module if a domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used. " \
                                          "It initialises the local observation information for PDAF-OMI for a single local analysis domain. This is used for non-isotropic localisation " \
                                          "where the localisation radius is different in each direction. See https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE and https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l#Settingsfornon-isotropiclocalization."
docstrings['omi_init_dim_obs_l_noniso_locweights'] = "The function has to be called in `init_dim_obs_l_OBTYPE` in each observation module if a domain-localized filter (LESTKF/LETKF/LNETF/LSEIK)is used. " \
                                          "It initialises the local observation information for PDAF-OMI for a single local analysis domain. This is used for non-isotropic localisation and different weight functions for horizontal and vertical directions. " \
                                          "See https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE and https://pdaf.awi.de/trac/wiki/PDAFomi_init_dim_obs_l#Settingdifferentweightfunctsforhorizontalandverticaldirections."
docstrings['omi_localize_covar_iso'] = "The function has to be called in `localize_covar_OBTYPE` in each observation module. It applies the covariance localisation in stochastic EnKF. This is used for isotropic localisation " \
                                       "where the localisation radius is the same in all directions. See https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar"
docstrings['omi_localize_covar_noniso'] = "The function has to be called in `localize_covar_OBTYPE` in each observation module. It applies the covariance localisation in stochastic EnKF. This is used for non-isotropic localisation " \
                                       "where the localisation radius is different. See https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar"
docstrings['omi_localize_covar_noniso_locweights'] = "The function has to be called in `localize_covar_OBTYPE` in each observation module. It applies the covariance localisation in stochastic EnKF. This is used for non-isotropic localisation with different weight function for horizontal and vertical directions. " \
                                       "where the localisation radius is different. See https://pdaf.awi.de/trac/wiki/PDAFomi_localize_covar"

docstrings['omi_omit_by_inno_l_cb'] = "The function is called during the analysis step on each local analysis domain. " \
                                      "It checks the size of the innovation and sets the observation error to a high value " \
                                      "if the squared innovation exceeds a limit relative to the observation error variance." \
                                      "This function is an internal PDAF-OMI function that is used as a call-back function. " \
                                      "This could be used to modify the observation vector when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."

docstrings['omi_omit_by_inno_cb'] = "The function is called during the analysis step of a global filter. " \
                                    "It checks the size of the innovation and " \
                                    "sets the observation error to a high value " \
                                    "if the squared innovation exceeds a limit relative to the observation error variance.""This function is called in the update routine of local filters and write statistics on locally used and excluded observations." \
                                    "This function is an internal PDAF-OMI function that is used as a call-back function. " \
                                    "This could be used to modify the observation vector when OMI is used with `pyPDAF.PDAF.assimilate_xxx` instead of `pyPDAF.PDAF.omi_assimilate_xxx`."

docstrings['omi_set_localization'] = "This function sets localization information " \
                                     "(locweight, cradius, sradius) in OMI, " \
                                     "and allocates local arrays for cradius and sradius, i.e. `obs_l`. " \
                                     "This variant is for isotropic localization. " \
                                     "The function is used by user-supplied implementations of " \
                                     "`pyPDAF.PDAF.omi_init_dim_obs_l_iso`. "
docstrings['omi_set_localization_noniso'] = "This function sets localization information " \
                                            "(locweight, cradius, sradius) in OMI, " \
                                            "and allocates local arrays for cradius and sradius, i.e. `obs_l`. " \
                                            "This variant is for non-isotropic localization. " \
                                            "The function is used by user-supplied implementations of " \
                                            "`pyPDAF.PDAF.omi_init_dim_obs_l_noniso`. "
docstrings['omi_set_dim_obs_l'] = "This function initialises number local observations. " \
                                  "It also returns number of local observations up to the current observation type. " \
                                  "It is used by a user-supplied implementations of " \
                                  "`pyPDAF.PDAF.omi_init_dim_obs_l_xxx`."
docstrings['omi_store_obs_l_index'] = "This function stores the mapping index " \
                                      "between the global and local observation vectors, " \
                                      "the distance and the cradius and sradius " \
                                      "for a single observations in OMI. " \
                                      "This variant is for non-factorised localisation. " \
                                      "The function is used by user-supplied implementations of " \
                                      "`pyPDAF.PDAF.omi_init_dim_obs_l_iso` or `pyPDAF.PDAF.omi_init_dim_obs_l_noniso`. "
docstrings['omi_store_obs_l_index_vdist'] = "This function stores the mapping index " \
                                            "between the global and local observation vectors, " \
                                            "the distance and the cradius and sradius " \
                                            "for a single observations in OMI. " \
                                            "This variant is for 2D+1D factorised localisation.\n    " \
                                            "The function is used by user-supplied implementations of " \
                                            "`pyPDAF.PDAF.omi_init_dim_obs_l_noniso_locweights`."

docstrings['omi_assimilate_3dvar_nondiagR'] =  "Using 3DVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                    "This is a deterministic filtering scheme. " \
                                    "This function should be called at each model time step. \n    \n    " \
                                    "The function is a combination of `pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`\n    " \
                                    "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                    "in the following sequence: \n    " \
                                    "1. py__collect_state_pdaf\n    " \
                                    "2. py__prepoststep_state_pdaf\n    " \
                                    "3. py__init_dim_obs_pdaf\n    " \
                                    "4. py__obs_op_pdaf\n    " \
                                    "Starting the iterative optimisation:\n    " \
                                    "5. py__cvt_pdaf\n    " \
                                    "6. py__obs_op_lin_pdaf\n    " \
                                    "7. py__prodRinvA_pdaf\n    " \
                                    "8. py__obs_op_adj_pdaf\n    " \
                                    "9. py__cvt_adj_pdaf\n    " \
                                    "10. core DA algorithm\n    " \
                                    "After the iterations: \n    " \
                                    "11. py__cvt_pdaf\n    " \
                                    "12. py__prepoststep_state_pdaf\n    " \
                                    "13. py__distribute_state_pdaf\n    " \
                                    "14. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_en3dvar_estkf_nondiagR'] = "Using 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                                      "The background error covariance matrix is estimated by ensemble. " \
                                                      "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                      "An ESTKF is used to generate ensemble perturbations. " \
                                                      "This function should be called at each model time step. \n    \n    " \
                                                      "The function is a combination of `pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`\n    " \
                                                      "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                      "in the following sequence: \n    " \
                                                      "1. py__collect_state_pdaf\n    " \
                                                      "2. py__prepoststep_state_pdaf\n    " \
                                                      "3. py__init_dim_obs_pdaf\n    " \
                                                      "4. py__obs_op_pdaf\n    " \
                                                      "Starting the iterative optimisation:\n    " \
                                                      "5. py__cvt_ens_pdaf\n    " \
                                                      "6. py__obs_op_lin_pdaf\n    " \
                                                      "7. py__prodRinvA_pdaf\n    " \
                                                      "8. py__obs_op_adj_pdaf\n    " \
                                                      "9. py__cvt_adj_ens_pdaf\n    " \
                                                      "10. core 3DEnVar algorithm\n    " \
                                                      "After the iterations: \n    " \
                                                      "11. py__cvt_ens_pdaf\n    " \
                                                      "Perform ESTKF: " \
                                                      "12. py__init_dim_obs_pdaf\n    " \
                                                      "13. py__obs_op_pdaf (for ensemble mean\n    " \
                                                      "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                      "15. py__prodRinvA_pdaf\n    " \
                                                      "16. core ESTKF algorithm\n    " \
                                                      "17. py__prepoststep_state_pdaf\n    " \
                                                      "18. py__distribute_state_pdaf\n    " \
                                                      "19. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_en3dvar_lestkf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency.\n    " \
                                                       "I.e., `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                       "Using 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                                       "The background error covariance matrix is estimated by ensemble. " \
                                                       "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                       "An LESTKF is used to generate ensemble perturbations. " \
                                                       "This function should be called at each model time step. \n    \n    " \
                                                       "The function is a combination of `pyPDAF.PDAF.omi_put_state_en3dvar_lestkf_nondiagR`\n    " \
                                                       "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                       "in the following sequence: \n    " \
                                                       "1. py__collect_state_pdaf\n    " \
                                                       "2. py__prepoststep_state_pdaf\n    " \
                                                       "3. py__init_dim_obs_pdaf\n    " \
                                                       "4. py__obs_op_pdaf\n    " \
                                                       "Starting the iterative optimisation:\n    " \
                                                       "5. py__cvt_ens_pdaf\n    " \
                                                       "6. py__obs_op_lin_pdaf\n    " \
                                                       "7. py__prodRinvA_pdaf\n    " \
                                                       "8. py__obs_op_adj_pdaf\n    " \
                                                       "9. py__cvt_adj_ens_pdaf\n    " \
                                                       "10. core DA algorithm\n    " \
                                                       "After the iterations: \n    " \
                                                       "11. py__cvt_ens_pdaf\n    " \
                                                       "Perform LESTKF: \n    " \
                                                       "12. py__init_n_domains_p_pdaf\n    " \
                                                       "13. py__init_dim_obs_pdaf\n    " \
                                                       "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                       "loop over each local domain:\n    " \
                                                       "15. py__init_dim_l_pdaf\n    " \
                                                       "16. py__init_dim_obs_l_pdaf\n    " \
                                                       "17. py__g2l_state_pdaf (localise mean ensemble in observation space)\n    " \
                                                       "18. py__prodRinvA_l_pdaf\n    " \
                                                       "19. core DA algorithm\n    " \
                                                       "20. py__l2g_state_pdaf\n    " \
                                                       "21. py__prepoststep_state_pdaf\n    " \
                                                       "22. py__distribute_state_pdaf\n    " \
                                                       "23. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_enkf_nondiagR'] = "Using stochastic EnKF (ensemble " \
                                             "Kalman filter) for DA with non-diagonal observation error covariance matrix. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function is a combination of `pyPDAF.PDAF.omi_put_state_enkf_nondiagR` " \
                                             "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_dim_obs_pdaf\n    " \
                                             "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                             "5. py__add_obs_err_pdaf\n    " \
                                             "6. py__init_obscovar_pdaf\n    " \
                                             "7. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "8. core DA algorithm\n    " \
                                             "9. py__prepoststep_state_pdaf\n    " \
                                             "10. py__distribute_state_pdaf\n    " \
                                             "11. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_global_nondiagR'] = "Using global filters for DA except for 3DVars with non-diagonal observation error covariance matrix.\n    " \
                                               "The function is a combination of `pyPDAF.PDAF.omi_put_state_global_nondiagR` " \
                                               "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_dim_obs_pdaf\n    " \
                                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                               "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "6. py__prodRinvA_pdaf\n    " \
                                               "7. core DA algorithm\n    " \
                                               "8. py__prepoststep_state_pdaf\n    " \
                                               "9. py__distribute_state_pdaf\n    " \
                                               "10. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_hyb3dvar_estkf_nondiagR'] = "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                          "Here, the background error covariance is hybridised by a static background error covariance, " \
                                          "and a flow-dependent background error covariance estimated from ensemble. " \
                                          "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                          "ESTKF in this implementation. " \
                                          "This function should be called at each model time step. \n    \n    " \
                                          "The function is a combination of `pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`\n    " \
                                          "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                          "in the following sequence: \n    " \
                                          "1. py__collect_state_pdaf\n    " \
                                          "2. py__prepoststep_state_pdaf\n    " \
                                          "3. py__init_dim_obs_pdaf\n    " \
                                          "4. py__obs_op_pdaf\n    " \
                                          "Starting the iterative optimisation:\n    " \
                                          "5. py__cvt_pdaf\n    " \
                                          "6. py__cvt_ens_pdaf\n    " \
                                          "7. py__obs_op_lin_pdaf\n    " \
                                          "8. py__prodRinvA_pdaf\n    " \
                                          "9. py__obs_op_adj_pdaf\n    " \
                                          "10. py__cvt_adj_pdaf\n    " \
                                          "11. py__cvt_adj_ens_pdaf\n    " \
                                          "12. core 3DEnVar algorithm\n    " \
                                          "After the iterations: \n    " \
                                          "13. py__cvt_pdaf\n    " \
                                          "14. py__cvt_ens_pdaf\n    " \
                                          "Perform ESTKF: " \
                                          "15. py__init_dim_obs_pdaf\n    " \
                                          "16. py__obs_op_pdaf (for ensemble mean\n    " \
                                          "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                          "18. py__prodRinvA_pdaf\n    " \
                                          "19. core ESTKF algorithm\n    " \
                                          "20. py__prepoststep_state_pdaf\n    " \
                                          "21. py__distribute_state_pdaf\n    " \
                                          "22. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_hyb3dvar_lestkf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency.\n    " \
                                                        "I.e., `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                        "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                                        "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                        "and a flow-dependent background error covariance estimated from ensemble. " \
                                                        "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                        "LESTKF in this implementation. " \
                                                        "This function should be called at each model time step. \n    \n    " \
                                                        "The function is a combination of `pyPDAF.PDAF.omi_put_state_hyb3dvar_lestkf_nondiagR`\n    " \
                                                        "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                        "in the following sequence: \n    " \
                                                        "1. py__collect_state_pdaf\n    " \
                                                        "2. py__prepoststep_state_pdaf\n    " \
                                                        "3. py__init_dim_obs_pdaf\n    " \
                                                        "4. py__obs_op_pdaf\n    " \
                                                        "Starting the iterative optimisation:\n    " \
                                                        "5. py__cvt_pdaf\n    " \
                                                        "6. py__cvt_ens_pdaf\n    " \
                                                        "7. py__obs_op_lin_pdaf\n    " \
                                                        "8. py__prodRinvA_pdaf\n    " \
                                                        "9. py__obs_op_adj_pdaf\n    " \
                                                        "10. py__cvt_adj_pdaf\n    " \
                                                        "11. py__cvt_adj_ens_pdaf\n    " \
                                                        "12. core DA algorithm\n    " \
                                                        "After the iterations: \n    " \
                                                        "13. py__cvt_pdaf\n    " \
                                                        "14. py__cvt_ens_pdaf\n    " \
                                                        "Perform LESTKF: \n    " \
                                                        "15. py__init_n_domains_p_pdaf\n    " \
                                                        "16. py__init_dim_obs_pdaf\n    " \
                                                        "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                                        "loop over each local domain:\n    " \
                                                        "18. py__init_dim_l_pdaf\n    " \
                                                        "19. py__init_dim_obs_l_pdaf\n    " \
                                                        "20. py__g2l_state_pdaf\n    " \
                                                        "21. py__init_obs_l_pdaf\n    "\
                                                        "22. py__prodRinvA_l_pdaf\n    " \
                                                        "23. core DA algorithm\n    " \
                                                        "24. py__l2g_state_pdaf\n    " \
                                                        "25. py__prepoststep_state_pdaf\n    " \
                                                        "26. py__distribute_state_pdaf\n    " \
                                                        "27. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_lenkf_nondiagR'] = "Using stochastic EnKF (ensemble " \
                                              "Kalman filter) with covariance localisation " \
                                              "for DA with non-diagonal observation error covariance matrix. " \
                                              "This is the only scheme for covariance localisation in PDAF. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_lenkf_nondiagR` " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "5. py__localize_pdaf\n    " \
                                              "6. py__add_obs_err_pdaf\n    " \
                                              "7. py__init_obscovar_pdaf\n    " \
                                              "8. py__obs_op_pdaf (repeated to reduce storage\n    " \
                                              "9. core DA algorith\n    " \
                                              "10. py__prepoststep_state_pdaf\n    " \
                                              "11. py__distribute_state_pdaf\n    " \
                                              "12. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_lknetf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                               "I.e., `pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`. \n    " \
                                               "Using LKNETF for DA with non-diagonal observation error covariance matrix. " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "The function is a combination of `pyPDAF.PDAF.omi_put_state_lknetf_nondiagR` " \
                                               "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_n_domains_p_pdaf\n    " \
                                               "4. py__init_dim_obs_pdaf\n    " \
                                               "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "loop over each local domain:\n    " \
                                               "6. py__init_dim_l_pdaf\n    " \
                                               "7. py__init_dim_obs_l_pdaf\n    " \
                                               "8. py__g2l_state_pdaf\n    " \
                                               "9. py__prodRinvA_pdaf\n    " \
                                               "10. py__likelihood_l_pdaf\n    " \
                                               "11. core DA algorithm\n    " \
                                               "12. py__l2g_state_pdaf\n    " \
                                               "13. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                               "14. py__likelihood_hyb_l_pda\n    " \
                                               "15. py__prodRinvA_hyb_l_pdaf\n    " \
                                               "16. py__prepoststep_state_pdaf\n    " \
                                               "17. py__distribute_state_pdaf\n    " \
                                               "18. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_lnetf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                              "I.e., `pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`. \n    " \
                                              "Using LNETF for DA with non-diagonal observation error covariance matrix. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_lnetf_nondiagR` " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_n_domains_p_pdaf\n    " \
                                              "4. py__init_dim_obs_pdaf\n    " \
                                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "loop over each local domain:\n    " \
                                              "6. py__init_dim_l_pdaf\n    " \
                                              "7. py__init_dim_obs_l_pdaf\n    " \
                                              "8. py__g2l_state_pdaf\n    " \
                                              "9. py__likelihood_l_pdaf\n    " \
                                              "10. core DA algorithm\n    " \
                                              "11. py__l2g_state_pdaf\n    " \
                                              "12. py__prepoststep_state_pdaf\n    " \
                                              "13. py__distribute_state_pdaf\n    " \
                                              "14. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_local_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                              "I.e., `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                              "Using domain localised filters for DA with non-diagonal observation error covariance matrix. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function is a combination of `pyPDAF.PDAF.omi_put_state_local_nondiagR` " \
                                              "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_n_domains_p_pdaf\n    " \
                                              "4. py__init_dim_obs_pdaf\n    " \
                                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "loop over each local domain:\n    " \
                                              "6. py__init_dim_l_pdaf\n    " \
                                              "7. py__init_dim_obs_l_pdaf\n    " \
                                              "8. py__g2l_state_pdaf\n    " \
                                              "9. py__init_obs_l_pdaf\n    "\
                                              "10. py__prodRinvA_l_pdaf\n    " \
                                              "11. core DA algorithm\n    " \
                                              "12. py__l2g_state_pdaf\n    " \
                                              "13. py__prepoststep_state_pdaf\n    " \
                                              "14. py__distribute_state_pdaf\n    " \
                                              "15. py__next_observation_pdaf\n    "
docstrings['omi_assimilate_nonlin_nondiagR'] = "Using nonlinear filters (particle filter, NETF) for DA except for 3DVars with non-diagonal observation error covariance matrix.\n    " \
                                               "The function is a combination of `pyPDAF.PDAF.put_state_nonlin_nondiagR` " \
                                               "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_dim_obs_pdaf\n    " \
                                               "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                               "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "6. py__likelihood_pdaf\n    " \
                                               "7. core DA algorithm\n    " \
                                               "8. py__prepoststep_state_pdaf\n    " \
                                               "9. py__distribute_state_pdaf\n    " \
                                               "10. py__next_observation_pdaf\n    "

docstrings['omi_put_state_3dvar_nondiagR'] = "Using 3DVar for DA with non-diagonal observation error covariance matrix\n    " \
                                             "without post-processing and analysis distribution to forecsat without OMI. " \
                                             "This is a deterministic filtering scheme. " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_dim_obs_pdaf\n    " \
                                             "4. py__obs_op_pdaf\n    " \
                                             "Starting the iterative optimisation:\n    " \
                                             "5. py__cvt_pdaf\n    " \
                                             "6. py__obs_op_lin_pdaf\n    " \
                                             "7. py__prodRinvA_pdaf\n    " \
                                             "8. py__obs_op_adj_pdaf\n    " \
                                             "9. py__cvt_adj_pdaf\n    " \
                                             "10. core DA algorithm\n    " \
                                             "After the iterations: \n    " \
                                             "11. py__cvt_pdaf\n    "
docstrings['omi_put_state_en3dvar_estkf_nondiagR'] = "Using 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                     "without post-processing and analysis distribution to forecsat without OMI. " \
                                                     "The background error covariance matrix is estimated by ensemble. " \
                                                     "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                     "An ESTKF is used to generate ensemble perturbations. " \
                                                     "This function is usually used in 'flexible' parallelisation. " \
                                                     "i.e., the ensemble size is larger than the available number of processes. " \
                                                     "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                     "to the model after this function. " \
                                                     "This function should be called at each model time step. \n    \n    " \
                                                     "The function executes the user-supplied function " \
                                                     "in the following sequence: \n    " \
                                                     "1. py__collect_state_pdaf\n    " \
                                                     "2. py__prepoststep_state_pdaf\n    " \
                                                     "3. py__init_dim_obs_pdaf\n    " \
                                                     "4. py__obs_op_pdaf\n    " \
                                                     "Starting the iterative optimisation:\n    " \
                                                     "5. py__cvt_ens_pdaf\n    " \
                                                     "6. py__obs_op_lin_pdaf\n    " \
                                                     "7. py__prodRinvA_pdaf\n    " \
                                                     "8. py__obs_op_adj_pdaf\n    " \
                                                     "9. py__cvt_adj_ens_pdaf\n    " \
                                                     "10. core 3DEnVar algorithm\n    " \
                                                     "After the iterations: \n    " \
                                                     "11. py__cvt_ens_pdaf\n    " \
                                                     "Perform ESTKF: " \
                                                     "12. py__init_dim_obs_pdaf\n    " \
                                                     "13. py__obs_op_pdaf (for ensemble mean\n    " \
                                                     "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                     "15. py__prodRinvA_pdaf\n    " \
                                                     "16. core ESTKF algorithm\n    "
docstrings['omi_put_state_en3dvar_lestkf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency.\n    " \
                                                      "I.e., `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                      "Using 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                      "without post-processing and analysis distribution to forecsat without OMI. " \
                                                      "The background error covariance matrix is estimated by ensemble. " \
                                                      "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                      "An LESTKF is used to generate ensemble perturbations. " \
                                                      "This function is usually used in 'flexible' parallelisation. " \
                                                      "i.e., the ensemble size is larger than the available number of processes. " \
                                                      "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                      "to the model after this function. " \
                                                      "This function should be called at each model time step. \n    \n    " \
                                                      "The function executes the user-supplied function " \
                                                      "in the following sequence: \n    " \
                                                      "1. py__collect_state_pdaf\n    " \
                                                      "2. py__prepoststep_state_pdaf\n    " \
                                                      "3. py__init_dim_obs_pdaf\n    " \
                                                      "4. py__obs_op_pdaf\n    " \
                                                      "Starting the iterative optimisation:\n    " \
                                                      "5. py__cvt_ens_pdaf\n    " \
                                                      "6. py__obs_op_lin_pdaf\n    " \
                                                      "7. py__prodRinvA_pdaf\n    " \
                                                      "8. py__obs_op_adj_pdaf\n    " \
                                                      "9. py__cvt_adj_ens_pdaf\n    " \
                                                      "10. core DA algorithm\n    " \
                                                      "After the iterations: \n    " \
                                                      "11. py__cvt_ens_pdaf\n    " \
                                                      "Perform LESTKF: \n    " \
                                                      "12. py__init_n_domains_p_pdaf\n    " \
                                                      "13. py__init_dim_obs_pdaf\n    " \
                                                      "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                      "loop over each local domain:\n    " \
                                                      "15. py__init_dim_l_pdaf\n    " \
                                                      "16. py__init_dim_obs_l_pdaf\n    " \
                                                      "17. py__g2l_state_pdaf (localise mean ensemble in observation space)\n    " \
                                                      "18. py__prodRinvA_l_pdaf\n    " \
                                                      "19. core DA algorithm\n    " \
                                                      "20. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_enkf_nondiagR'] =  "Using stochastic EnKF (ensemble " \
                                             "Kalman filter) for DA with non-diagonal observation error covariance matrix. " \
                                             "without post-processing and analysis distribution to forecsat without OMI. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_dim_obs_pdaf\n    " \
                                             "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                             "5. py__add_obs_err_pdaf\n    " \
                                             "6. py__init_obscovar_pdaf\n    " \
                                             "7. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "8. core DA algorithm\n    "
docstrings['omi_put_state_global_nondiagR'] = "Using global filters for DA except for 3DVars with non-diagonal observation error covariance matrix\n    " \
                                              "without post-processing and analysis distribution to forecsat without OMI. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "This function is usually used in 'flexible' parallelisation. " \
                                              "i.e., the ensemble size is larger than the available number of processes. " \
                                              "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                              "to the model after this function. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "6. py__prodRinvA_pdaf\n    " \
                                              "7. core DA algorithm\n    "
docstrings['omi_put_state_hyb3dvar_estkf_nondiagR'] = "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                      "without post-processing and analysis distribution to forecsat without OMI. " \
                                                      "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                      "and a flow-dependent background error covariance estimated from ensemble. " \
                                                      "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                      "ESTKF in this implementation. " \
                                                      "This function should be called at each model time step. \n    \n    " \
                                                      "This function is usually used in 'flexible' parallelisation. " \
                                                      "i.e., the ensemble size is larger than the available number of processes. " \
                                                      "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                      "to the model after this function. " \
                                                      "This function should be called at each model time step. \n    \n    " \
                                                      "The function executes the user-supplied function " \
                                                      "in the following sequence: \n    " \
                                                      "1. py__collect_state_pdaf\n    " \
                                                      "2. py__prepoststep_state_pdaf\n    " \
                                                      "3. py__init_dim_obs_pdaf\n    " \
                                                      "4. py__obs_op_pdaf\n    " \
                                                      "Starting the iterative optimisation:\n    " \
                                                      "5. py__cvt_pdaf\n    " \
                                                      "6. py__cvt_ens_pdaf\n    " \
                                                      "7. py__obs_op_lin_pdaf\n    " \
                                                      "8. py__prodRinvA_pdaf\n    " \
                                                      "9. py__obs_op_adj_pdaf\n    " \
                                                      "10. py__cvt_adj_pdaf\n    " \
                                                      "11. py__cvt_adj_ens_pdaf\n    " \
                                                      "12. core 3DEnVar algorithm\n    " \
                                                      "After the iterations: \n    " \
                                                      "13. py__cvt_pdaf\n    " \
                                                      "14. py__cvt_ens_pdaf\n    " \
                                                      "Perform ESTKF: " \
                                                      "15. py__init_dim_obs_pdaf\n    " \
                                                      "16. py__obs_op_pdaf (for ensemble mean\n    " \
                                                      "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                                      "18. py__prodRinvA_pdaf\n    " \
                                                      "19. core ESTKF algorithm\n    "
docstrings['omi_put_state_hyb3dvar_lestkf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency.\n    " \
                                                       "I.e., `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                       "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                       "without post-processing and analysis distribution to forecsat without OMI. " \
                                                       "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                       "and a flow-dependent background error covariance estimated from ensemble. " \
                                                       "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                       "LESTKF in this implementation. " \
                                                       "This function should be called at each model time step. \n    \n    " \
                                                       "This function is usually used in 'flexible' parallelisation. " \
                                                       "i.e., the ensemble size is larger than the available number of processes. " \
                                                       "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                       "to the model after this function. " \
                                                       "This function should be called at each model time step. \n    \n    " \
                                                       "The function executes the user-supplied function " \
                                                       "in the following sequence: \n    " \
                                                       "1. py__collect_state_pdaf\n    " \
                                                       "2. py__prepoststep_state_pdaf\n    " \
                                                       "3. py__init_dim_obs_pdaf\n    " \
                                                       "4. py__obs_op_pdaf\n    " \
                                                       "Starting the iterative optimisation:\n    " \
                                                       "5. py__cvt_pdaf\n    " \
                                                       "6. py__cvt_ens_pdaf\n    " \
                                                       "7. py__obs_op_lin_pdaf\n    " \
                                                       "8. py__prodRinvA_pdaf\n    " \
                                                       "9. py__obs_op_adj_pdaf\n    " \
                                                       "10. py__cvt_adj_pdaf\n    " \
                                                       "11. py__cvt_adj_ens_pdaf\n    " \
                                                       "12. core DA algorithm\n    " \
                                                       "After the iterations: \n    " \
                                                       "13. py__cvt_pdaf\n    " \
                                                       "14. py__cvt_ens_pdaf\n    " \
                                                       "Perform LESTKF: \n    " \
                                                       "15. py__init_n_domains_p_pdaf\n    " \
                                                       "16. py__init_dim_obs_pdaf\n    " \
                                                       "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                                       "loop over each local domain:\n    " \
                                                       "18. py__init_dim_l_pdaf\n    " \
                                                       "19. py__init_dim_obs_l_pdaf\n    " \
                                                       "20. py__g2l_state_pdaf\n    " \
                                                       "21. py__init_obs_l_pdaf\n    "\
                                                       "22. py__prodRinvA_l_pdaf\n    " \
                                                       "23. core DA algorithm\n    " \
                                                       "24. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_lenkf_nondiagR'] = "Using stochastic EnKF (ensemble " \
                                             "Kalman filter) with covariance localisation " \
                                             "for DA with non-diagonal observation error covariance matrix " \
                                             "without post-processing and analysis distribution to forecsat without OMI. " \
                                             "This is the only scheme for covariance localisation in PDAF. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_dim_obs_pdaf\n    " \
                                             "4. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "5. py__localize_pdaf\n    " \
                                             "6. py__add_obs_err_pdaf\n    " \
                                             "7. py__init_obscovar_pdaf\n    " \
                                             "8. py__obs_op_pdaf (repeated to reduce storage\n    " \
                                             "9. core DA algorith\n    "
docstrings['omi_put_state_lknetf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                              "I.e., `pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`. \n    " \
                                              "Using LKNETF for DA with non-diagonal observation error covariance matrix " \
                                              "without post-processing and analysis distribution to forecsat without OMI. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "This function is usually used in 'flexible' parallelisation. " \
                                              "i.e., the ensemble size is larger than the available number of processes. " \
                                              "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                              "to the model after this function. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_n_domains_p_pdaf\n    " \
                                              "4. py__init_dim_obs_pdaf\n    " \
                                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "loop over each local domain:\n    " \
                                              "6. py__init_dim_l_pdaf\n    " \
                                              "7. py__init_dim_obs_l_pdaf\n    " \
                                              "8. py__g2l_state_pdaf\n    " \
                                              "9. py__prodRinvA_pdaf\n    " \
                                              "10. py__likelihood_l_pdaf\n    " \
                                              "11. core DA algorithm\n    " \
                                              "12. py__l2g_state_pdaf\n    " \
                                              "13. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                              "14. py__likelihood_hyb_l_pda\n    " \
                                              "15. py__prodRinvA_hyb_l_pdaf\n    "
docstrings['omi_put_state_lnetf_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                             "I.e., `pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`. \n    " \
                                             "Using LNETF for DA with non-diagonal observation error covariance matrix " \
                                             "without post-processing and analysis distribution to forecsat without OMI. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function is a combination of `pyPDAF.PDAF.omi_put_state_lnetf_nondiagR` " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_n_domains_p_pdaf\n    " \
                                             "4. py__init_dim_obs_pdaf\n    " \
                                             "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "loop over each local domain:\n    " \
                                             "6. py__init_dim_l_pdaf\n    " \
                                             "7. py__init_dim_obs_l_pdaf\n    " \
                                             "8. py__g2l_state_pdaf\n    " \
                                             "9. py__likelihood_l_pdaf\n    " \
                                             "10. core DA algorithm\n    " \
                                             "11. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_local_nondiagR'] = "It is recommended to use local module for fewer user-supplied functions and improved efficiency. " \
                                             "I.e., `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                             "Using domain localised filters for DA with non-diagonal observation error covariance matrix " \
                                             "without post-processing and analysis distribution to forecsat without OMI. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "This function is usually used in 'flexible' parallelisation. " \
                                             "i.e., the ensemble size is larger than the available number of processes. " \
                                             "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                             "to the model after this function. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_n_domains_p_pdaf\n    " \
                                             "4. py__init_dim_obs_pdaf\n    " \
                                             "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "loop over each local domain:\n    " \
                                             "6. py__init_dim_l_pdaf\n    " \
                                             "7. py__init_dim_obs_l_pdaf\n    " \
                                             "8. py__g2l_state_pdaf\n    " \
                                             "9. py__init_obs_l_pdaf\n    "\
                                             "10. py__prodRinvA_l_pdaf\n    " \
                                             "11. core DA algorithm\n    " \
                                             "12. py__l2g_state_pdaf\n    "
docstrings['omi_put_state_nonlin_nondiagR'] = "Using nonlinear filters (particle filter, NETF) for DA except for 3DVars with non-diagonal observation error covariance matrix\n    " \
                                              "without post-processing and analysis distribution to forecsat without OMI. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "This function is usually used in 'flexible' parallelisation. " \
                                              "i.e., the ensemble size is larger than the available number of processes. " \
                                              "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                              "to the model after this function. " \
                                              "This function should be called at each model time step. \n    \n    " \
                                              "The function executes the user-supplied function " \
                                              "in the following sequence: \n    " \
                                              "1. py__collect_state_pdaf\n    " \
                                              "2. py__prepoststep_state_pdaf\n    " \
                                              "3. py__init_dim_obs_pdaf\n    " \
                                              "4. py__obs_op_pdaf (for ensemble mean\n    " \
                                              "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                              "6. py__likelihood_pdaf\n    " \
                                              "7. core DA algorithm\n    "

docstrings['local_set_indices'] = "Set index vector to map local state vector to global state vectors. " \
                                  "This is called in the user-supplied function `py__init_dim_l_pdaf`."
docstrings['local_set_increment_weights'] = "This function initialises a PDAF_internal local array " \
                                            "of increment weights. The weights are applied in " \
                                            "in PDAF_local_l2g_cb where the local state vector\n    " \
                                            "is weighted by given weights. " \
                                            "These can e.g. be used to apply a vertical localisation."
docstrings['local_clear_increment_weights'] = "This function deallocates the local increment weight vector " \
                                              "in `pyPDAF.PDAF.local_set_increment_weights` if it is allocated"
docstrings['local_g2l_cb'] = "Project a global to a local state vector for the localized filters.\n    " \
                             "This is the full callback function to be used internally. The mapping " \
                             "is done using the index vector id_lstate_in_pstate that is initialised " \
                             "in `pyPDAF.PDAF.local_set_indices`."
docstrings['local_l2g_cb'] = "Initialise elements of a global state vector from a local state vector.\n    " \
                             "This is the full callback function to be used internally. The mapping " \
                             "is done using the index vector `id_lstate_in_pstate` that is initialised " \
                             "in `pyPDAF.PDAF.local_set_indices`. \n    \n    " \
                             "To exclude any element of the local state vector from the initialisation" \
                             "one can set the corresponding index value to 0."

docstrings['local_assimilate_en3dvar_lestkf'] =  "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                                 "I.e., `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf` or `pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                 "Using 3DEnVar for DA without OMI.\n    " \
                                                 "The background error covariance matrix is estimated by ensemble. " \
                                                 "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                 "An LESTKF is used to generate ensemble perturbations. " \
                                                 "This function should be called at each model time step. \n    \n    " \
                                                 "The function is a combination of `pyPDAF.PDAF.local_put_state_en3dvar_lestkf`\n    " \
                                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                 "in the following sequence: \n    " \
                                                 "1. py__collect_state_pdaf\n    " \
                                                 "2. py__prepoststep_state_pdaf\n    " \
                                                 "3. py__init_dim_obs_pdaf\n    " \
                                                 "4. py__obs_op_pdaf\n    " \
                                                 "5. py__init_obs_pdaf\n    " \
                                                 "Starting the iterative optimisation:\n    " \
                                                 "6. py__cvt_ens_pdaf\n    " \
                                                 "7. py__obs_op_lin_pdaf\n    " \
                                                 "8. py__prodRinvA_pdaf\n    " \
                                                 "9. py__obs_op_adj_pdaf\n    " \
                                                 "10. py__cvt_adj_ens_pdaf\n    " \
                                                 "11. core DA algorithm\n    " \
                                                 " After the iterations: n    " \
                                                 "12. py__cvt_ens_pdaf\n    " \
                                                 "Perform LESTKF: \n    " \
                                                 "13. py__init_n_domains_p_pdaf\n    " \
                                                 "14. py__init_dim_obs_pdaf\n    " \
                                                 "15. py__obs_op_pdaf (for each ensemble member\n    " \
                                                 "16. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                                 "17. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                                 "loop over each local domain:\n    " \
                                                 "18. py__init_dim_l_pdaf\n    " \
                                                 "19. py__init_dim_obs_l_pdaf\n    " \
                                                 "20. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                                 "21. py__init_obs_l_pdaf\n    "\
                                                 "22. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                                 "23. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                                 "24. py__prodRinvA_l_pdaf\n    " \
                                                 "25. core DA algorithm\n    " \
                                                 "26. py__prepoststep_state_pdaf\n    " \
                                                 "27. py__distribute_state_pdaf\n    " \
                                                 "28. py__next_observation_pdaf\n    "
docstrings['local_assimilate_hyb3dvar_lestkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                                 "I.e., `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf` or `pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                 "Using Hybrid 3DEnVar for DA without OMI.\n    " \
                                                 "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                 "and a flow-dependent background error covariance estimated from ensemble. " \
                                                 "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                 "LESTKF in this implementation. " \
                                                 "This function should be called at each model time step. \n    \n    " \
                                                 "The function is a combination of `pyPDAF.PDAF.local_put_state_hyb3dvar_lestkf`\n    " \
                                                 "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                 "in the following sequence: \n    " \
                                                 "1. py__collect_state_pdaf\n    " \
                                                 "2. py__prepoststep_state_pdaf\n    " \
                                                 "3. py__init_dim_obs_pdaf\n    " \
                                                 "4. py__obs_op_pdaf\n    " \
                                                 "5. py__init_obs_pdaf\n    " \
                                                 "Starting the iterative optimisation:\n    " \
                                                 "6. py__cvt_pdaf\n    " \
                                                 "7. py__cvt_ens_pdaf\n    " \
                                                 "8. py__obs_op_lin_pdaf\n    " \
                                                 "9. py__prodRinvA_pdaf\n    " \
                                                 "10. py__obs_op_adj_pdaf\n    " \
                                                 "11. py__cvt_adj_pdaf\n    " \
                                                 "12. py__cvt_adj_ens_pdaf\n    " \
                                                 "13. core DA algorithm\n    " \
                                                 " After the iterations: n    " \
                                                 "14. py__cvt_pdaf\n    " \
                                                 "15. py__cvt_ens_pdaf\n    " \
                                                 "Perform LESTKF: \n    " \
                                                 "16. py__init_n_domains_p_pdaf\n    " \
                                                 "17. py__init_dim_obs_pdaf\n    " \
                                                 "18. py__obs_op_pdaf (for each ensemble member\n    " \
                                                 "19. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in `pyPDAF.PDAF.init`))\n    " \
                                                 "20. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                                 "loop over each local domain:\n    " \
                                                 "21. py__init_dim_l_pdaf\n    " \
                                                 "22. py__init_dim_obs_l_pdaf\n    " \
                                                 "23. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                                 "24. py__init_obs_l_pdaf\n    "\
                                                 "25. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                                 "26. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)\n    "\
                                                 "27. py__prodRinvA_l_pdaf\n    " \
                                                 "28. core DA algorithm\n    " \
                                                 "29. py__prepoststep_state_pdaf\n    " \
                                                 "30. py__distribute_state_pdaf\n    " \
                                                 "31. py__next_observation_pdaf\n    "
docstrings['local_assimilate_lestkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                        "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                        "Using Local ESTKF (error space transform " \
                                        "Kalman filter) for DA without OMI. " \
                                        "This is a domain localisation method. " \
                                        "This function should be called at each model time step. " \
                                        "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                        "The function is a combination of `pyPDAF.PDAF.local_put_state_lestkf` " \
                                        "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                        "in the following sequence: \n    " \
                                        "1. py__collect_state_pdaf\n    " \
                                        "2. py__prepoststep_state_pdaf\n    " \
                                        "3. py__init_n_domains_p_pdaf\n    " \
                                        "4. py__init_dim_obs_pdaf\n    " \
                                        "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                        "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                        "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                        "loop over each local domain:\n    " \
                                        "8. py__init_dim_l_pdaf\n    " \
                                        "9. py__init_dim_obs_l_pdaf\n    " \
                                        "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                        "11. py__init_obs_l_pdaf\n    "\
                                        "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                        "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                        "14. py__prodRinvA_l_pdaf\n    " \
                                        "15. core DA algorithm\n    " \
                                        "16. py__prepoststep_state_pdaf\n    " \
                                        "17. py__distribute_state_pdaf (18. py__next_observation_pdaf \n    "
docstrings['local_assimilate_letkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                       "Using local ensemble transform Kalman filter for DA without OMI. " \
                                       "This is a domain localisation method. " \
                                       "This function should be called at each model time step. \n    \n    " \
                                       "The function is a combination of `pyPDAF.PDAF.local_put_state_letkf` " \
                                       "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_n_domains_p_pdaf\n    " \
                                       "4. py__init_dim_obs_pdaf\n    " \
                                       "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                       "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                       "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                       "loop over each local domain:\n    " \
                                       "8. py__init_dim_l_pdaf\n    " \
                                       "9. py__init_dim_obs_l_pdaf\n    " \
                                       "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                       "11. py__init_obs_l_pdaf\n    "\
                                       "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                       "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                       "14. py__prodRinvA_l_pdaf\n    " \
                                       "15. core DA algorithm\n    " \
                                       "16. py__prepoststep_state_pdaf\n    " \
                                       "17. py__distribute_state_pdaf\n    " \
                                       "18. py__next_observation_pdaf \\n    "
docstrings['local_assimilate_lknetf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                        "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`. \n    " \
                                        "This function will is a hybridised LETKF and LNETF " \
                                        "for DA without OMI. The LNETF computes the distribution up to " \
                                        "the second moment similar to KF but using a nonlinear weighting similar to " \
                                        "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                        "The hybridisation with LETKF is expected to lead to improved performance for " \
                                        "quasi-Gaussian problems. " \
                                        "The function should be called at each model step. \n    \n    " \
                                        "The function is a combination of `pyPDAF.PDAF.local_put_state_lknetf` " \
                                        "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                        "in the following sequence: \n    " \
                                        "1. py__collect_state_pdaf \\n    " \
                                        "2. py__prepoststep_state_pdaf \\n    " \
                                        "3. py__init_n_domains_p_pdaf \\n    " \
                                        "4. py__init_dim_obs_pdaf \\n    " \
                                        "5. py__obs_op_pdaf (for each ensemble member)\\n    " \
                                        "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in `pyPDAF.PDAF.init`\n    " \
                                        "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                        "loop over each local domain:\n    " \
                                        "8. py__init_dim_l_pdaf \\n    " \
                                        "9. py__init_dim_obs_l_pdaf \\n    " \
                                        "10. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                        "11. py__init_obs_l_pdaf \\n    "\
                                        "12. py__init_obsvar_l_pdaf \n(only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                        "13. py__prodRinvA_pdaf \\n    " \
                                        "14. py__likelihood_l_pdaf \\n    " \
                                        "15. core DA algorithm \\n    " \
                                        "16. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member)\\n    " \
                                        "17. py__likelihood_hyb_l_pdaf\\n    " \
                                        "18. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                        "19. py__prodRinvA_hyb_l_pdaf \\n    " \
                                        "20. py__prepoststep_state_pdaf \\n    " \
                                        "21. py__distribute_state_pdaf \\n    " \
                                        "22. py__next_observation_pdaf\n    "
docstrings['local_assimilate_lnetf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`. \n    " \
                                       "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                       "for DA without OMI. The nonlinear filter computes the distribution up to " \
                                       "the second moment similar to KF but using a nonlinear weighting similar to " \
                                       "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                       "This function should be called at each model time step. \n    \n    " \
                                       "The function is a combination of `pyPDAF.PDAF.local_put_state_lnetf` " \
                                       "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_n_domains_p_pdaf\n    " \
                                       "4. py__init_dim_obs_pdaf\n    " \
                                       "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                       "loop over each local domain:\n    " \
                                       "6. py__init_dim_l_pdaf\n    " \
                                       "7. py__init_dim_obs_l_pdaf\n    " \
                                       "8. py__init_obs_l_pdaf\n    "\
                                       "9. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                       "10. py__likelihood_l_pdaf\n    " \
                                       "11. core DA algorithm\n    " \
                                       "12. py__prepoststep_state_pdaf\n    " \
                                       "13. py__distribute_state_pdaf\n    " \
                                       "14. py__next_observation_pdaf\n    "
docstrings['local_assimilate_lseik'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.localomi_assimilate` or `pyPDAF.PDAF.localomi_assimilate_nondiagR`. \n    " \
                                       "Using local singular evolutive interpolated Kalman filter for DA without OMI. " \
                                       "This is a domain localisation method. " \
                                       "This function should be called at each model time step.\n    \n    " \
                                       "The function is a combination of `pyPDAF.PDAF.local_put_state_lseik` " \
                                       "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_n_domains_p_pdaf\n    " \
                                       "4. py__init_dim_obs_pdaf\n    " \
                                       "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                       "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                       "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                       "loop over each local domain:\n    " \
                                       "8. py__init_dim_l_pdaf\n    " \
                                       "9. py__init_dim_obs_l_pdaf\n    " \
                                       "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                       "11. py__init_obs_l_pdaf\n    "\
                                       "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                       "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                       "14. py__prodRinvA_l_pdaf\n    " \
                                       "15. core DA algorithm\n    " \
                                       "16. py__prepoststep_state_pdaf\n    " \
                                       "17. py__distribute_state_pdaf\n    " \
                                       "18. py__next_observation_pdaf\n    "

docstrings['local_put_state_en3dvar_lestkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                               "I.e., `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf` or `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                               "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                               "The background error covariance matrix is estimated by ensemble. " \
                                               "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                               "An LESTKF is used to generate ensemble perturbations. " \
                                               "This function is usually used in 'flexible' parallelisation. " \
                                               "i.e., the ensemble size is larger than the available number of processes. " \
                                               "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                               "to the model after this function. " \
                                               "This function should be called at each model time step. \n    \n    " \
                                               "The function executes the user-supplied function " \
                                               "in the following sequence: \n    " \
                                               "1. py__collect_state_pdaf\n    " \
                                               "2. py__prepoststep_state_pdaf\n    " \
                                               "3. py__init_dim_obs_pdaf\n    " \
                                               "4. py__obs_op_pdaf\n    " \
                                               "5. py__init_obs_pdaf\n    " \
                                               "Starting the iterative optimisation:\n    " \
                                               "6. py__cvt_ens_pdaf\n    " \
                                               "7. py__obs_op_lin_pdaf\n    " \
                                               "8. py__prodRinvA_pdaf\n    " \
                                               "9. py__obs_op_adj_pdaf\n    " \
                                               "10. py__cvt_adj_ens_pdaf\n    " \
                                               "11. core DA algorithm\n    " \
                                               " After the iterations: n    " \
                                               "12. py__cvt_ens_pdaf\n    " \
                                               "Perform LESTKF: \n    " \
                                               "13. py__init_n_domains_p_pdaf\n    " \
                                               "14. py__init_dim_obs_pdaf\n    " \
                                               "15. py__obs_op_pdaf (for each ensemble member\n    " \
                                               "16. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                               "17. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                               "loop over each local domain:\n    " \
                                               "18. py__init_dim_l_pdaf\n    " \
                                               "19. py__init_dim_obs_l_pdaf\n    " \
                                               "20. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                               "21. py__init_obs_l_pdaf\n    "\
                                               "22. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                               "23. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                               "24. py__prodRinvA_pdaf\n    " \
                                               "25. core DA algorithm\n    "
docstrings['local_put_state_hyb3dvar_lestkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency.\n    " \
                                                "I.e., `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf` or `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`. \n    \n    \n    " \
                                                "Using Hybrid 3DEnVar for DA without post-processing and analysis distribution to forecsat without OMI.\n    " \
                                                "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                "and a flow-dependent background error covariance estimated from ensemble. " \
                                                "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                "LESTKF in this implementation. \n    \n    " \
                                                "This function is usually used in 'flexible' parallelisation. " \
                                                "i.e., the ensemble size is larger than the available number of processes. " \
                                                "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                "to the model after this function. " \
                                                "This function should be called at each model time step. \n    \n    " \
                                                "The function executes the user-supplied function " \
                                                "in the following sequence: \n    " \
                                                "1. py__collect_state_pdaf\n    " \
                                                "2. py__prepoststep_state_pdaf\n    " \
                                                "3. py__init_dim_obs_pdaf\n    " \
                                                "4. py__obs_op_pdaf\n    " \
                                                "5. py__init_obs_pdaf\n    " \
                                                "Starting the iterative optimisation:\n    " \
                                                "6. py__cvt_pdaf\n    " \
                                                "7. py__cvt_ens_pdaf\n    " \
                                                "8. py__obs_op_lin_pdaf\n    " \
                                                "9. py__prodRinvA_pdaf\n    " \
                                                "10. py__obs_op_adj_pdaf\n    " \
                                                "11. py__cvt_adj_pdaf\n    " \
                                                "12. py__cvt_adj_ens_pdaf\n    " \
                                                "13. core DA algorithm\n    " \
                                                "After the iterations: n    " \
                                                "14. py__cvt_pdaf\n    " \
                                                "15. py__cvt_ens_pdaf\n    " \
                                                "Perform LESTKF: \n    " \
                                                "16. py__init_n_domains_p_pdaf\n    " \
                                                "17. py__init_dim_obs_pdaf\n    " \
                                                "18. py__obs_op_pdaf (for each ensemble member\n    " \
                                                "19. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in `pyPDAF.PDAF.init`))\n    " \
                                                "20. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                                "loop over each local domain:\n    " \
                                                "21. py__init_dim_l_pdaf\n    " \
                                                "22. py__init_dim_obs_l_pdaf\n    " \
                                                "23. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                                "24. py__init_obs_l_pdaf\n    "\
                                                "25. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                                "26. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used)\n    "\
                                                "27. py__prodRinvA_pdaf\n    " \
                                                "28. core DA algorithm\n    "
docstrings['local_put_state_lestkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                       "Using Local ESTKF (error space transform " \
                                       "Kalman filter) for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                       "This is a domain localisation method. " \
                                       "This function is usually used in 'flexible' parallelisation. " \
                                       "i.e., the ensemble size is larger than the available number of processes. " \
                                       "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                       "to the model after this function. " \
                                       "This function should be called at each model time step. " \
                                       "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                       "The function executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_n_domains_p_pdaf\n    " \
                                       "4. py__init_dim_obs_pdaf\n    " \
                                       "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                       "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                       "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                       "loop over each local domain:\n    " \
                                       "8. py__init_dim_l_pdaf\n    " \
                                       "9. py__init_dim_obs_l_pdaf\n    " \
                                       "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                       "11. py__init_obs_l_pdaf\n    "\
                                       "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                       "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                       "14. py__prodRinvA_l_pdaf\n    " \
                                       "15. core DA algorithm\n    "
docstrings['local_put_state_letkf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                      "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                      "Using local ensemble transform Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                      "This is a domain localisation method. " \
                                      "This function is usually used in 'flexible' parallelisation. " \
                                      "i.e., the ensemble size is larger than the available number of processes. " \
                                      "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                      "to the model after this function. " \
                                      "This function should be called at each model time step. \n    \n    " \
                                      "The function executes the user-supplied function " \
                                      "in the following sequence: \n    " \
                                      "1. py__collect_state_pdaf\n    " \
                                      "2. py__prepoststep_state_pdaf\n    " \
                                      "3. py__init_n_domains_p_pdaf\n    " \
                                      "4. py__init_dim_obs_pdaf\n    " \
                                      "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                      "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                      "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                      "loop over each local domain:\n    " \
                                      "8. py__init_dim_l_pdaf\n    " \
                                      "9. py__init_dim_obs_l_pdaf\n    " \
                                      "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                      "11. py__init_obs_l_pdaf\n    "\
                                      "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                      "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                      "14. py__prodRinvA_l_pdaf\n    " \
                                      "15. core DA algorithm\n    "
docstrings['local_put_state_lknetf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                       "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`. \n    " \
                                       "This function will is a hybridised LETKF and LNETF " \
                                       "for DA without post-processing and analysis distribution to forecsat without OMI. The LNETF computes the distribution up to " \
                                       "the second moment similar to KF but using a nonlinear weighting similar to " \
                                       "particle filter. This leads to an equal weights assumption for prior ensemble. " \
                                       "The hybridisation with LETKF is expected to lead to improved performance for " \
                                       "quasi-Gaussian problems.  \n    \n    " \
                                       "This function is usually used in 'flexible' parallelisation. " \
                                       "i.e., the ensemble size is larger than the available number of processes. " \
                                       "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                       "to the model after this function. " \
                                       "This function should be called at each model time step. \n    \n    " \
                                       "The function executes the user-supplied function " \
                                       "in the following sequence: \n    " \
                                       "1. py__collect_state_pdaf\n    " \
                                       "2. py__prepoststep_state_pdaf\n    " \
                                       "3. py__init_n_domains_p_pdaf\n    " \
                                       "4. py__init_dim_obs_pdaf\n    " \
                                       "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                       "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                       "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                       "loop over each local domain:\n    " \
                                       "8. py__init_dim_l_pdaf\n    " \
                                       "9. py__init_dim_obs_l_pdaf\n    " \
                                       "10. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                       "11. py__init_obs_l_pdaf\n    "\
                                       "12. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                       "13. py__prodRinvA_pdaf\n    " \
                                       "14. py__likelihood_l_pdaf\n    " \
                                       "15. core DA algorithm\n    " \
                                       "16. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                       "17. py__likelihood_hyb_l_pda\n    " \
                                       "18. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                       "19. py__prodRinvA_hyb_l_pdaf\n    "
docstrings['local_put_state_lnetf'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                      "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`. \n    " \
                                      "This function will use Local Nonlinear Ensemble Transform Filter (LNETF) " \
                                      "for DA without post-processing and analysis distribution to forecsat without OMI. The nonlinear filter computes the distribution up to " \
                                      "the second moment similar to KF but using a nonlinear weighting similar to " \
                                      "particle filter. This leads to an equal weights assumption for prior ensemble. \n    \n    " \
                                      "This function is usually used in 'flexible' parallelisation. " \
                                      "i.e., the ensemble size is larger than the available number of processes. " \
                                      "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                      "to the model after this function. " \
                                      "This function should be called at each model time step. \n    \n    " \
                                      "The function executes the user-supplied function " \
                                      "in the following sequence: \n    " \
                                      "1. py__collect_state_pdaf\n    " \
                                      "2. py__prepoststep_state_pdaf\n    " \
                                      "3. py__init_n_domains_p_pdaf\n    " \
                                      "4. py__init_dim_obs_pdaf\n    " \
                                      "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                      "loop over each local domain:\n    " \
                                      "6. py__init_dim_l_pdaf\n    " \
                                      "7. py__init_dim_obs_l_pdaf\n    " \
                                      "8. py__init_obs_l_pdaf\n    "\
                                      "9. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                      "10. py__likelihood_l_pdaf\n    " \
                                      "11. core DA algorithm\n    "
docstrings['local_put_state_lseik'] = "It is recommended to use OMI functionalities for fewer user-supplied functions and improved efficiency. " \
                                      "I.e., `pyPDAF.PDAF.localomi_put_state` or `pyPDAF.PDAF.localomi_put_state_nondiagR`. \n    " \
                                      "Using local singular evolutive interpolated Kalman filter for DA without post-processing and analysis distribution to forecsat without OMI. " \
                                      "This is a domain localisation method. " \
                                      "This function is usually used in 'flexible' parallelisation. " \
                                      "i.e., the ensemble size is larger than the available number of processes. " \
                                      "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                      "to the model after this function. " \
                                      "This function should be called at each model time step. \n    \n    " \
                                      "The function executes the user-supplied function " \
                                      "in the following sequence: \n    " \
                                      "1. py__collect_state_pdaf\n    " \
                                      "2. py__prepoststep_state_pdaf\n    " \
                                      "3. py__init_n_domains_p_pdaf\n    " \
                                      "4. py__init_dim_obs_pdaf\n    " \
                                      "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                      "6. py__init_obs_pdaf (if global adaptive forgetting factor is used (type_forget=1 in pyPDAF.PDAF.init\n    " \
                                      "7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used\n    " \
                                      "loop over each local domain:\n    " \
                                      "8. py__init_dim_l_pdaf\n    " \
                                      "9. py__init_dim_obs_l_pdaf\n    " \
                                      "10. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                      "11. py__init_obs_l_pdaf\n    "\
                                      "12. py__g2l_obs_pdaf (localise each ensemble member in observation space\n    " \
                                      "13. py__init_obsvar_l_pdaf (only called if local adaptive forgetting factor (type_forget=2) is used\n    "\
                                      "14. py__prodRinvA_l_pdaf\n    " \
                                      "15. core DA algorithm\n    "

docstrings['localomi_assimilate'] = "Using domain localised filters for DA with diagonal observation error covariance matrix. " \
                                    "This function should be called at each model time step. \n    \n    " \
                                    "The function is a combination of `pyPDAF.PDAF.localomi_put_state` " \
                                    "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                    "in the following sequence: \n    " \
                                    "1. py__collect_state_pdaf\n    " \
                                    "2. py__prepoststep_state_pdaf\n    " \
                                    "3. py__init_n_domains_p_pdaf\n    " \
                                    "4. py__init_dim_obs_pdaf\n    " \
                                    "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                    "loop over each local domain:\n    " \
                                    "6. py__init_dim_l_pdaf\n    " \
                                    "7. py__init_dim_obs_l_pdaf\n    " \
                                    "8. py__init_obs_l_pdaf\n    "\
                                    "9. core DA algorithm\n    " \
                                    "10. py__prepoststep_state_pdaf\n    " \
                                    "11. py__distribute_state_pdaf\n    " \
                                    "12. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_en3dvar_lestkf'] = "Using 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                                   "The background error covariance matrix is estimated by ensemble. " \
                                                   "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                   "An LESTKF is used to generate ensemble perturbations. " \
                                                   "This function should be called at each model time step. \n    \n    " \
                                                   "The function is a combination of `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
                                                   "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                   "in the following sequence: \n    " \
                                                   "1. py__collect_state_pdaf\n    " \
                                                   "2. py__prepoststep_state_pdaf\n    " \
                                                   "3. py__init_dim_obs_pdaf\n    " \
                                                   "4. py__obs_op_pdaf\n    " \
                                                   "Starting the iterative optimisation:\n    " \
                                                   "5. py__cvt_ens_pdaf\n    " \
                                                   "6. py__obs_op_lin_pdaf\n    " \
                                                   "7. py__obs_op_adj_pdaf\n    " \
                                                   "8. py__cvt_adj_ens_pdaf\n    " \
                                                   "9. core DA algorithm\n    " \
                                                   "After the iterations: \n    " \
                                                   "10. py__cvt_ens_pdaf\n    " \
                                                   "Perform LESTKF: \n    " \
                                                   "11. py__init_n_domains_p_pdaf\n    " \
                                                   "12. py__init_dim_obs_pdaf\n    " \
                                                   "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                                   "loop over each local domain:\n    " \
                                                   "14. py__init_dim_l_pdaf\n    " \
                                                   "15. py__init_dim_obs_l_pdaf (localise mean ensemble in observation space)\n    " \
                                                   "16. core DA algorithm\n    " \
                                                   "17. py__prepoststep_state_pdaf\n    " \
                                                   "18. py__distribute_state_pdaf\n    " \
                                                   "19. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_hyb3dvar_lestkf'] =  "Using Hybrid 3DEnVar for DA with diagonal observation error covariance matrix.\n    " \
                                                     "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                     "and a flow-dependent background error covariance estimated from ensemble. " \
                                                     "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                     "LESTKF in this implementation. " \
                                                     "This function should be called at each model time step. \n    \n    " \
                                                     "The function is a combination of `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                                     "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                     "in the following sequence: \n    " \
                                                     "1. py__collect_state_pdaf\n    " \
                                                     "2. py__prepoststep_state_pdaf\n    " \
                                                     "3. py__init_dim_obs_pdaf\n    " \
                                                     "4. py__obs_op_pdaf\n    " \
                                                     "Starting the iterative optimisation:\n    " \
                                                     "5. py__cvt_pdaf\n    " \
                                                     "6. py__cvt_ens_pdaf\n    " \
                                                     "7. py__obs_op_lin_pdaf\n    " \
                                                     "8. py__obs_op_adj_pdaf\n    " \
                                                     "9. py__cvt_adj_pdaf\n    " \
                                                     "10. py__cvt_adj_ens_pdaf\n    " \
                                                     "11. core DA algorithm\n    " \
                                                     "After the iterations: n    " \
                                                     "12. py__cvt_pdaf\n    " \
                                                     "13. py__cvt_ens_pdaf\n    " \
                                                     "Perform LESTKF: \n    " \
                                                     "14. py__init_n_domains_p_pdaf\n    " \
                                                     "15. py__init_dim_obs_pdaf\n    " \
                                                     "16. py__obs_op_pdaf (for each ensemble member\n    " \
                                                     "loop over each local domain:\n    " \
                                                     "17. py__init_dim_l_pdaf\n    " \
                                                     "18. py__init_dim_obs_l_pdaf\n    " \
                                                     "19. py__init_obs_l_pdaf\n    "\
                                                     "20. core DA algorithm\n    " \
                                                     "21. py__prepoststep_state_pdaf\n    " \
                                                     "22. py__distribute_state_pdaf\n    " \
                                                     "23. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_en3dvar_lestkf_nondiagR'] = "Using 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                                            "The background error covariance matrix is estimated by ensemble. " \
                                                            "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                            "An LESTKF is used to generate ensemble perturbations. " \
                                                            "This function should be called at each model time step. \n    \n    " \
                                                            "The function is a combination of `pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`\n    " \
                                                            "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                            "in the following sequence: \n    " \
                                                            "1. py__collect_state_pdaf\n    " \
                                                            "2. py__prepoststep_state_pdaf\n    " \
                                                            "3. py__init_dim_obs_pdaf\n    " \
                                                            "4. py__obs_op_pdaf\n    " \
                                                            "Starting the iterative optimisation:\n    " \
                                                            "5. py__cvt_ens_pdaf\n    " \
                                                            "6. py__obs_op_lin_pdaf\n    " \
                                                            "7. py__prodRinvA_pdaf\n    " \
                                                            "8. py__obs_op_adj_pdaf\n    " \
                                                            "9. py__cvt_adj_ens_pdaf\n    " \
                                                            "10. core DA algorithm\n    " \
                                                            "After the iterations: \n    " \
                                                            "11. py__cvt_ens_pdaf\n    " \
                                                            "Perform LESTKF: \n    " \
                                                            "12. py__init_n_domains_p_pdaf\n    " \
                                                            "13. py__init_dim_obs_pdaf\n    " \
                                                            "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                            "loop over each local domain:\n    " \
                                                            "15. py__init_dim_l_pdaf\n    " \
                                                            "16. py__init_dim_obs_l_pdaf (localise mean ensemble in observation space)\n    " \
                                                            "17. py__prodRinvA_l_pdaf\n    " \
                                                            "18. core DA algorithm\n    " \
                                                            "19. py__prepoststep_state_pdaf\n    " \
                                                            "20. py__distribute_state_pdaf\n    " \
                                                            "21. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_hyb3dvar_lestkf_nondiagR'] = "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix.\n    " \
                                                             "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                             "and a flow-dependent background error covariance estimated from ensemble. " \
                                                             "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                             "LESTKF in this implementation. " \
                                                             "This function should be called at each model time step. \n    \n    " \
                                                             "The function is a combination of `pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`\n    " \
                                                             "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                             "in the following sequence: \n    " \
                                                             "1. py__collect_state_pdaf\n    " \
                                                             "2. py__prepoststep_state_pdaf\n    " \
                                                             "3. py__init_dim_obs_pdaf\n    " \
                                                             "4. py__obs_op_pdaf\n    " \
                                                             "Starting the iterative optimisation:\n    " \
                                                             "5. py__cvt_pdaf\n    " \
                                                             "6. py__cvt_ens_pdaf\n    " \
                                                             "7. py__obs_op_lin_pdaf\n    " \
                                                             "8. py__prodRinvA_pdaf\n    " \
                                                             "9. py__obs_op_adj_pdaf\n    " \
                                                             "10. py__cvt_adj_pdaf\n    " \
                                                             "11. py__cvt_adj_ens_pdaf\n    " \
                                                             "12. core DA algorithm\n    " \
                                                             "After the iterations: \n    " \
                                                             "13. py__cvt_pdaf\n    " \
                                                             "14. py__cvt_ens_pdaf\n    " \
                                                             "Perform LESTKF: \n    " \
                                                             "15. py__init_n_domains_p_pdaf\n    " \
                                                             "16. py__init_dim_obs_pdaf\n    " \
                                                             "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                                             "loop over each local domain:\n    " \
                                                             "18. py__init_dim_l_pdaf\n    " \
                                                             "19. py__init_dim_obs_l_pdaf\n    " \
                                                             "20. py__init_obs_l_pdaf\n    "\
                                                             "21. py__prodRinvA_l_pdaf\n    " \
                                                             "22. core DA algorithm\n    " \
                                                             "23. py__prepoststep_state_pdaf\n    " \
                                                             "24. py__distribute_state_pdaf\n    " \
                                                             "25. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_lknetf_nondiagR'] = "Using LKNETF for DA with non-diagonal observation error covariance matrix. " \
                                                    "This function should be called at each model time step. \n    \n    " \
                                                    "The function is a combination of `pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR` " \
                                                    "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                    "in the following sequence: \n    " \
                                                    "1. py__collect_state_pdaf\n    " \
                                                    "2. py__prepoststep_state_pdaf\n    " \
                                                    "3. py__init_n_domains_p_pdaf\n    " \
                                                    "4. py__init_dim_obs_pdaf\n    " \
                                                    "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                                    "loop over each local domain:\n    " \
                                                    "6. py__init_dim_l_pdaf\n    " \
                                                    "7. py__init_dim_obs_l_pdaf\n    " \
                                                    "8. py__prodRinvA_pdaf\n    " \
                                                    "9. py__likelihood_l_pdaf\n    " \
                                                    "10. core DA algorithm\n    " \
                                                    "11. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                                    "12. py__likelihood_hyb_l_pda\n    " \
                                                    "13. py__prodRinvA_hyb_l_pdaf\n    " \
                                                    "14. py__prepoststep_state_pdaf\n    " \
                                                    "15. py__distribute_state_pdaf\n    " \
                                                    "16. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_lnetf_nondiagR'] = "Using LNETF for DA with non-diagonal observation error covariance matrix. " \
                                                   "This function should be called at each model time step. \n    \n    " \
                                                   "The function is a combination of `pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR` " \
                                                   "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                                   "in the following sequence: \n    " \
                                                   "1. py__collect_state_pdaf\n    " \
                                                   "2. py__prepoststep_state_pdaf\n    " \
                                                   "3. py__init_n_domains_p_pdaf\n    " \
                                                   "4. py__init_dim_obs_pdaf\n    " \
                                                   "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                                   "loop over each local domain:\n    " \
                                                   "6. py__init_dim_l_pdaf\n    " \
                                                   "7. py__init_dim_obs_l_pdaf\n    " \
                                                   "8. py__likelihood_l_pdaf\n    " \
                                                   "9. core DA algorithm\n    " \
                                                   "10. py__prepoststep_state_pdaf\n    " \
                                                   "11. py__distribute_state_pdaf\n    " \
                                                   "12. py__next_observation_pdaf\n    "
docstrings['localomi_assimilate_nondiagR'] = "Using domain localised filters for DA with non-diagonal observation error covariance matrix. " \
                                             "This function should be called at each model time step. \n    \n    " \
                                             "The function is a combination of `pyPDAF.PDAF.localomi_put_state_local_nondiagR` " \
                                             "and `pyPDAF.PDAF.get_state`, and executes the user-supplied function " \
                                             "in the following sequence: \n    " \
                                             "1. py__collect_state_pdaf\n    " \
                                             "2. py__prepoststep_state_pdaf\n    " \
                                             "3. py__init_n_domains_p_pdaf\n    " \
                                             "4. py__init_dim_obs_pdaf\n    " \
                                             "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                             "loop over each local domain:\n    " \
                                             "6. py__init_dim_l_pdaf\n    " \
                                             "7. py__init_dim_obs_l_pdaf\n    " \
                                             "8. py__init_obs_l_pdaf\n    "\
                                             "9. py__prodRinvA_l_pdaf\n    " \
                                             "10. core DA algorithm\n    " \
                                             "11. py__prepoststep_state_pdaf\n    " \
                                             "12. py__distribute_state_pdaf\n    " \
                                             "13. py__next_observation_pdaf\n    "

docstrings['localomi_put_state'] = "Using domain localised filters for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix. " \
                                   "This is a domain localisation method. " \
                                   "This function is usually used in 'flexible' parallelisation. " \
                                   "i.e., the ensemble size is larger than the available number of processes. " \
                                   "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                   "to the model after this function. " \
                                   "This function should be called at each model time step. " \
                                   "The LESTKF is a more efficient equivalent to the LETKF. \n    \n    " \
                                   "The function executes the user-supplied function " \
                                   "in the following sequence: \n    " \
                                   "1. py__collect_state_pdaf\n    " \
                                   "2. py__prepoststep_state_pdaf\n    " \
                                   "3. py__init_n_domains_p_pdaf\n    " \
                                   "4. py__init_dim_obs_pdaf\n    " \
                                   "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                   "loop over each local domain:\n    " \
                                   "6. py__init_dim_l_pdaf\n    " \
                                   "7. py__init_dim_obs_l_pdaf\n    " \
                                   "8. core DA algorithm\n    " \
                                   "9. py__prepoststep_state_pdaf\n    " \
                                   "10. py__distribute_state_pdaf\n    " \
                                   "11. py__next_observation_pdaf\n    "
docstrings['localomi_put_state_en3dvar_lestkf'] = "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                                  "The background error covariance matrix is estimated by ensemble. " \
                                                  "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                  "An LESTKF is used to generate ensemble perturbations. " \
                                                  "This function is usually used in 'flexible' parallelisation. " \
                                                  "i.e., the ensemble size is larger than the available number of processes. " \
                                                  "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                  "to the model after this function. " \
                                                  "This function should be called at each model time step. \n    \n    " \
                                                  "The function executes the user-supplied function " \
                                                  "in the following sequence: \n    " \
                                                  "1. py__collect_state_pdaf\n    " \
                                                  "2. py__prepoststep_state_pdaf\n    " \
                                                  "3. py__init_dim_obs_pdaf\n    " \
                                                  "4. py__obs_op_pdaf\n    " \
                                                  "Starting the iterative optimisation:\n    " \
                                                  "5. py__cvt_ens_pdaf\n    " \
                                                  "6. py__obs_op_lin_pdaf\n    " \
                                                  "7. py__obs_op_adj_pdaf\n    " \
                                                  "8. py__cvt_adj_ens_pdaf\n    " \
                                                  "9. core DA algorithm\n    " \
                                                  "After the iterations: \n    " \
                                                  "10. py__cvt_ens_pdaf\n    " \
                                                  "Perform LESTKF: \n    " \
                                                  "11. py__init_n_domains_p_pdaf\n    " \
                                                  "12. py__init_dim_obs_pdaf\n    " \
                                                  "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                                  "loop over each local domain:\n    " \
                                                  "14. py__init_dim_l_pdaf\n    " \
                                                  "15. py__init_dim_obs_l_pdaf\n    " \
                                                  "16. core DA algorithm\n    "
docstrings['localomi_put_state_hyb3dvar_lestkf'] = "Using 3DEnVar for DA without post-processing and analysis distribution to forecsat with diagonal observation error covariance matrix.\n    " \
                                                   "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                   "and a flow-dependent background error covariance estimated from ensemble. " \
                                                   "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                   "LESTKF in this implementation. \n    \n    " \
                                                   "This function is usually used in 'flexible' parallelisation. " \
                                                   "i.e., the ensemble size is larger than the available number of processes. " \
                                                   "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                   "to the model after this function. " \
                                                   "This function should be called at each model time step. \n    \n    " \
                                                   "The function executes the user-supplied function " \
                                                   "in the following sequence: \n    " \
                                                   "1. py__collect_state_pdaf\n    " \
                                                   "2. py__prepoststep_state_pdaf\n    " \
                                                   "3. py__init_dim_obs_pdaf\n    " \
                                                   "4. py__obs_op_pdaf\n    " \
                                                   "Starting the iterative optimisation:\n    " \
                                                   "5. py__cvt_ens_pdaf\n    " \
                                                   "6. py__obs_op_lin_pdaf\n    " \
                                                   "7. py__obs_op_adj_pdaf\n    " \
                                                   "8. py__cvt_adj_ens_pdaf\n    " \
                                                   "9. core DA algorithm\n    " \
                                                   "After the iterations: \n    " \
                                                   "10. py__cvt_ens_pdaf\n    " \
                                                   "Perform LESTKF: \n    " \
                                                   "11. py__init_n_domains_p_pdaf\n    " \
                                                   "12. py__init_dim_obs_pdaf\n    " \
                                                   "13. py__obs_op_pdaf (for each ensemble member\n    " \
                                                   "loop over each local domain:\n    " \
                                                   "14. py__init_dim_l_pdaf\n    " \
                                                   "15. py__init_dim_obs_l_pdaf (localise mean ensemble in observation space)\n    " \
                                                   "16. core DA algorith\n    "
docstrings['localomi_put_state_en3dvar_lestkf_nondiagR'] = "Using 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                           "without post-processing and analysis distribution to forecsat without OMI. " \
                                                           "The background error covariance matrix is estimated by ensemble. " \
                                                           "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                           "An LESTKF is used to generate ensemble perturbations. " \
                                                           "This function is usually used in 'flexible' parallelisation. " \
                                                           "i.e., the ensemble size is larger than the available number of processes. " \
                                                           "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                           "to the model after this function. " \
                                                           "This function should be called at each model time step. \n    \n    " \
                                                           "The function executes the user-supplied function " \
                                                           "in the following sequence: \n    " \
                                                           "1. py__collect_state_pdaf\n    " \
                                                           "2. py__prepoststep_state_pdaf\n    " \
                                                           "3. py__init_dim_obs_pdaf\n    " \
                                                           "4. py__obs_op_pdaf\n    " \
                                                           "Starting the iterative optimisation:\n    " \
                                                           "5. py__cvt_ens_pdaf\n    " \
                                                           "6. py__obs_op_lin_pdaf\n    " \
                                                           "7. py__prodRinvA_pdaf\n    " \
                                                           "8. py__obs_op_adj_pdaf\n    " \
                                                           "9. py__cvt_adj_ens_pdaf\n    " \
                                                           "10. core DA algorithm\n    " \
                                                           "After the iterations: \n    " \
                                                           "11. py__cvt_ens_pdaf\n    " \
                                                           "Perform LESTKF: \n    " \
                                                           "12. py__init_n_domains_p_pdaf\n    " \
                                                           "13. py__init_dim_obs_pdaf\n    " \
                                                           "14. py__obs_op_pdaf (for each ensemble member\n    " \
                                                           "loop over each local domain:\n    " \
                                                           "15. py__init_dim_l_pdaf\n    " \
                                                           "16. py__init_dim_obs_l_pdaf (localise mean ensemble in observation space)\n    " \
                                                           "17. py__prodRinvA_l_pdaf\n    " \
                                                           "18. core DA algorithm\n    "
docstrings['localomi_put_state_hyb3dvar_lestkf_nondiagR'] = "Using Hybrid 3DEnVar for DA with non-diagonal observation error covariance matrix\n    " \
                                                            "without post-processing and analysis distribution to forecsat without OMI. " \
                                                            "Here, the background error covariance is hybridised by a static background error covariance, " \
                                                            "and a flow-dependent background error covariance estimated from ensemble. " \
                                                            "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                            "LESTKF in this implementation. " \
                                                            "This function should be called at each model time step. \n    \n    " \
                                                            "This function is usually used in 'flexible' parallelisation. " \
                                                            "i.e., the ensemble size is larger than the available number of processes. " \
                                                            "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                            "to the model after this function. " \
                                                            "This function should be called at each model time step. \n    \n    " \
                                                            "The function executes the user-supplied function " \
                                                            "in the following sequence: \n    " \
                                                            "1. py__collect_state_pdaf\n    " \
                                                            "2. py__prepoststep_state_pdaf\n    " \
                                                            "3. py__init_dim_obs_pdaf\n    " \
                                                            "4. py__obs_op_pdaf\n    " \
                                                            "Starting the iterative optimisation:\n    " \
                                                            "5. py__cvt_pdaf\n    " \
                                                            "6. py__cvt_ens_pdaf\n    " \
                                                            "7. py__obs_op_lin_pdaf\n    " \
                                                            "8. py__prodRinvA_pdaf\n    " \
                                                            "9. py__obs_op_adj_pdaf\n    " \
                                                            "10. py__cvt_adj_pdaf\n    " \
                                                            "11. py__cvt_adj_ens_pdaf\n    " \
                                                            "12. core DA algorithm\n    " \
                                                            "After the iterations: \n    " \
                                                            "13. py__cvt_pdaf\n    " \
                                                            "14. py__cvt_ens_pdaf\n    " \
                                                            "Perform LESTKF: \n    " \
                                                            "15. py__init_n_domains_p_pdaf\n    " \
                                                            "16. py__init_dim_obs_pdaf\n    " \
                                                            "17. py__obs_op_pdaf (for each ensemble member\n    " \
                                                            "loop over each local domain:\n    " \
                                                            "18. py__init_dim_l_pdaf\n    " \
                                                            "19. py__init_dim_obs_l_pdaf\n    " \
                                                            "20. py__init_obs_l_pdaf\n    "\
                                                            "21. py__prodRinvA_l_pdaf\n    " \
                                                            "22. core DA algorithm\n    "
docstrings['localomi_put_state_lknetf_nondiagR'] = "Using LKNETF for DA with non-diagonal observation error covariance matrix " \
                                                   "without post-processing and analysis distribution to forecsat without OMI. " \
                                                   "This function should be called at each model time step. \n    \n    " \
                                                   "This function is usually used in 'flexible' parallelisation. " \
                                                   "i.e., the ensemble size is larger than the available number of processes. " \
                                                   "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                   "to the model after this function. " \
                                                   "This function should be called at each model time step. \n    \n    " \
                                                   "The function executes the user-supplied function " \
                                                   "in the following sequence: \n    " \
                                                   "1. py__collect_state_pdaf\n    " \
                                                   "2. py__prepoststep_state_pdaf\n    " \
                                                   "3. py__init_n_domains_p_pdaf\n    " \
                                                   "4. py__init_dim_obs_pdaf\n    " \
                                                   "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                                   "loop over each local domain:\n    " \
                                                   "6. py__init_dim_l_pdaf\n    " \
                                                   "7. py__init_dim_obs_l_pdaf\n    " \
                                                   "8. py__prodRinvA_pdaf\n    " \
                                                   "9. py__likelihood_l_pdaf\n    " \
                                                   "10. core DA algorithm\n    " \
                                                   "11. py__obs_op_pdaf (only called with `HKN` and `HNK` options called for each ensemble member\n    " \
                                                   "12. py__likelihood_hyb_l_pda\n    " \
                                                   "13. py__prodRinvA_hyb_l_pdaf\n    "
docstrings['localomi_put_state_lnetf_nondiagR'] = "Using LNETF for DA with non-diagonal observation error covariance matrix " \
                                                  "without post-processing and analysis distribution to forecsat without OMI. " \
                                                  "This function should be called at each model time step. \n    \n    " \
                                                  "The function is a combination of `pyPDAF.PDAF.omi_put_state_lnetf_nondiagR` " \
                                                  "This function is usually used in 'flexible' parallelisation. " \
                                                  "i.e., the ensemble size is larger than the available number of processes. " \
                                                  "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                                  "to the model after this function. " \
                                                  "This function should be called at each model time step. \n    \n    " \
                                                  "The function executes the user-supplied function " \
                                                  "in the following sequence: \n    " \
                                                  "1. py__collect_state_pdaf\n    " \
                                                  "2. py__prepoststep_state_pdaf\n    " \
                                                  "3. py__init_n_domains_p_pdaf\n    " \
                                                  "4. py__init_dim_obs_pdaf\n    " \
                                                  "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                                  "loop over each local domain:\n    " \
                                                  "6. py__init_dim_l_pdaf\n    " \
                                                  "7. py__init_dim_obs_l_pdaf\n    " \
                                                  "8. py__likelihood_l_pdaf\n    " \
                                                  "9. core DA algorithm\n    "
docstrings['localomi_put_state_nondiagR'] = "Using domain localised filters for DA with non-diagonal observation error covariance matrix " \
                                            "without post-processing and analysis distribution to forecsat without OMI. " \
                                            "This function should be called at each model time step. \n    \n    " \
                                            "This function is usually used in 'flexible' parallelisation. " \
                                            "i.e., the ensemble size is larger than the available number of processes. " \
                                            "A `pyPDAF.PDAF.get_state` function should be used to post-process and distribute the ensemble " \
                                            "to the model after this function. " \
                                            "This function should be called at each model time step. \n    \n    " \
                                            "The function executes the user-supplied function " \
                                            "in the following sequence: \n    " \
                                            "1. py__collect_state_pdaf\n    " \
                                            "2. py__prepoststep_state_pdaf\n    " \
                                            "3. py__init_n_domains_p_pdaf\n    " \
                                            "4. py__init_dim_obs_pdaf\n    " \
                                            "5. py__obs_op_pdaf (for each ensemble member\n    " \
                                            "loop over each local domain:\n    " \
                                            "6. py__init_dim_l_pdaf\n    " \
                                            "7. py__init_dim_obs_l_pdaf\n    " \
                                            "9. py__init_obs_l_pdaf\n    "\
                                            "10. py__prodRinvA_l_pdaf\n    " \
                                            "11. core DA algorithm\n    "
