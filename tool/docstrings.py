docstrings = dict()
docstrings['assimilate_3dvar'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_3dvar`\n    "\
                                 "or :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "3DVar DA for a single step without OMI.\n    " \
                                 "When 3DVar is used, the background error covariance matrix\n    "\
                                 "has to be modelled for cotrol variable transformation.\n    " \
                                 "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                 "This function should be called at each model time step.\n\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_3dvar`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                 "User-supplied functions are executed in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pdaf\n    " \
                                 "    5. py__init_obs_pdaf\n    " \
                                 "    6. Iterative optimisation:\n    " \
                                 "        1. py__cvt_pdaf\n    " \
                                 "        2. py__obs_op_lin_pdaf\n    " \
                                 "        3. py__prodRinvA_pdaf\n    " \
                                 "        4. py__obs_op_adj_pdaf\n    " \
                                 "        5. py__cvt_adj_pdaf\n    " \
                                 "        6. core DA algorithm\n    " \
                                 "    7. py__cvt_pdaf\n    " \
                                 "    8. py__prepoststep_state_pdaf\n    " \
                                 "    9. py__distribute_state_pdaf\n    " \
                                 "    10. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_3dvar`\n    " \
                                 "   and :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`"
docstrings['assimilate_en3dvar_estkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`\n    "\
                                         "or :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`.\n\n    "\
                                         "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                         "3DEnVar for a single DA step.\n    " \
                                         "The background error covariance matrix is estimated by an ensemble.\n    " \
                                         "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                         "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                         "This function should be called at each model time step.\n\n    " \
                                         "The function is a combination of :func:`pyPDAF.PDAF.put_state_en3dvar_estkf`\n    " \
                                         "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                         "The user-supplied functions are executed in the following sequence:\n    " \
                                         "    1. py__collect_state_pdaf\n    " \
                                         "    2. py__prepoststep_state_pdaf\n    " \
                                         "    3. py__init_dim_obs_pdaf\n    " \
                                         "    4. py__obs_op_pdaf\n    " \
                                         "    5. py__init_obs_pdaf\n    " \
                                         "    6. the iterative optimisation:\n    " \
                                         "        1. py__cvt_ens_pdaf\n    " \
                                         "        2. py__obs_op_lin_pdaf\n    " \
                                         "        3. py__prodRinvA_pdaf\n    " \
                                         "        4. py__obs_op_adj_pdaf\n    " \
                                         "        5. py__cvt_adj_ens_pdaf\n    " \
                                         "        6. core 3DEnVar algorithm\n    " \
                                         "    7. py__cvt_ens_pdaf\n    " \
                                         "    8. ESTKF:\n    " \
                                         "        1. py__init_dim_obs_pdaf\n    " \
                                         "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                         "        3. py__init_obs_pdaf\n    " \
                                         "        4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                         "        5. py__init_obsvar_pdaf\n    " \
                                         "           (only relevant for adaptive forgetting factor schemes)\n    " \
                                         "        6. py__prodRinvA_pdaf\n    " \
                                         "        7. core ESTKF algorithm\n    " \
                                         "    9. py__prepoststep_state_pdaf\n    " \
                                         "    10. py__distribute_state_pdaf\n    " \
                                         "    11. py__next_observation_pdaf\n" \
                                         "\n    " \
                                         ".. deprecated:: 1.0.0\n\n    " \
                                         "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`\n    " \
                                         "   and :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`"
docstrings['assimilate_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    "\
                                          "or :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`.\n\n    "\
                                          "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                          "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF.\n    " \
                                          "The background error covariance matrix is estimated by ensemble.\n    " \
                                          "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                          "An LESTKF is used to generate ensemble perturbations.\n    " \
                                          "This function should be called at each model time step.\n\n    " \
                                          "The function is a combination of :func:`pyPDAF.PDAF.put_state_en3dvar_lestkf`\n    " \
                                          "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                          "The user-supplied function are executed in the following sequence:\n    " \
                                          "    1. py__collect_state_pdaf\n    " \
                                          "    2. py__prepoststep_state_pdaf\n    " \
                                          "    3. py__init_dim_obs_pdaf\n    " \
                                          "    4. py__obs_op_pdaf\n    " \
                                          "    5. py__init_obs_pdaf\n    " \
                                          "    6. Starting the iterative optimisation:\n    " \
                                          "        1. py__cvt_ens_pdaf\n    " \
                                          "        2. py__obs_op_lin_pdaf\n    " \
                                          "        3. py__prodRinvA_pdaf\n    " \
                                          "        4. py__obs_op_adj_pdaf\n    " \
                                          "        5. py__cvt_adj_ens_pdaf\n    " \
                                          "        6. core DA algorithm\n    " \
                                          "    7. py__cvt_ens_pdaf\n    " \
                                          "    8. Perform LESTKF:\n    " \
                                          "        1. py__init_n_domains_p_pdaf\n    " \
                                          "        2. py__init_dim_obs_pdaf\n    " \
                                          "        3. py__obs_op_pdaf\n    "\
                                          "           (for each ensemble member)\n    " \
                                          "        4. py__init_obs_pdaf\n    "\
                                          "           (if global adaptive forgetting factor is used\n    "\
                                          "           `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                          "        5. py__init_obsvar_pdaf\n    "\
                                          "           (if global adaptive forgetting factor is used)\n    " \
                                          "        6. loop over each local domain:\n    " \
                                          "            1. py__init_dim_l_pdaf\n    " \
                                          "            2. py__init_dim_obs_l_pdaf\n    " \
                                          "            3. py__g2l_state_pdaf\n    " \
                                          "            4. py__g2l_obs_pdaf\n    "\
                                          "               (localise mean ensemble in observation space)\n    " \
                                          "            5. py__init_obs_l_pdaf\n    "\
                                          "            6. py__g2l_obs_pdaf\n    " \
                                          "               (localise each ensemble member in observation space)\n    " \
                                          "            7. py__init_obsvar_l_pdaf\n    " \
                                          "               (only called if local adaptive forgetting factor\n    " \
                                          "               `type_forget=2` is used)\n    "\
                                          "            8. py__prodRinvA_l_pdaf\n    " \
                                          "            9. core DA algorithm\n    " \
                                          "            10. py__l2g_state_pdaf\n    " \
                                          "    9. py__prepoststep_state_pdaf\n    " \
                                          "    10. py__distribute_state_pdaf\n    " \
                                          "    11. py__next_observation_pdaf\n" \
                                          "\n    " \
                                          ".. deprecated:: 1.0.0\n\n    " \
                                          "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    " \
                                          "   and :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`"
docstrings['assimilate_hyb3dvar_estkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf`\n    "\
                                          "or :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`.\n\n    "\
                                          "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                          "Hybrid 3DEnVar for a single DA step where\n    " \
                                          "the background error covariance is hybridised by a static background error covariance,\n    " \
                                          "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                          "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                          "ESTKF in this implementation.\n    " \
                                          "This function should be called at each model time step.\n\n    " \
                                          "The function is a combination of :func:`pyPDAF.PDAF.put_state_hyb3dvar_estkf`\n    " \
                                          "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                          "The user-supplied functions are executed in the following sequence:\n    " \
                                          "    1. py__collect_state_pdaf\n    " \
                                          "    2. py__prepoststep_state_pdaf\n    " \
                                          "    3. py__init_dim_obs_pdaf\n    " \
                                          "    4. py__obs_op_pdaf\n    " \
                                          "    5. py__init_obs_pdaf\n    " \
                                          "    6. the iterative optimisation:\n    " \
                                          "        1. py__cvt_pdaf\n    " \
                                          "        2. py__cvt_ens_pdaf\n    " \
                                          "        3. py__obs_op_lin_pdaf\n    " \
                                          "        4. py__prodRinvA_pdaf\n    " \
                                          "        5. py__obs_op_adj_pdaf\n    " \
                                          "        6. py__cvt_adj_pdaf\n    " \
                                          "        7. py__cvt_adj_ens_pdaf\n    " \
                                          "        8. core 3DEnVar algorithm\n    " \
                                          "    7. py__cvt_pdaf\n    " \
                                          "    8. py__cvt_ens_pdaf\n    " \
                                          "    9. Perform ESTKF:\n    " \
                                          "        1. py__init_dim_obs_pdaf\n    " \
                                          "        2. py__obs_op_pdaf\n    "\
                                          "           (for ensemble mean)\n    " \
                                          "        3. py__init_obs_pdaf\n    " \
                                          "        4. py__obs_op_pdaf\n    " \
                                          "           (for each ensemble member)\n    " \
                                          "        5. py__init_obsvar_pdaf\n    " \
                                          "           (only relevant for adaptive forgetting factor schemes)\n    " \
                                          "        6. py__prodRinvA_pdaf\n    " \
                                          "        7. core ESTKF algorithm\n    " \
                                          "    10. py__prepoststep_state_pdaf\n    " \
                                          "    11. py__distribute_state_pdaf\n    " \
                                          "    12. py__next_observation_pdaf\n" \
                                          "\n    " \
                                          ".. deprecated:: 1.0.0\n\n    " \
                                          "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf`\n    " \
                                          "   and :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`"
docstrings['assimilate_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    "\
                                           "or :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                           "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                           "Hybrid 3DEnVar for a single DA step where\n    " \
                                           "the background error covariance is hybridised by a static background error covariance,\n    " \
                                           "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                           "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                           "LESTKF in this implementation.\n    " \
                                           "This function should be called at each model time step.\n\n    " \
                                           "The function is a combination of :func:`pyPDAF.PDAF.put_state_hyb3dvar_lestkf`\n    " \
                                           "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                           "The user-supplied functions are executed in the following sequence:\n    " \
                                           "    1. py__collect_state_pdaf\n    " \
                                           "    2. py__prepoststep_state_pdaf\n    " \
                                           "    3. py__init_dim_obs_pdaf\n    " \
                                           "    4. py__obs_op_pdaf\n    " \
                                           "    5. py__init_obs_pdaf\n    " \
                                           "    6. The iterative optimisation:\n    " \
                                           "        1. py__cvt_pdaf\n    " \
                                           "        2. py__cvt_ens_pdaf\n    " \
                                           "        3. py__obs_op_lin_pdaf\n    " \
                                           "        4. py__prodRinvA_pdaf\n    " \
                                           "        5. py__obs_op_adj_pdaf\n    " \
                                           "        6. py__cvt_adj_pdaf\n    " \
                                           "        7. py__cvt_adj_ens_pdaf\n    " \
                                           "        8. core DA algorithm\n    " \
                                           "    7. py__cvt_pdaf\n    " \
                                           "    8. py__cvt_ens_pdaf\n    " \
                                           "    9. Perform LESTKF:\n    " \
                                           "        1. py__init_n_domains_p_pdaf\n    " \
                                           "        2. py__init_dim_obs_pdaf\n    " \
                                           "        3. py__obs_op_pdaf\n    " \
                                           "           (for each ensemble member)\n    " \
                                           "        4. py__init_obs_pdaf\n    " \
                                           "           (if global adaptive forgetting factor `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                           "        5. py__init_obsvar_pdaf\n    " \
                                           "           (if global adaptive forgetting factor is used)\n    " \
                                           "        6. loop over each local domain:\n    " \
                                           "            1. py__init_dim_l_pdaf\n    " \
                                           "            2. py__init_dim_obs_l_pdaf\n    " \
                                           "            3. py__g2l_state_pdaf\n    " \
                                           "            4. py__g2l_obs_pdaf\n    "\
                                           "               (localise mean ensemble in observation space)\n    " \
                                           "            5. py__init_obs_l_pdaf\n    "\
                                           "            6. py__g2l_obs_pdaf\n    " \
                                           "               (localise each ensemble member in observation space)\n    " \
                                           "            7. py__init_obsvar_l_pdaf\n    " \
                                           "               (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                           "            8. py__prodRinvA_l_pdaf\n    " \
                                           "            9. core DA algorithm\n    " \
                                           "            10. py__l2g_state_pdaf\n    " \
                                           "    10. py__prepoststep_state_pdaf\n    " \
                                           "    11. py__distribute_state_pdaf\n    " \
                                           "    12. py__next_observation_pdaf\n" \
                                           "\n    " \
                                           ".. deprecated:: 1.0.0\n\n    " \
                                           "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    " \
                                           "   and :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`"
docstrings['assimilate_enkf'] =  "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                                 "or :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Stochastic EnKF (ensemble " \
                                 "Kalman filter) [1]_ for a single DA step without OMI. " \
                                 "This function should be called at each model time step. \n\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_enkf` " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                 "This function executes the user-supplied functions in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                 "    5. py__add_obs_err_pdaf\n    " \
                                 "    6. py__init_obs_pdaf\n    " \
                                 "    7. py__init_obscovar_pdaf\n    " \
                                 "    8. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    9. core DA algorithm\n    " \
                                 "    10. py__prepoststep_state_pdaf\n    " \
                                 "    11. py__distribute_state_pdaf\n    " \
                                 "    12. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                 "   and :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Evensen, G. (1994), \n    "\
                                 "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                 "       using Monte Carlo methods to forecast error statistics,\n    "\
                                 "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572."
docstrings['assimilate_estkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                 "or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR` instead of this function.\n\n    " \
                                 "OMI functions need fewer user-supplied functions and improve DA efficiency.\n\n    " \
                                 "This function calls ESTKF (error space transform Kalman filter) [1]_.\n    " \
                                 "The ESTKF is a more efficient equivalent to the ETKF.\n\n    " \
                                 "The function should be called at each model time step.\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_estkf`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                 "User-supplied functions are executed in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                 "    5. py__init_obs_pdaf\n    " \
                                 "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                 "    8. py__prodRinvA_pdaf\n    " \
                                 "    9. core DA algorithm\n    " \
                                 "    10. py__prepoststep_state_pdaf\n    " \
                                 "    11. py__distribute_state_pdaf\n    " \
                                 "    12. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                 "   and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                 "       A unification of ensemble square root Kalman filters. \n    " \
                                 "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['assimilate_etkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`.\n\n    "\
                                "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "Using ETKF (ensemble transform " \
                                "Kalman filter) [1]_ for a single DA step without OMI. The implementation is baed on [2]_.\n\n    " \
                                "This function should be called at each model time step.\n    " \
                                "The function is a combination of :func:`pyPDAF.PDAF.put_state_etkf` " \
                                "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                "This function executes the user-supplied function in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                "    5. py__init_obs_pdaf\n    " \
                                "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                "    8. py__prodRinvA_pdaf\n    " \
                                "    9. core DA algorithm\n    " \
                                "    10. py__prepoststep_state_pdaf\n    " \
                                "    11. py__distribute_state_pdaf\n    " \
                                "    12. py__next_observation_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Bishop, C. H., B. J. Etherton, and S. J. Majumdar (2001)\n    "\
                                "       Adaptive Sampling with the Ensemble Transform Kalman Filter.\n    "\
                                "       Part I: Theoretical Aspects. Mon. Wea. Rev., 129, 420–436,\n    "\
                                "       doi: 10.1175/1520-0493(2001)129<0420:ASWTET>2.0.CO;2. \n    " \
                                ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                "       A unification of ensemble square root Kalman filters. \n    " \
                                "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['assimilate_seek'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "This function will use singular evolutive extended Kalman filter [1]_ for a single DA step.\n    " \
                                "This is a deterministic Kalman filter.\n    " \
                                "The function should be called at each model step.\n\n    " \
                                "The function is a combination of :func:`pyPDAF.PDAF.put_state_seek`\n    " \
                                "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                "This function executes the user-supplied functions in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                "    5. py__init_obs_pdaf\n    " \
                                "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    7. py__prodRinvA_pdaf\n    " \
                                "    8. core DA algorithm\n    " \
                                "    9. py__prepoststep_state_pdaf\n    " \
                                "    10. py__distribute_state_pdaf\n    " \
                                "    11. py__next_observation_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                 ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                 "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                 "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['assimilate_seik'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "This function will use singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n    " \
                                "The function should be called at each model step.\n\n    " \
                                "The function is a combination of :func:`pyPDAF.PDAF.put_state_seik`\n    " \
                                "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                "The function executes the user-supplied functions in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                "    5. py__init_obs_pdaf\n    " \
                                "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                "    8. py__prodRinvA_pdaf\n    " \
                                "    9. core DA algorithm\n    " \
                                "    10. py__prepoststep_state_pdaf\n    " \
                                "    11. py__distribute_state_pdaf\n    " \
                                "    12. py__next_observation_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                 ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                 "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                 "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['assimilate_netf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "This function will use Nonlinear Ensemble Transform Filter (NETF) [1]_ \n    " \
                                "for a single DA step. The nonlinear filter computes the distribution up to\n    " \
                                "the second moment similar to KF but using a nonlinear weighting similar to\n    " \
                                "particle filter. This leads to an equal weights assumption for prior ensemble.\n    " \
                                "The function should be called at each model step.\n\n    " \
                                "The function is a combination of :func:`pyPDAF.PDAF.put_state_netf`\n    " \
                                "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                "This function executes the user-supplied function in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__init_obs_pdaf\n    " \
                                "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    6. py__likelihood_pdaf\n    " \
                                "    7. core DA algorithm\n    " \
                                "    8. py__prepoststep_state_pdaf\n    " \
                                "    9. py__distribute_state_pdaf\n    " \
                                "    10. py__next_observation_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                "       A second-order exact ensemble square root filter\n    " \
                                "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['assimilate_pf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_global`\n    "\
                              "or :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`.\n\n    "\
                              "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                              "This function will use particle filter for a single DA step.\n    " \
                              "This is a fully nonlinear filter, and may require a high number of ensemble members.\n    " \
                              "A review of particle filter can be found at [1]_.\n    " \
                              "The function should be called at each model step.\n\n    " \
                              "The function is a combination of :func:`pyPDAF.PDAF.put_state_pf`\n    " \
                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                              "This function executes the user-supplied functions in the following sequence:\n    " \
                              "    1. py__collect_state_pdaf\n    " \
                              "    2. py__prepoststep_state_pdaf\n    " \
                              "    3. py__init_dim_obs_pdaf\n    " \
                              "    4. py__init_obs_pdaf\n    " \
                              "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                              "    6. py__likelihood_pdaf\n    " \
                              "    7. core DA algorithm\n    " \
                              "    8. py__prepoststep_state_pdaf\n    " \
                              "    9. py__distribute_state_pdaf\n    " \
                              "    10. py__next_observation_pdaf\n" \
                              "\n    " \
                              ".. deprecated:: 1.0.0\n\n    " \
                              "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_global`\n    " \
                              "   and :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`" \
                              "\n\n    " \
                              "References\n    " \
                              "----------\n    " \
                              ".. [1] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                              "       Particle filters for high‐dimensional geoscience applications:\n    "\
                              "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['assimilate_lenkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_assimilate_lenkf`\n    "\
                                 "or :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Stochastic EnKF (ensemble Kalman filter) with covariance localisation [1]_\n    " \
                                 "for a single DA step without OMI.\n\n    " \
                                 "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                 "This function should be called at each model time step.\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_lenkf`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                 "The user-supplied function is executed in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    5. py__localize_pdaf\n    " \
                                 "    6. py__add_obs_err_pdaf\n    " \
                                 "    7. py__init_obs_pdaf\n    " \
                                 "    8. py__init_obscovar_pdaf\n    " \
                                 "    9. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                 "    10. core DA algorith\n    " \
                                 "    11. py__prepoststep_state_pdaf\n    " \
                                 "    12. py__distribute_state_pdaf\n    " \
                                 "    13. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_lenkf`\n    " \
                                 "   and :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                 "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                 "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                 "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['assimilate_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                  "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                  "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                  "Local ESTKF (error space transform " \
                                  "Kalman filter) [1]_ for a single DA step without OMI.\n    " \
                                  "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                  "This function should be called at each model time step.\n    " \
                                  "The function is a combination of :func:`pyPDAF.PDAF.put_state_lestkf`\n    " \
                                  "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                  "User-supplied functions are executed in the following sequence:\n    " \
                                  "    1. py__collect_state_pdaf\n    "\
                                  "    2. py__prepoststep_state_pdaf\n    "\
                                  "    3. py__init_n_domains_p_pdaf\n    "\
                                  "    4. py__init_dim_obs_pdaf\n    "\
                                  "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                  "    6. py__init_obs_pdaf\n    " \
                                  "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                  "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                  "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                  "    8. loop over each local domain:\n    " \
                                  "        1. py__init_dim_l_pdaf\n    "\
                                  "        2. py__init_dim_obs_l_pdaf\n    "\
                                  "        3. py__g2l_state_pdaf\n    "\
                                  "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                  "        5. py__init_obs_l_pdaf\n    " \
                                  "        6. py__g2l_obs_pdaf\n    "\
                                  "           (localise each ensemble member in observation space)\n    "\
                                  "        7. py__init_obsvar_l_pdaf\n    " \
                                  "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                  "        8. py__prodRinvA_l_pdaf\n    "\
                                  "        9. core DA algorithm\n    " \
                                  "        10. py__l2g_state_pdaf\n    "\
                                  "    9. py__prepoststep_state_pdaf\n    "\
                                  "    10. py__distribute_state_pdaf\n    "\
                                  "    11. py__next_observation_pdaf\n" \
                                  "\n    " \
                                  ".. deprecated:: 1.0.0\n\n    " \
                                  "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                  "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                  "\n\n    " \
                                  "References\n    " \
                                  "----------\n    " \
                                  ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                  "       A unification of ensemble square root Kalman filters. \n    " \
                                  "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['assimilate_letkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                 "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                 "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Local ensemble transform Kalman filter (LETKF) [1]_ " \
                                 "for a single DA step without OMI. Implementation is based on [2]_.\n    " \
                                 "Note that the LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                 "This function should be called at each model time step.\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_letkf`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                 "User-supplied functions are executed in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    "\
                                 "    2. py__prepoststep_state_pdaf\n    "\
                                 "    3. py__init_n_domains_p_pdaf\n    "\
                                 "    4. py__init_dim_obs_pdaf\n    "\
                                 "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                 "    6. py__init_obs_pdaf\n    " \
                                 "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                 "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                 "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                 "    8. loop over each local domain:\n    " \
                                 "        1. py__init_dim_l_pdaf\n    "\
                                 "        2. py__init_dim_obs_l_pdaf\n    "\
                                 "        3. py__g2l_state_pdaf\n    "\
                                 "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                 "        5. py__init_obs_l_pdaf\n    " \
                                 "        6. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    "\
                                 "        7. py__init_obsvar_l_pdaf\n    " \
                                 "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                 "        8. py__prodRinvA_l_pdaf\n    "\
                                 "        9. core DA algorithm\n    " \
                                 "        10. py__l2g_state_pdaf\n    "\
                                 "    9. py__prepoststep_state_pdaf\n    "\
                                 "    10. py__distribute_state_pdaf\n    "\
                                 "    11. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                 "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007).\n    "\
                                 "       Efficient data assimilation for spatiotemporal chaos:\n    "\
                                 "       A local ensemble transform Kalman filter. \n    "\
                                 "       Physica D: Nonlinear Phenomena, 230(1-2), 112-126.\n    " \
                                 ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                 "       A unification of ensemble square root Kalman filters. \n    " \
                                 "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['assimilate_lseik'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                 "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Local singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n    " \
                                 "This function should be called at each model time step.\n\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_lseik` " \
                                 "and :func:`pyPDAF.PDAF.get_state`\n\n    "\
                                 "This function  executes the user-supplied functions in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_n_domains_p_pdaf\n    " \
                                 "    4. py__init_dim_obs_pdaf\n    " \
                                 "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    6. py__init_obs_pdaf\n    "\
                                 "       (if global adaptive forgetting factor `type_forget=1`\n    " \
                                 "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                 "    7. py__init_obsvar_pdaf\n    "\
                                 "       (if global adaptive forgetting factor is used)\n    " \
                                 "    8. loop over each local domain:\n    " \
                                 "        1. py__init_dim_l_pdaf\n    " \
                                 "        2. py__init_dim_obs_l_pdaf\n    " \
                                 "        3. py__g2l_state_pdaf\n    " \
                                 "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                 "        5. py__init_obs_l_pdaf\n    "\
                                 "        6. py__g2l_obs_pdaf\n    "\
                                 "           (localise each ensemble member in observation space)\n    " \
                                 "        7. py__init_obsvar_l_pdaf\n    "\
                                 "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                 "        8. py__prodRinvA_l_pdaf\n    " \
                                 "        9. core DA algorithm\n    " \
                                 "        10. py__l2g_state_pdaf\n    " \
                                 "    9. py__prepoststep_state_pdaf\n    " \
                                 "    10. py__distribute_state_pdaf\n    " \
                                 "    11. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                 "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                 "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                 "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['assimilate_lnetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                 "or :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Local Nonlinear Ensemble Transform Filter (LNETF) [1]_ for a single DA step.\n    " \
                                 "The nonlinear filter computes the distribution up to\n    " \
                                 "the second moment similar to Kalman filters but it uses a nonlinear weighting similar to\n    " \
                                 "particle filters. This leads to an equal weights assumption for the prior ensemble at each step.\n    " \
                                 "This function should be called at each model time step.\n\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.put_state_lnetf`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                 "This function executes the user-supplied function in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_n_domains_p_pdaf\n    " \
                                 "    4. py__init_dim_obs_pdaf\n    " \
                                 "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    6. loop over each local domain:\n    " \
                                 "        1. py__init_dim_l_pdaf\n    " \
                                 "        2. py__init_dim_obs_l_pdaf\n    " \
                                 "        3. py__g2l_state_pdaf\n    " \
                                 "        4. py__init_obs_l_pdaf\n    "\
                                 "        5. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                 "        6. py__likelihood_l_pdaf\n    " \
                                 "        7. core DA algorithm\n    " \
                                 "        8. py__l2g_state_pdaf\n    " \
                                 "    7. py__prepoststep_state_pdaf\n    " \
                                 "    8. py__distribute_state_pdaf\n    " \
                                 "    9. py__next_observation_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                 "   and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                 "       A second-order exact ensemble square root filter\n    " \
                                 "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                 "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['assimilate_lknetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                  "or :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`.\n\n    "\
                                  "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                  "A hybridised LETKF and LNETF [1]_ for a single DA step.\n    " \
                                  "The LNETF computes the distribution up to\n    " \
                                  "the second moment similar to Kalman filters but using a nonlinear weighting similar to\n    " \
                                  "particle filters. This leads to an equal weights assumption for the prior ensemble.\n    " \
                                  "The hybridisation with LETKF is expected to lead to improved performance for\n    " \
                                  "quasi-Gaussian problems.\n    " \
                                  "The function should be called at each model step.\n\n    " \
                                  "The function is a combination of :func:`pyPDAF.PDAF.put_state_lknetf`\n    " \
                                  "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                  "This function executes the user-supplied function in the following sequence:\n    " \
                                  "    1. py__collect_state_pdaf\n    " \
                                  "    2. py__prepoststep_state_pdaf\n    " \
                                  "    3. py__init_n_domains_p_pdaf\n    " \
                                  "    4. py__init_dim_obs_pdaf\n    " \
                                  "    5. py__obs_op_pdaf\n    "\
                                  "       (for each ensemble member)\n    " \
                                  "    6. py__init_obs_pdaf\n    " \
                                  "       (if global adaptive forgetting factor `type_forget=1`\n    "\
                                  "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                  "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                  "    8. loop over each local domain:\n    " \
                                  "        1. py__init_dim_l_pdaf\n    " \
                                  "        2. py__init_dim_obs_l_pdaf\n    " \
                                  "        3. py__g2l_state_pdaf\n    " \
                                  "        4. py__g2l_obs_pdaf\n    "\
                                  "           (localise each ensemble member in observation space)\n    " \
                                  "        5. py__init_obs_l_pdaf\n    "\
                                  "        6. py__init_obsvar_l_pdaf\n    "\
                                  "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                  "        7. py__prodRinvA_pdaf\n    " \
                                  "        8. py__likelihood_l_pdaf\n    " \
                                  "        9. core DA algorithm\n    " \
                                  "        10. py__l2g_state_pdaf\n    " \
                                  "    9. py__obs_op_pdaf\n    " \
                                  "       (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                  "    10. py__likelihood_hyb_l_pda\n    " \
                                  "    11. py__init_obsvar_l_pdaf\n    " \
                                  "        (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                  "    12. py__prodRinvA_hyb_l_pdaf\n    " \
                                  "    13. py__prepoststep_state_pdaf\n    " \
                                  "    14. py__distribute_state_pdaf\n    " \
                                  "    15. py__next_observation_pdaf\n" \
                                  "\n    " \
                                  ".. deprecated:: 1.0.0\n\n    " \
                                  "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                  "   and :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`" \
                                  "\n\n    " \
                                  "References\n    " \
                                  "----------\n    " \
                                  ".. [1] Nerger, L.. (2022) \n    " \
                                  "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                  "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['assimilate_prepost'] = "It is used to preprocess and postprocess of the ensemble.\n\n    " \
                                   "No DA is performed in this function.\n    " \
                                   "Compared to :func:`pyPDAF.PDAF.prepost`, this function sets assimilation flag, \n    " \
                                   "which means that it is acted as an assimilation in PDAF.\n\n    " \
                                   "The function is a combination of :func:`pyPDAF.PDAF.put_state_prepost`\n    " \
                                   "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                   "This function executes the user-supplied functions in the following sequence: \n    " \
                                   "    1. py__collect_state_pdaf\n    " \
                                   "    2. py__prepoststep_state_pdaf (preprocess, step < 0)\n    " \
                                   "    3. py__prepoststep_state_pdaf (postprocess, step > 0)\n    " \
                                   "    4. py__distribute_state_pdaf\n    " \
                                   "    5. py__next_observation_pdaf"
docstrings['generate_obs'] = "Generation of synthetic observations based on given error statistics and observation operator.\n\n    " \
                             "When diagonal observation error covariance matrix is used,\n    " \
                             "it is recommended to use :func:`pyPDAF.PDAF.omi_generate_obs` functionalities\n    "\
                             "for fewer user-supplied functions and improved efficiency.\n\n    " \
                             "The generated synthetic observations are based on each member of model forecast.\n    " \
                             "Therefore, an ensemble of observations can be obtained. In a typical experiment,\n    "\
                             "one may only need one ensemble member.\n    " \
                             "The implementation strategy is similar to an assimilation step. This means that, \n    " \
                             "one can reuse many user-supplied functions for assimilation and observation generation.\n\n    " \
                             "The function is a combination of :func:`pyPDAF.PDAF.put_state_generate_obs`\n    " \
                             "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                             "This function executes the user-supplied function in the following sequence:\n    " \
                             "    1. py__collect_state_pdaf\n    " \
                             "    2. py__prepoststep_state_pdaf\n    " \
                             "    3. py__init_dim_obs_pdaf\n    " \
                             "    4. py__obs_op_pda\n    " \
                             "    5. py__init_obserr_f_pdaf\n    " \
                             "    6. py__get_obs_f_pdaf\n    " \
                             "    7. py__prepoststep_state_pdaf\n    " \
                             "    8. py__distribute_state_pdaf\n    " \
                             "    9. py__next_observation_pdaf"

docstrings['put_state_3dvar'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_3dvar`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "3DVar DA for a single DA step.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_3dvar`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "When 3DVar is used, the background error covariance matrix\n    "\
                                "has to be modelled for cotrol variable transformation.\n    " \
                                "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                "This function should be called at each model time step.\n\n    " \
                                "User-supplied functions are executed in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf\n    " \
                                "    5. py__init_obs_pdaf\n    " \
                                "    6. Iterative optimisation:\n    " \
                                "        1. py__cvt_pdaf\n    " \
                                "        2. py__obs_op_lin_pdaf\n    " \
                                "        3. py__prodRinvA_pdaf\n    " \
                                "        4. py__obs_op_adj_pdaf\n    " \
                                "        5. py__cvt_adj_pdaf\n    " \
                                "        6. core DA algorithm\n    " \
                                "    7. py__cvt_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_3dvar`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`"
docstrings['put_state_en3dvar_estkf'] =  "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    "\
                                         "or :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`.\n\n    "\
                                         "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                         "3DEnVar for a single DA step.\n\n    " \
                                         "Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_estkf`, this function has no :func:`get_state` call.\n    " \
                                         "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                         "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                         "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                         "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                         "function call to ensure the sequential DA.\n\n    " \
                                         "The background error covariance matrix is estimated by an ensemble.\n    " \
                                         "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                         "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                         "This function should be called at each model time step.\n\n    " \
                                         "The user-supplied functions are executed in the following sequence:\n    " \
                                         "    1. py__collect_state_pdaf\n    " \
                                         "    2. py__prepoststep_state_pdaf\n    " \
                                         "    3. py__init_dim_obs_pdaf\n    " \
                                         "    4. py__obs_op_pdaf\n    " \
                                         "    5. py__init_obs_pdaf\n    " \
                                         "    6. the iterative optimisation:\n    " \
                                         "        1. py__cvt_ens_pdaf\n    " \
                                         "        2. py__obs_op_lin_pdaf\n    " \
                                         "        3. py__prodRinvA_pdaf\n    " \
                                         "        4. py__obs_op_adj_pdaf\n    " \
                                         "        5. py__cvt_adj_ens_pdaf\n    " \
                                         "        6. core 3DEnVar algorithm\n    " \
                                         "    7. py__cvt_ens_pdaf\n    " \
                                         "    8. ESTKF:\n    " \
                                         "        1. py__init_dim_obs_pdaf\n    " \
                                         "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                         "        3. py__init_obs_pdaf\n    " \
                                         "        4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                         "        5. py__init_obsvar_pdaf\n    " \
                                         "           (only relevant for adaptive forgetting factor schemes)\n    " \
                                         "        6. py__prodRinvA_pdaf\n    " \
                                         "        7. core ESTKF algorithm\n" \
                                         "\n    " \
                                         ".. deprecated:: 1.0.0\n\n    " \
                                         "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    " \
                                         "   and :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`"
docstrings['put_state_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    "\
                                         "or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`.\n\n    "\
                                         "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                         "3DEnVar for a single DA step without post-processing, distributing analysis, and setting next observation step,\n     "\
                                         "where the ensemble anomaly is generated by LESTKF.\n\n    " \
                                         "Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                         "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                         "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                         "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                         "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                         "function call to ensure the sequential DA.\n\n    " \
                                         "The background error covariance matrix is estimated by ensemble.\n    " \
                                         "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                         "An LESTKF is used to generate ensemble perturbations.\n    " \
                                         "This function should be called at each model time step.\n\n    " \
                                         "The user-supplied function are executed in the following sequence:\n    " \
                                         "    1. py__collect_state_pdaf\n    " \
                                         "    2. py__prepoststep_state_pdaf\n    " \
                                         "    3. py__init_dim_obs_pdaf\n    " \
                                         "    4. py__obs_op_pdaf\n    " \
                                         "    5. py__init_obs_pdaf\n    " \
                                         "    6. Starting the iterative optimisation:\n    " \
                                         "        1. py__cvt_ens_pdaf\n    " \
                                         "        2. py__obs_op_lin_pdaf\n    " \
                                         "        3. py__prodRinvA_pdaf\n    " \
                                         "        4. py__obs_op_adj_pdaf\n    " \
                                         "        5. py__cvt_adj_ens_pdaf\n    " \
                                         "        6. core DA algorithm\n    " \
                                         "    7. py__cvt_ens_pdaf\n    " \
                                         "    8. Perform LESTKF:\n    " \
                                         "        1. py__init_n_domains_p_pdaf\n    " \
                                         "        2. py__init_dim_obs_pdaf\n    " \
                                         "        3. py__obs_op_pdaf\n    "\
                                         "           (for each ensemble member)\n    " \
                                         "        4. py__init_obs_pdaf\n    "\
                                         "           (if global adaptive forgetting factor is used\n    "\
                                         "           `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                         "        5. py__init_obsvar_pdaf\n    "\
                                         "           (if global adaptive forgetting factor is used)\n    " \
                                         "        6. loop over each local domain:\n    " \
                                         "            1. py__init_dim_l_pdaf\n    " \
                                         "            2. py__init_dim_obs_l_pdaf\n    " \
                                         "            3. py__g2l_state_pdaf\n    " \
                                         "            4. py__g2l_obs_pdaf\n    "\
                                         "               (localise mean ensemble in observation space)\n    " \
                                         "            5. py__init_obs_l_pdaf\n    "\
                                         "            6. py__g2l_obs_pdaf\n    " \
                                         "               (localise each ensemble member in observation space)\n    " \
                                         "            7. py__init_obsvar_l_pdaf\n    " \
                                         "               (only called if local adaptive forgetting factor\n    " \
                                         "               `type_forget=2` is used)\n    "\
                                         "            8. py__prodRinvA_l_pdaf\n    " \
                                         "            9. core DA algorithm\n    " \
                                         "            10. py__l2g_state_pdaf\n" \
                                         "\n    " \
                                         ".. deprecated:: 1.0.0\n\n    " \
                                         "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
                                         "   and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`"
docstrings['put_state_hyb3dvar_estkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    "\
                                         "or :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`.\n\n    "\
                                         "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                         "Hybrid 3DEnVar for a single DA step where\n    " \
                                         "the background error covariance is hybridised by a static background error covariance,\n    " \
                                         "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                         "Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_estkf`, this function has no :func:`get_state` call.\n    " \
                                         "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                         "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                         "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                         "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                         "function call to ensure the sequential DA.\n\n    " \
                                         "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                         "ESTKF in this implementation.\n    " \
                                         "This function should be called at each model time step.\n\n    " \
                                         "The user-supplied functions are executed in the following sequence:\n    " \
                                         "    1. py__collect_state_pdaf\n    " \
                                         "    2. py__prepoststep_state_pdaf\n    " \
                                         "    3. py__init_dim_obs_pdaf\n    " \
                                         "    4. py__obs_op_pdaf\n    " \
                                         "    5. py__init_obs_pdaf\n    " \
                                         "    6. the iterative optimisation:\n    " \
                                         "        1. py__cvt_pdaf\n    " \
                                         "        2. py__cvt_ens_pdaf\n    " \
                                         "        3. py__obs_op_lin_pdaf\n    " \
                                         "        4. py__prodRinvA_pdaf\n    " \
                                         "        5. py__obs_op_adj_pdaf\n    " \
                                         "        6. py__cvt_adj_pdaf\n    " \
                                         "        7. py__cvt_adj_ens_pdaf\n    " \
                                         "        8. core 3DEnVar algorithm\n    " \
                                         "    7. py__cvt_pdaf\n    " \
                                         "    8. py__cvt_ens_pdaf\n    " \
                                         "    9. Perform ESTKF:\n    " \
                                         "        1. py__init_dim_obs_pdaf\n    " \
                                         "        2. py__obs_op_pdaf\n    "\
                                         "           (for ensemble mean)\n    " \
                                         "        3. py__init_obs_pdaf\n    " \
                                         "        4. py__obs_op_pdaf\n    " \
                                         "           (for each ensemble member)\n    " \
                                         "        5. py__init_obsvar_pdaf\n    " \
                                         "           (only relevant for adaptive forgetting factor schemes)\n    " \
                                         "        6. py__prodRinvA_pdaf\n    " \
                                         "        7. core ESTKF algorithm\n" \
                                         "\n    " \
                                         ".. deprecated:: 1.0.0\n\n    " \
                                         "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    " \
                                         "   and :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`"
docstrings['put_state_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    "\
                                          "or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                          "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                          "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                          "without post-processing, distributing analysis, and setting next observation step, where\n    " \
                                          "the background error covariance is hybridised by a static background error covariance,\n    " \
                                          "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                          "Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                          "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                          "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                          "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                          "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                          "function call to ensure the sequential DA.\n\n    " \
                                          "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                          "LESTKF in this implementation.\n    " \
                                          "This function should be called at each model time step.\n\n    " \
                                          "The user-supplied functions are executed in the following sequence:\n    " \
                                          "    1. py__collect_state_pdaf\n    " \
                                          "    2. py__prepoststep_state_pdaf\n    " \
                                          "    3. py__init_dim_obs_pdaf\n    " \
                                          "    4. py__obs_op_pdaf\n    " \
                                          "    5. py__init_obs_pdaf\n    " \
                                          "    6. The iterative optimisation:\n    " \
                                          "        1. py__cvt_pdaf\n    " \
                                          "        2. py__cvt_ens_pdaf\n    " \
                                          "        3. py__obs_op_lin_pdaf\n    " \
                                          "        4. py__prodRinvA_pdaf\n    " \
                                          "        5. py__obs_op_adj_pdaf\n    " \
                                          "        6. py__cvt_adj_pdaf\n    " \
                                          "        7. py__cvt_adj_ens_pdaf\n    " \
                                          "        8. core DA algorithm\n    " \
                                          "    7. py__cvt_pdaf\n    " \
                                          "    8. py__cvt_ens_pdaf\n    " \
                                          "    9. Perform LESTKF:\n    " \
                                          "        1. py__init_n_domains_p_pdaf\n    " \
                                          "        2. py__init_dim_obs_pdaf\n    " \
                                          "        3. py__obs_op_pdaf\n    " \
                                          "           (for each ensemble member)\n    " \
                                          "        4. py__init_obs_pdaf\n    " \
                                          "           (if global adaptive forgetting factor `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                          "        5. py__init_obsvar_pdaf\n    " \
                                          "           (if global adaptive forgetting factor is used)\n    " \
                                          "        6. loop over each local domain:\n    " \
                                          "            1. py__init_dim_l_pdaf\n    " \
                                          "            2. py__init_dim_obs_l_pdaf\n    " \
                                          "            3. py__g2l_state_pdaf\n    " \
                                          "            4. py__g2l_obs_pdaf\n    "\
                                          "               (localise mean ensemble in observation space)\n    " \
                                          "            5. py__init_obs_l_pdaf\n    "\
                                          "            6. py__g2l_obs_pdaf\n    " \
                                          "               (localise each ensemble member in observation space)\n    " \
                                          "            7. py__init_obsvar_l_pdaf\n    " \
                                          "               (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                          "            8. py__prodRinvA_l_pdaf\n    " \
                                          "            9. core DA algorithm\n    " \
                                          "            10. py__l2g_state_pdaf\n" \
                                          "\n    " \
                                          ".. deprecated:: 1.0.0\n\n    " \
                                          "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                          "   and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`"
docstrings['put_state_enkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                                 "or :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Stochastic EnKF (ensemble " \
                                 "Kalman filter) [1]_ for a single DA step without OMI.\n\n    " \
                                 "Compared to :func:`pyPDAF.PDAF.assimilate_enkf`, this function has no :func:`get_state` call.\n    " \
                                 "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                 "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                 "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                 "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                 "function call to ensure the sequential DA.\n\n    " \
                                 "This function should be called at each model time step. \n\n    " \
                                 "This function executes the user-supplied functions in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                 "    5. py__add_obs_err_pdaf\n    " \
                                 "    6. py__init_obs_pdaf\n    " \
                                 "    7. py__init_obscovar_pdaf\n    " \
                                 "    8. py__obs_op_pdaf (for each ensemble member)\n    " \
                                 "    9. core DA algorithm\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                                 "   and :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Evensen, G. (1994), \n    "\
                                 "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                 "       using Monte Carlo methods to forecast error statistics,\n    "\
                                 "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572."
docstrings['put_state_estkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                                "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR` instead of this function.\n\n    " \
                                "OMI functions need fewer user-supplied functions and improve DA efficiency.\n\n    " \
                                "This function calls ESTKF (error space transform Kalman filter) [1]_.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_estkf`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "The ESTKF is a more efficient equivalent to the ETKF.\n\n    " \
                                "The function should be called at each model time step.\n\n    " \
                                "User-supplied functions are executed in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                "    5. py__init_obs_pdaf\n    " \
                                "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                "    8. py__prodRinvA_pdaf\n    " \
                                "    9. core DA algorithm\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`." \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                "       A unification of ensemble square root Kalman filters. \n    " \
                                "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_etkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                               "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
                               "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                               "Using ETKF (ensemble transform " \
                               "Kalman filter) [1]_ for a single DA step without OMI. The implementation is baed on [2]_.\n\n    " \
                               "Compared to :func:`pyPDAF.PDAF.assimilate_etkf`, this function has no :func:`get_state` call.\n    " \
                               "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                               "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                               "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                               "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                               "function call to ensure the sequential DA.\n\n    " \
                               "This function should be called at each model time step.\n\n    " \
                               "This function executes the user-supplied function in the following sequence:\n    " \
                               "    1. py__collect_state_pdaf\n    " \
                               "    2. py__prepoststep_state_pdaf\n    " \
                               "    3. py__init_dim_obs_pdaf\n    " \
                               "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                               "    5. py__init_obs_pdaf\n    " \
                               "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                               "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                               "    8. py__prodRinvA_pdaf\n    " \
                               "    9. core DA algorithm\n" \
                               "\n    " \
                               ".. deprecated:: 1.0.0\n\n    " \
                               "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                               "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`" \
                               "\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Bishop, C. H., B. J. Etherton, and S. J. Majumdar (2001)\n    "\
                               "       Adaptive Sampling with the Ensemble Transform Kalman Filter.\n    "\
                               "       Part I: Theoretical Aspects. Mon. Wea. Rev., 129, 420–436,\n    "\
                               "       doi: 10.1175/1520-0493(2001)129<0420:ASWTET>2.0.CO;2. \n    " \
                               ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                               "       A unification of ensemble square root Kalman filters. \n    " \
                               "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_seek'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                               "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
                               "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                               "This function will use singular evolutive extended Kalman filter [1]_ for a single DA step.\n\n    " \
                               "Compared to :func:`pyPDAF.PDAF.assimilate_seek`, this function has no :func:`get_state` call.\n    " \
                               "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                               "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                               "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                               "function call to ensure the sequential DA.\n\n    " \
                               "This is a deterministic Kalman filter.\n    " \
                               "The function should be called at each model step.\n\n    " \
                               "This function executes the user-supplied functions in the following sequence:\n    " \
                               "    1. py__collect_state_pdaf\n    " \
                               "    2. py__prepoststep_state_pdaf\n    " \
                               "    3. py__init_dim_obs_pdaf\n    " \
                               "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                               "    5. py__init_obs_pdaf\n    " \
                               "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                               "    7. py__prodRinvA_pdaf\n    " \
                               "    8. core DA algorithm\n" \
                               "\n    " \
                               ".. deprecated:: 1.0.0\n\n    " \
                               "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                               "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`" \
                               "\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                               "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                               "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_seik'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                               "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
                               "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                               "This function will use singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n\n    " \
                               "Compared to :func:`pyPDAF.PDAF.assimilate_seik`, this function has no :func:`get_state` call.\n    " \
                               "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                               "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                               "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                               "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                               "function call to ensure the sequential DA.\n\n    " \
                               "The function should be called at each model step.\n\n    " \
                               "The function executes the user-supplied functions in the following sequence:\n    " \
                               "    1. py__collect_state_pdaf\n    " \
                               "    2. py__prepoststep_state_pdaf\n    " \
                               "    3. py__init_dim_obs_pdaf\n    " \
                               "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                               "    5. py__init_obs_pdaf\n    " \
                               "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                               "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                               "    8. py__prodRinvA_pdaf\n    " \
                               "    9. core DA algorithm\n" \
                               "\n    " \
                               ".. deprecated:: 1.0.0\n\n    " \
                               "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                               "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`" \
                               "\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                               "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                               "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_netf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                               "or :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`.\n\n    "\
                               "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                               "This function will use Nonlinear Ensemble Transform Filter (NETF) [1]_ \n    " \
                               "for a single DA step.\n\n    "\
                               "Compared to :func:`pyPDAF.PDAF.assimilate_netf`, this function has no :func:`get_state` call.\n    " \
                               "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                               "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                               "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                               "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                               "function call to ensure the sequential DA.\n\n    " \
                               "The nonlinear filter computes the distribution up to\n    " \
                               "the second moment similar to KF but using a nonlinear weighting similar to\n    " \
                               "particle filter. This leads to an equal weights assumption for prior ensemble.\n    " \
                               "The function should be called at each model step.\n\n    " \
                               "This function executes the user-supplied function in the following sequence:\n    " \
                               "    1. py__collect_state_pdaf\n    " \
                               "    2. py__prepoststep_state_pdaf\n    " \
                               "    3. py__init_dim_obs_pdaf\n    " \
                               "    4. py__init_obs_pdaf\n    " \
                               "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                               "    6. py__likelihood_pdaf\n    " \
                               "    7. core DA algorithm\n" \
                               "\n    " \
                               ".. deprecated:: 1.0.0\n\n    " \
                               "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                               "   and :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`" \
                               "\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                               "       A second-order exact ensemble square root filter\n    " \
                               "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                               "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['put_state_pf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
                              "or :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`.\n\n    "\
                              "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                              "This function will use particle filter for a single DA step.\n\n    " \
                              "Compared to :func:`pyPDAF.PDAF.assimilate_pf`, this function has no :func:`get_state` call.\n    " \
                              "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                              "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                              "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                              "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                              "function call to ensure the sequential DA.\n\n    " \
                              "This is a fully nonlinear filter, and may require a high number of ensemble members.\n    " \
                              "A review of particle filter can be found at [1]_.\n    " \
                              "The function should be called at each model step.\n\n    " \
                              "This function executes the user-supplied functions in the following sequence:\n    " \
                              "    1. py__collect_state_pdaf\n    " \
                              "    2. py__prepoststep_state_pdaf\n    " \
                              "    3. py__init_dim_obs_pdaf\n    " \
                              "    4. py__init_obs_pdaf\n    " \
                              "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                              "    6. py__likelihood_pdaf\n    " \
                              "    7. core DA algorithm\n" \
                              "\n    " \
                              ".. deprecated:: 1.0.0\n\n    " \
                              "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                              "   and :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`" \
                              "\n\n    " \
                              "References\n    " \
                              "----------\n    " \
                              ".. [1] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                              "       Particle filters for high‐dimensional geoscience applications:\n    "\
                              "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['put_state_lenkf'] = "It is recommended to use :func:`pyPDAF.PDAF.omi_put_state_lenkf`\n    "\
                                "or :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "Stochastic EnKF (ensemble Kalman filter) with covariance localisation [1]_\n    " \
                                "for a single DA step without OMI.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_lenkf`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                "This function should be called at each model time step.\n\n    " \
                                "The user-supplied function is executed in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_dim_obs_pdaf\n    " \
                                "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    5. py__localize_pdaf\n    " \
                                "    6. py__add_obs_err_pdaf\n    " \
                                "    7. py__init_obs_pdaf\n    " \
                                "    8. py__init_obscovar_pdaf\n    " \
                                "    9. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                "    10. core DA algorith\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_lenkf`\n    " \
                                "   and :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['put_state_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                 "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                 "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "Local ESTKF (error space transform " \
                                 "Kalman filter) [1]_ for a single DA step without OMI.\n\n    " \
                                 "Compared to :func:`pyPDAF.PDAF.assimilate_lestkf`, this function has no :func:`get_state` call.\n    " \
                                 "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                 "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                 "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                 "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                 "function call to ensure the sequential DA.\n\n    " \
                                 "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                 "This function should be called at each model time step.\n\n    " \
                                 "User-supplied functions are executed in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    "\
                                 "    2. py__prepoststep_state_pdaf\n    "\
                                 "    3. py__init_n_domains_p_pdaf\n    "\
                                 "    4. py__init_dim_obs_pdaf\n    "\
                                 "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                 "    6. py__init_obs_pdaf\n    " \
                                 "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                 "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                 "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                 "    8. loop over each local domain:\n    " \
                                 "        1. py__init_dim_l_pdaf\n    "\
                                 "        2. py__init_dim_obs_l_pdaf\n    "\
                                 "        3. py__g2l_state_pdaf\n    "\
                                 "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                 "        5. py__init_obs_l_pdaf\n    " \
                                 "        6. py__g2l_obs_pdaf\n    "\
                                 "           (localise each ensemble member in observation space)\n    "\
                                 "        7. py__init_obsvar_l_pdaf\n    " \
                                 "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                 "        8. py__prodRinvA_l_pdaf\n    "\
                                 "        9. core DA algorithm\n    " \
                                 "        10. py__l2g_state_pdaf\n"\
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                 "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                 "       A unification of ensemble square root Kalman filters. \n    " \
                                 "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_letkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "Local ensemble transform Kalman filter (LETKF) [1]_\n    " \
                                "for a single DA step without OMI. Implementation is based on [2]_.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_letkf`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "Note that the LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                "This function should be called at each model time step.\n\n    " \
                                "User-supplied functions are executed in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    "\
                                "    2. py__prepoststep_state_pdaf\n    "\
                                "    3. py__init_n_domains_p_pdaf\n    "\
                                "    4. py__init_dim_obs_pdaf\n    "\
                                "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                "    6. py__init_obs_pdaf\n    " \
                                "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                "    8. loop over each local domain:\n    " \
                                "        1. py__init_dim_l_pdaf\n    "\
                                "        2. py__init_dim_obs_l_pdaf\n    "\
                                "        3. py__g2l_state_pdaf\n    "\
                                "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                "        5. py__init_obs_l_pdaf\n    " \
                                "        6. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    "\
                                "        7. py__init_obsvar_l_pdaf\n    " \
                                "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                "        8. py__prodRinvA_l_pdaf\n    "\
                                "        9. core DA algorithm\n    " \
                                "        10. py__l2g_state_pdaf\n"\
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007).\n    "\
                                "       Efficient data assimilation for spatiotemporal chaos:\n    "\
                                "       A local ensemble transform Kalman filter. \n    "\
                                "       Physica D: Nonlinear Phenomena, 230(1-2), 112-126.\n    " \
                                ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                "       A unification of ensemble square root Kalman filters. \n    " \
                                "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_lseik'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "Local singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_lseik`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "This function should be called at each model time step.\n\n    " \
                                "This function  executes the user-supplied functions in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_n_domains_p_pdaf\n    " \
                                "    4. py__init_dim_obs_pdaf\n    " \
                                "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    6. py__init_obs_pdaf\n    "\
                                "       (if global adaptive forgetting factor `type_forget=1`\n    " \
                                "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                "    7. py__init_obsvar_pdaf\n    "\
                                "       (if global adaptive forgetting factor is used)\n    " \
                                "    8. loop over each local domain:\n    " \
                                "        1. py__init_dim_l_pdaf\n    " \
                                "        2. py__init_dim_obs_l_pdaf\n    " \
                                "        3. py__g2l_state_pdaf\n    " \
                                "        4. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                "        5. py__init_obs_l_pdaf\n    "\
                                "        6. py__g2l_obs_pdaf\n    "\
                                "           (localise each ensemble member in observation space)\n    " \
                                "        7. py__init_obsvar_l_pdaf\n    "\
                                "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                "        8. py__prodRinvA_l_pdaf\n    " \
                                "        9. core DA algorithm\n    " \
                                "        10. py__l2g_state_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_lnetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                "or :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`.\n\n    "\
                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                "Local Nonlinear Ensemble Transform Filter (LNETF) [1]_ for a single DA step.\n\n    " \
                                "Compared to :func:`pyPDAF.PDAF.assimilate_lnetf`, this function has no :func:`get_state` call.\n    " \
                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                "function call to ensure the sequential DA.\n\n    " \
                                "The nonlinear filter computes the distribution up to\n    " \
                                "the second moment similar to Kalman filters but it uses a nonlinear weighting similar to\n    " \
                                "particle filters. This leads to an equal weights assumption for the prior ensemble at each step.\n    " \
                                "This function should be called at each model time step.\n\n    " \
                                "This function executes the user-supplied function in the following sequence:\n    " \
                                "    1. py__collect_state_pdaf\n    " \
                                "    2. py__prepoststep_state_pdaf\n    " \
                                "    3. py__init_n_domains_p_pdaf\n    " \
                                "    4. py__init_dim_obs_pdaf\n    " \
                                "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                "    6. loop over each local domain:\n    " \
                                "        1. py__init_dim_l_pdaf\n    " \
                                "        2. py__init_dim_obs_l_pdaf\n    " \
                                "        3. py__g2l_state_pdaf\n    " \
                                "        4. py__init_obs_l_pdaf\n    "\
                                "        5. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                "        6. py__likelihood_l_pdaf\n    " \
                                "        7. core DA algorithm\n    " \
                                "        8. py__l2g_state_pdaf\n" \
                                "\n    " \
                                ".. deprecated:: 1.0.0\n\n    " \
                                "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`" \
                                "\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                "       A second-order exact ensemble square root filter\n    " \
                                "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['put_state_lknetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                 "or :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`.\n\n    "\
                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                 "A hybridised LETKF and LNETF [1]_ for a single DA step.\n\n    " \
                                 "Compared to :func:`pyPDAF.PDAF.assimilate_lknetf`, this function has no :func:`get_state` call.\n    " \
                                 "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                 "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                 "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                 "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                 "function call to ensure the sequential DA.\n\n    " \
                                 "The LNETF computes the distribution up to\n    " \
                                 "the second moment similar to Kalman filters but using a nonlinear weighting similar to\n    " \
                                 "particle filters. This leads to an equal weights assumption for the prior ensemble.\n    " \
                                 "The hybridisation with LETKF is expected to lead to improved performance for\n    " \
                                 "quasi-Gaussian problems.\n    " \
                                 "The function should be called at each model step.\n\n    " \
                                 "This function executes the user-supplied function in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_n_domains_p_pdaf\n    " \
                                 "    4. py__init_dim_obs_pdaf\n    " \
                                 "    5. py__obs_op_pdaf\n    "\
                                 "       (for each ensemble member)\n    " \
                                 "    6. py__init_obs_pdaf\n    " \
                                 "       (if global adaptive forgetting factor `type_forget=1`\n    "\
                                 "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                 "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                 "    8. loop over each local domain:\n    " \
                                 "        1. py__init_dim_l_pdaf\n    " \
                                 "        2. py__init_dim_obs_l_pdaf\n    " \
                                 "        3. py__g2l_state_pdaf\n    " \
                                 "        4. py__g2l_obs_pdaf\n    "\
                                 "           (localise each ensemble member in observation space)\n    " \
                                 "        5. py__init_obs_l_pdaf\n    "\
                                 "        6. py__init_obsvar_l_pdaf\n    "\
                                 "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                 "        7. py__prodRinvA_pdaf\n    " \
                                 "        8. py__likelihood_l_pdaf\n    " \
                                 "        9. core DA algorithm\n    " \
                                 "        10. py__l2g_state_pdaf\n    " \
                                 "    9. py__obs_op_pdaf\n    " \
                                 "       (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                 "    10. py__likelihood_hyb_l_pda\n    " \
                                 "    11. py__init_obsvar_l_pdaf\n    " \
                                 "        (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                 "    12. py__prodRinvA_hyb_l_pdaf\n" \
                                 "\n    " \
                                 ".. deprecated:: 1.0.0\n\n    " \
                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                 "   and :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`" \
                                 "\n\n    " \
                                 "References\n    " \
                                 "----------\n    " \
                                 ".. [1] Nerger, L.. (2022) \n" \
                                 "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n" \
                                 "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['put_state_prepost'] = "It is used to preprocess and postprocess of the ensemble.\n\n    " \
                                  "No DA is performed in this function.\n    " \
                                  "Compared to :func:`pyPDAF.PDAF.assimilate_prepost`, this function does not set assimilation flag, \n    " \
                                  "and does not distribute the processed ensemble to the model field.\n    " \
                                  "This function also does not set the next assimilation step as :func:`pyPDAF.PDAF.assimilate_prepost`\n    " \
                                  "because it does not call :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                  "This function executes the user-supplied functions in the following sequence: \n    " \
                                  "    1. py__collect_state_pdaf\n    " \
                                  "    2. py__prepoststep_state_pdaf (preprocess, step < 0)"
docstrings['put_state_generate_obs'] = "Generation of synthetic observations based on given error statistics and observation operator\n    " \
                                       "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                       "When diagonal observation error covariance matrix is used,\n    " \
                                       "it is recommended to use :func:`pyPDAF.PDAF.omi_generate_obs` functionalities\n    "\
                                       "for fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "The generated synthetic observations are based on each member of model forecast.\n    " \
                                       "Therefore, an ensemble of observations can be obtained. In a typical experiment,\n    "\
                                       "one may only need one ensemble member.\n\n    " \
                                       "Compared to :func:`pyPDAF.PDAF.generate_obs`, this function has no :func:`get_state` call.\n    " \
                                       "This means that the next DA step will not be assigned by user-supplied functions.\n    " \
                                       "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                       "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                       "function call to ensure the sequential DA.\n\n    " \
                                       "The implementation strategy is similar to an assimilation step. This means that, \n    " \
                                       "one can reuse many user-supplied functions for assimilation and observation generation.\n\n    " \
                                       "This function executes the user-supplied function in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    " \
                                       "    2. py__prepoststep_state_pdaf\n    " \
                                       "    3. py__init_dim_obs_pdaf\n    " \
                                       "    4. py__obs_op_pda\n    " \
                                       "    5. py__init_obserr_f_pdaf\n    " \
                                       "    6. py__get_obs_f_pdaf"

docstrings['deallocate'] = "Finalise the PDAF systems including freeing all memory used by PDAF."

docstrings['diag_effsample'] = "Calculating the effective sample size of a particle filter.\n\n    " \
                               "Based on [1]_, it is defined as the inverse of the sum of the squared particle filter weights:\n    " \
                               r":math:`N_{eff} = \frac{1}{\sum_{i=1}^{N} w_i^2}` where :math:`w_i` is the weight of particle with index i.""\n    " \
                               r"and :math:`N` is the number of particles.""\n\n    " \
                               r"If the :math:`N_{eff}=N`, all weights are identical, and the filter has no influence on the analysis.""\n    "  \
                               r"If :math:`N_{eff}=0`, the filter is collapsed.""\n\n    " \
                               "This is typically called during the analysis step of a particle filter,\n    "\
                               "e.g. in the analysis step of NETF and LNETF.\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Doucet, A., de Freitas, N., Gordon, N. (2001). \n    "\
                               "       An Introduction to Sequential Monte Carlo Methods. \n    "\
                               "       In: Doucet, A., de Freitas, N., Gordon, N. (eds) \n    "\
                               "       Sequential Monte Carlo Methods in Practice.\n    "\
                               "       Statistics for Engineering and Information Science.\n    "\
                               "       Springer, New York, NY. https://doi.org/10.1007/978-1-4757-3437-9_1"
docstrings['diag_ensstats']  = "Computing the skewness and kurtosis of the ensemble of a given element of the state vector.\n\n    " \
                               "The definition used for kurtosis follows that used by [1]_.\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Lawson, W. G., & Hansen, J. A. (2004).\n    "\
                               "       Implications of stochastic and deterministic filters as ensemble-based\n    "\
                               "       data assimilation methods in varying regimes of error growth.\n    "\
                               "       Monthly weather review, 132(8), 1966-1981."
docstrings['diag_histogram'] = "Computing the rank histogram of an ensemble.\n\n    " \
                               "A rank histogram is used to diagnose the reliability of the ensemble [1]_.\n    " \
                               "A perfectly reliable ensemble should have a uniform rank histogram.\n\n    " \
                               "The function can be called in the pre/poststep routine of PDAF\n    "\
                               "both before and after the analysis step to collect the histogram information.\n\n    " \
                               "References\n    " \
                               "----------\n    " \
                               ".. [1] Hamill, T. M. (2001).\n    " \
                               "       Interpretation of rank histograms for verifying ensemble forecasts.\n    " \
                               "       Monthly Weather Review, 129(3), 550-560."
docstrings['diag_CRPS_nompi'] = "A continuous rank probability score for an ensemble without using MPI parallelisation.\n\n    " \
                                "The implementation is based on [1]_.\n\n    " \
                                "References\n    " \
                                "----------\n    " \
                                ".. [1] Hersbach, H. (2000), \n    "\
                                "       Decomposition of the Continuous Ranked Probability Score for\n    " \
                                "       Ensemble Prediction Systems,\n    " \
                                "       Wea. Forecasting, 15, 559–570, doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2"
docstrings['diag_CRPS'] = "Obtain a continuous rank probability score for an ensemble.\n\n    " \
                          "The implementation is based on [1]_.\n\n    " \
                          "References\n    " \
                          "----------\n    " \
                          ".. [1] Hersbach, H. (2000), \n    "\
                          "       Decomposition of the Continuous Ranked Probability Score for\n    " \
                          "       Ensemble Prediction Systems,\n    " \
                          "       Wea. Forecasting, 15, 559–570, doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2"

docstrings['eofcovar'] = "EOF analysis of an ensemble of state vectors " \
                         "by singular value decomposition.\n\n    " \
                         "Typically, this function is used with :func:`pyPDAF.PDAF.sampleens`\n    " \
                         "to generate an ensemble of a chosen size \n    " \
                         "(up to the number of EOFs plus one).\n\n    " \
                         "Here, the function performs a singular value decomposition of the ensemble anomaly of the input matrix,\n    " \
                         "which is usually an ensemble formed by state vectors at multiple time steps.\n    " \
                         "The singular values and corresponding singular vectors can be used to\n    " \
                         "construct an error covariance matrix.\n    " \
                         "This can be used as the initial error covariance for the initial ensemble.\n\n    " \
                         "A multivariate scaling can be performed to ensure that " \
                         "all fields in the state vectors have unit variance.\n\n    " \
                         "It can be useful to store more EOFs than one finally\n    " \
                         "might want to use to have the flexibility to cary the ensemble size.\n\n    " \
                         "See also `PDAF webpage <https://pdaf.awi.de/trac/wiki/EnsembleGeneration>`_"

docstrings['gather_dim_obs_f'] = "Gathers the dimension of observation vector across multiple local domains/filter processors.\n\n    " \
                                 "This function is typically used in the user-supplied function of :func:`py__init_dim_obs_f_pdaf`\n    " \
                                 "when OMI functionality is not used. Otherwise, one can use :func:`pyPDAF.PDAF.omi_gather_obs`.\n\n    " \
                                 "This function does two things:\n    " \
                                 "    1. Receiving observation dimension on each local process and allocate arrays for observation vectors.\n    "\
                                 "    2. Gather the total dimension of observation across local process\n\n    " \
                                 "The second functionality is used for domain localised filters when local domains are distributed in different processors.\n    "\
                                 "This means :math:`npes_filter > 0`. This function must be used\n    " \
                                 "before :func:`pyPDAF.PDAF.gather_obs_f` or :func:`pyPDAF.PDAF.gather_obs_f2`."

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

docstrings['get_state'] = "Post-processing the analysis and distributing state vector back to the model.\n\n    " \
                          "This function also sets the next model step for assimilation, or end the entire assimilation.\n\n    " \
                          "The function executes the user-supplied function in the following sequence:\n\n    " \
                          "1. py__prepoststep_state_pdaf\n\n    " \
                          "2. py__distribute_state_pdaf\n\n    " \
                          "3. py__next_observation_pdaf"

docstrings['init'] = "This function initialises the PDAF system.\n\n    " \
                     "It is called once at the beginning of the assimilation.\n    " \
                     "The function specifies the type of DA methods, \n    " \
                     "parameters of the filters, the MPI communicators, and other parallel options.\n    " \
                     "The user-supplied function :func:`py__init_ens_pdaf`\n    " \
                     "provides an initial ensemble to the internal PDAF ensemble array.\n    " \
                     "The internal PDAF ensemble can be distribute to the model by\n    " \
                     ":func:`pyPDAF.PDAF.get_state`.\n\n    " \
                     "The filter options including `filtertype`, `subtype`, `param_int`, and `param_real`\n    " \
                     "are introduced in\n    " \
                     "`PDAF filter options wiki page <https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>`_.\n\n    " \
                     "The MPI communicators are defined in\n    " \
                     "`pyPDAF example page <https://github.com/yumengch/pyPDAF/blob/main/example/parallelisation.py>`_.\n    " \
                     "The example script is based on the parallelisation strategy of PDAF\n    " \
                     "which is available at\n    " \
                     "`PDAF parallelisation strategy wiki page <https://pdaf.awi.de/trac/wiki/ImplementationConceptOnline#Parallelizationofthedataassimilationprogram>`_\n    " \
                     "and `PDAF parallelisation adaptation wiki page <https://pdaf.awi.de/trac/wiki/AdaptParallelization>`_.\n    " \
                     "In most cases, the user does not need to change the parallelisation script.\n\n    " \

docstrings['local_weight'] = "The function is used for localisation in the analysis step of a filter " \
                             "and computes a weight according to the specified distance " \
                             "and the settings for the localising function. " \
                             "Typically the function is called in `py__prodRinvA_l_pdaf` " \
                             "in the domain-localised filters. " \
                             "Also, the function is typically called for the LEnKF " \
                             "in the `py__localize_covar_pdaf`. \n    " \
                             "This function is usually only used in user-codes that do not use PDAF-OMI."

docstrings['print_info'] = "Printing the wallclock time and memory measured by PDAF.\n\n    " \
                           "This is called at the end of the DA program.\n\n    " \
                           "The function displays the following information:\n    " \
                           "    - Memory required for the ensemble array, state vector, and transform matrix\n    " \
                           "    - Memory required by the analysis step\n    " \
                           "    - Memory required to perform the ensemble transformation"

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
docstrings['set_comm_pdaf'] = "Setting the MPI communicator used by PDAF.\n\n    " \
                              "Without using this function `MPI_COMM_WORLD` is used.\n    " \
                              "This function is very useful if a set of processors is dedicated for I/O or other operations."
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


docstrings['omi_init'] = "Allocating an array of `obs_f` derived types instances.\n\n    " \
                         "This function initialises the number of observation types.\n    " \
                         "This should be called at the start of the DA system after :func:`pyPDAF.PDAF.init`."
docstrings['omi_set_doassim'] = "Setting the `doassim` attribute of `obs_f`.\n\n    " \
                                "Properties of `obs_f` are typically set in user-supplied function\n    " \
                                "`py__init_dim_obs_pdaf`.\n\n    " \
                                "If `doassim` is set to 0,\n    " \
                                "the given type of observation is not assimilated in the DA system. "
docstrings['omi_set_disttype'] = "Setting the `doassim` attribute of `obs_f`.\n\n    " \
                                 "Properties of `obs_f` are typically set in user-supplied function\n    " \
                                 "`py__init_dim_obs_pdaf`.\n\n    " \
                                 "The `disttype` determines the way the distance\n    " \
                                 "between observation and model grid is calculated in OMI.\n    " \
                                 "See https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype."
docstrings['omi_set_ncoord'] = "This function sets the `ncoord` attribute of `obs_f` typically used in user-supplied function `py__init_dim_obs_pdaf`. " \
                               "This is the dimension of coordinates of the observation. "
docstrings['omi_set_id_obs_p'] = "Setting the `id_obs_p` attribute of `obs_f`.\n\n    " \
                                 "The function is typically used in user-supplied function `py__init_dim_obs_pdaf`.\n\n    " \
                                 "Here, `id_obs_p(nrows, dim_obs_p)` is a 2D array of integers.\n    " \
                                 "The value of `nrows` depends on the observation operator used for an observation.\n\n    " \
                                 "Examples:\n\n    " \
                                 "- `nrows=1`: observations are located on model grid point.\n    "\
                                 "  In this case, " \
                                 "`id_obs_p` stores the index of the state vector (starting from 1) corresponds to the observations,\n    " \
                                 "  e.g. `id_obs_p[0, j] = i` means that the location and variable of the `i`-th element of the state vector\n    " \
                                 "  is the same as the `j`-th observation.\n\n    " \
                                 "- `nrows=4`: each observation corresponds to 4 indices of elements in the state vector.\n    "\
                                 "   In this case,\n    " \
                                 "   the location of these elements is used to perform bi-linear interpolation\n    "\
                                 "   from model grid to observation location.\n    " \
                                 "   This information is used in the :func:`pyPDAF.PDAF.omi_obs_op_gridavg`\n    "\
                                 "   and :func:`pyPDAF.PDAF.omi_obs_op_interp_lin` functions.\n    " \
                                 "   When interpolation is needed,\n    "\
                                 "   the weighting of the interpolation is done\n    "\
                                 "   in the :func:`pyPDAF.PDAF.omi_get_interp_coeff_lin`,\n    " \
                                 "   :func:`pyPDAF.PDAF.omi_get_interp_coeff_lin1D`,\n    "\
                                 "   and :func:`pyPDAF.PDAF.omi_get_interp_coeff_tri` functions.\n    " \
                                 "   The details of interpolation setup can be found at\n    " \
                                 "   `PDAF wiki page " \
                                 "<https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_.\n"
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
docstrings['omi_obs_op_gridpoint'] = "A (partial) identity observation operator\n\n    " \
                                     "This observation operator is used when observations and model use the same grid. \n\n    " \
                                     "The observations operator selects state vectors where observations are present. \n\n    " \
                                     "The function is used in the user-supplied function `py__obs_op_pdaf`. \n\n    "
docstrings['omi_obs_op_gridavg'] = "Observation operator that average values on given model grid points.\n\n    " \
                                   "The averaged model grid points are specified in `id_obs_p` property of `obs_f`,\n    " \
                                   "which can be set in :func:`pyPDAF.PDAF.omi_set_id_obs_p`.\n\n    "  \
                                   "The function is used in the user-supplied function `py__obs_op_pdaf`. "
docstrings['omi_obs_op_interp_lin'] = "Observation operator that linearly interpolates model grid values to observation location.\n\n    " \
                                      "The grid points used by linear interpolation is specified in `id_obs_p` of `obs_f`,\n    " \
                                      "which can be set by :func:`pyPDAF.PDAF.omi_set_id_obs_p`.\n\n    " \
                                      "The function also requires `icoeff_p` attribute of `obs_f`,\n    " \
                                      "which can be set by :func:`pyPDAF.PDAF.omi_set_icoeff_p`\n\n    " \
                                      "The interpolation coefficient can be obtained by " \
                                      ":func:`pyPDAF.PDAF.omi_get_interp_coeff_lin1D`,\n    " \
                                      ":func:`pyPDAF.PDAF.omi_get_interp_coeff_lin`, and\n    " \
                                      ":func:`pyPDAF.PDAF.omi_get_interp_coeff_tri`\n\n    " \
                                      "The details of interpolation setup can be found at\n    `PDAF wiki page " \
                                      "<https://pdaf.awi.de/trac/wiki/OMI_observation_operators#Initializinginterpolationcoefficients>`_\n\n    " \
                                      "The function is used in the user-supplied function `py__obs_op_pdaf`. "
docstrings['omi_obs_op_adj_gridavg'] = "The adjoint observation operator of :func:`pyPDAF.PDAF.omi_obs_op_gridavg`."
docstrings['omi_obs_op_adj_gridpoint'] = "The adjoint observation operator of :func:`pyPDAF.PDAF.omi_obs_op_gridpoint`."
docstrings['omi_obs_op_adj_interp_lin'] = "The adjoint observation operator of :func:`pyPDAF.PDAF.omi_obs_op_interp_lin`."
docstrings['omi_get_interp_coeff_tri'] = "The coefficient for linear interpolation in 2D on unstructure triangular grid.\n\n    " \
                                         "The resulting coefficient is used in :func:`omi_obs_op_interp_lin`.\n\n    " \
                                         "This function is for triangular model grid interpolation coefficients " \
                                         "determined as barycentric coordinates."
docstrings['omi_get_interp_coeff_lin1D'] = "The coefficient for linear interpolation in 1D.\n\n    " \
                                           "The resulting coefficient is used in :func:`omi_obs_op_interp_lin`.\n\n    "
docstrings['omi_get_interp_coeff_lin'] = "The coefficient for linear interpolation up to 3D.\n\n    " \
                                         "The resulting coefficient is used in :func:`omi_obs_op_interp_lin`.\n\n    " \
                                         "See introduction in `PDAF-OMI wiki page \n    " \
                                         "<https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin>`_"

docstrings['omi_assimilate_3dvar'] = "3DVar DA for a single DA step using diagnoal observation error covariance matrix.\n\n    "\
                                     "See :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR` for non-diagonal observation error covariance matrix.\n\n    "\
                                     "When 3DVar is used, the background error covariance matrix\n    "\
                                     "has to be modelled for cotrol variable transformation.\n    " \
                                     "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                     "This function should be called at each model time step.\n\n    " \
                                     "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_3dvar`\n    " \
                                     "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                     "User-supplied functions are executed in the following sequence:\n    " \
                                     "    1. py__collect_state_pdaf\n    " \
                                     "    2. py__prepoststep_state_pdaf\n    " \
                                     "    3. py__init_dim_obs_pdaf\n    " \
                                     "    4. py__obs_op_pdaf\n    " \
                                     "    5. Iterative optimisation:\n    " \
                                     "        1. py__cvt_pdaf\n    " \
                                     "        2. py__obs_op_lin_pdaf\n    " \
                                     "        3. py__obs_op_adj_pdaf\n    " \
                                     "        4. py__cvt_adj_pdaf\n    " \
                                     "        5. core DA algorithm\n    " \
                                     "    6. py__cvt_pdaf\n    " \
                                     "    7. py__prepoststep_state_pdaf\n    " \
                                     "    8. py__distribute_state_pdaf\n    " \
                                     "    9. py__next_observation_pdaf"
docstrings['omi_assimilate_en3dvar_estkf'] = "3DEnVar for a single DA step using diagnoal observation error covariance matrix.\n\n    "\
                                             "See :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR` for using non-diagonal observation error covariance matirx.\n\n    " \
                                             "Here, the background error covariance matrix is estimated by an ensemble.\n    " \
                                             "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                             "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    " \
                                             "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                             "The user-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf\n    " \
                                             "    5. the iterative optimisation:\n    " \
                                             "        1. py__cvt_ens_pdaf\n    " \
                                             "        2. py__obs_op_lin_pdaf\n    " \
                                             "        3. py__obs_op_adj_pdaf\n    " \
                                             "        4. py__cvt_adj_ens_pdaf\n    " \
                                             "        5. core 3DEnVar algorithm\n    " \
                                             "    6. py__cvt_ens_pdaf\n    " \
                                             "    7. ESTKF:\n    " \
                                             "        1. py__init_dim_obs_pdaf\n    " \
                                             "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                             "        3. py__obs_op_pdaf (for each ensemble member)\n    " \
                                             "        4. core ESTKF algorithm\n    " \
                                             "    8. py__prepoststep_state_pdaf\n    " \
                                             "    9. py__distribute_state_pdaf\n    " \
                                             "    10. py__next_observation_pdaf"
docstrings['omi_assimilate_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    "\
                                              "or :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`.\n\n    "\
                                              "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                              "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using diagnoal observation error covariance matrix.\n    " \
                                              "The background error covariance matrix is estimated by ensemble.\n    " \
                                              "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                              "An LESTKF is used to generate ensemble perturbations.\n    " \
                                              "This function should be called at each model time step.\n\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_en3dvar_lestkf`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                              "The user-supplied function are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf\n    " \
                                              "    5. Starting the iterative optimisation:\n    " \
                                              "        1. py__cvt_ens_pdaf\n    " \
                                              "        2. py__obs_op_lin_pdaf\n    " \
                                              "        3. py__obs_op_adj_pdaf\n    " \
                                              "        4. py__cvt_adj_ens_pdaf\n    " \
                                              "        5. core DA algorithm\n    " \
                                              "    6. py__cvt_ens_pdaf\n    " \
                                              "    7. Perform LESTKF:\n    " \
                                              "        1. py__init_n_domains_p_pdaf\n    " \
                                              "        2. py__init_dim_obs_pdaf\n    " \
                                              "        3. py__obs_op_pdaf\n    "\
                                              "           (for each ensemble member)\n    " \
                                              "        4. loop over each local domain:\n    " \
                                              "            1. py__init_dim_l_pdaf\n    " \
                                              "            2. py__init_dim_obs_l_pdaf\n    " \
                                              "            3. py__g2l_state_pdaf\n    " \
                                              "            4. core DA algorithm\n    " \
                                              "            5. py__l2g_state_pdaf\n    " \
                                              "    8. py__prepoststep_state_pdaf\n    " \
                                              "    9. py__distribute_state_pdaf\n    " \
                                              "    10. py__next_observation_pdaf\n" \
                                              "\n    " \
                                              ".. deprecated:: 1.0.0\n\n    " \
                                              "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    " \
                                              "   and :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`"
docstrings['omi_assimilate_hyb3dvar_estkf'] = "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix.\n\n    "\
                                              "See :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR` for non-diagonal observation error covariance matrix.\n\n    " \
                                              "Here the background error covariance is hybridised by a static background error covariance,\n    " \
                                              "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                              "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                              "ESTKF in this implementation.\n    " \
                                              "This function should be called at each model time step.\n\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                              "The user-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf\n    " \
                                              "    5. the iterative optimisation:\n    " \
                                              "        1. py__cvt_pdaf\n    " \
                                              "        2. py__cvt_ens_pdaf\n    " \
                                              "        3. py__obs_op_lin_pdaf\n    " \
                                              "        4. py__obs_op_adj_pdaf\n    " \
                                              "        5. py__cvt_adj_pdaf\n    " \
                                              "        6. py__cvt_adj_ens_pdaf\n    " \
                                              "        7. core 3DEnVar algorithm\n    " \
                                              "    6. py__cvt_pdaf\n    " \
                                              "    7. py__cvt_ens_pdaf\n    " \
                                              "    8. Perform ESTKF:\n    " \
                                              "        1. py__init_dim_obs_pdaf\n    " \
                                              "        2. py__obs_op_pdaf\n    "\
                                              "           (for ensemble mean)\n    " \
                                              "        3. py__obs_op_pdaf\n    " \
                                              "           (for each ensemble member)\n    " \
                                              "        4. core ESTKF algorithm\n    " \
                                              "    9. py__prepoststep_state_pdaf\n    " \
                                              "    10. py__distribute_state_pdaf\n    " \
                                              "    11. py__next_observation_pdaf"
docstrings['omi_assimilate_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    "\
                                               "or :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                               "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                               "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix where\n    " \
                                               "the background error covariance is hybridised by a static background error covariance,\n    " \
                                               "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                               "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                               "LESTKF in this implementation.\n    " \
                                               "This function should be called at each model time step.\n\n    " \
                                               "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_lestkf`\n    " \
                                               "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                               "The user-supplied functions are executed in the following sequence:\n    " \
                                               "    1. py__collect_state_pdaf\n    " \
                                               "    2. py__prepoststep_state_pdaf\n    " \
                                               "    3. py__init_dim_obs_pdaf\n    " \
                                               "    4. py__obs_op_pdaf\n    " \
                                               "    5. The iterative optimisation:\n    " \
                                               "        1. py__cvt_pdaf\n    " \
                                               "        2. py__cvt_ens_pdaf\n    " \
                                               "        3. py__obs_op_lin_pdaf\n    " \
                                               "        4. py__obs_op_adj_pdaf\n    " \
                                               "        5. py__cvt_adj_pdaf\n    " \
                                               "        6. py__cvt_adj_ens_pdaf\n    " \
                                               "        7. core DA algorithm\n    " \
                                               "    6. py__cvt_pdaf\n    " \
                                               "    7. py__cvt_ens_pdaf\n    " \
                                               "    8. Perform LESTKF:\n    " \
                                               "        1. py__init_n_domains_p_pdaf\n    " \
                                               "        2. py__init_dim_obs_pdaf\n    " \
                                               "        3. py__obs_op_pdaf\n    " \
                                               "           (for each ensemble member)\n    " \
                                               "        4. loop over each local domain:\n    " \
                                               "            1. py__init_dim_l_pdaf\n    " \
                                               "            2. py__init_dim_obs_l_pdaf\n    " \
                                               "            3. py__g2l_state_pdaf\n    " \
                                               "            4. core DA algorithm\n    " \
                                               "            5. py__l2g_state_pdaf\n    " \
                                               "    9. py__prepoststep_state_pdaf\n    " \
                                               "    10. py__distribute_state_pdaf\n    " \
                                               "    11. py__next_observation_pdaf\n" \
                                               "\n    " \
                                               ".. deprecated:: 1.0.0\n\n    " \
                                               "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    " \
                                               "   and :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`"
docstrings['omi_assimilate_global'] =  "Global filters except for 3DVar for a single DA step using diagnoal observation error covariance matrix.\n\n    "\
                                       "See :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`, \n    " \
                                       "or :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`,\n    " \
                                       "or :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR` for non-diagonal observation error covariance matrix.\n\n    " \
                                       "Here, this function call is used for global stochastic EnKF [1]_, E(S)TKF [2]_, \n    " \
                                       "SEEK [2]_, SEIK [2]_, NETF [3]_, and particle filter [4]_.\n    " \
                                       "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                       "This function should be called at each model time step. \n\n    " \
                                       "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_global` " \
                                       "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                       "This function executes the user-supplied functions in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    " \
                                       "    2. py__prepoststep_state_pdaf\n    " \
                                       "    3. py__init_dim_obs_pdaf\n    " \
                                       "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                       "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                       "    6. core DA algorithm\n    " \
                                       "    7. py__prepoststep_state_pdaf\n    " \
                                       "    8. py__distribute_state_pdaf\n    " \
                                       "    9. py__next_observation_pdaf" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Evensen, G. (1994), \n    "\
                                       "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                       "       using Monte Carlo methods to forecast error statistics,\n    "\
                                       "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572.\n    " \
                                       ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                       "       A unification of ensemble square root Kalman filters. \n    " \
                                       "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                       ".. [3] Tödter, J., and B. Ahrens, 2015:\n    "\
                                       "       A second-order exact ensemble square root filter\n    " \
                                       "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                       "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                       ".. [4] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                                       "       Particle filters for high‐dimensional geoscience applications:\n    "\
                                       "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['omi_assimilate_lenkf'] = "Covariance localised stochastic EnKF for a single DA step using diagnoal observation error covariance matrix.\n\n    " \
                                     "See :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR` for non-diagnoal observation error covariance matrix.\n\n    "\
                                     "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                     "The implementation is based on [1]_.\n\n    " \
                                     "This function should be called at each model time step.\n    " \
                                     "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_lenkf`\n    " \
                                     "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                     "The user-supplied function is executed in the following sequence:\n    " \
                                     "    1. py__collect_state_pdaf\n    " \
                                     "    2. py__prepoststep_state_pdaf\n    " \
                                     "    3. py__init_dim_obs_pdaf\n    " \
                                     "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                     "    5. py__localize_pdaf\n    " \
                                     "    6. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                     "    7. core DA algorith\n    " \
                                     "    8. py__prepoststep_state_pdaf\n    " \
                                     "    9. py__distribute_state_pdaf\n    " \
                                     "    10. py__next_observation_pdaf" \
                                     "\n\n    " \
                                     "References\n    " \
                                     "----------\n    " \
                                     ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                     "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                     "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                     "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['omi_assimilate_local'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                     "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`,\n    " \
                                     "or :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`,\n    " \
                                     "or :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`.\n\n    " \
                                     "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                     "Domain local filters for a single DA step using diagnoal observation error covariance matrix.\n    " \
                                     "Here, this function call is used for LE(S)TKF [1]_, \n    " \
                                     "LSEIK [1]_, LNETF [2]_, and LKNETF [3]_.\n    " \
                                     "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                     "This function should be called at each model time step.\n    " \
                                     "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_local`\n    " \
                                     "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                     "User-supplied functions are executed in the following sequence:\n    " \
                                     "    1. py__collect_state_pdaf\n    "\
                                     "    2. py__prepoststep_state_pdaf\n    "\
                                     "    3. py__init_n_domains_p_pdaf\n    "\
                                     "    4. py__init_dim_obs_pdaf\n    "\
                                     "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                     "    6. loop over each local domain:\n    " \
                                     "        1. py__init_dim_l_pdaf\n    "\
                                     "        2. py__init_dim_obs_l_pdaf\n    "\
                                     "        3. py__g2l_state_pdaf\n    "\
                                     "        4. core DA algorithm\n    " \
                                     "        5. py__l2g_state_pdaf\n    "\
                                     "    7. py__prepoststep_state_pdaf\n    "\
                                     "    8. py__distribute_state_pdaf\n    "\
                                     "    9. py__next_observation_pdaf\n" \
                                     "\n    " \
                                     ".. deprecated:: 1.0.0\n\n    " \
                                     "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                     "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`,\n    " \
                                     "   and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`,\n    " \
                                     "   and :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`." \
                                     "\n\n    " \
                                     "References\n    " \
                                     "----------\n    " \
                                     ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                     "       A unification of ensemble square root Kalman filters. \n    " \
                                     "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                     ".. [2] Tödter, J., and B. Ahrens, 2015:\n    "\
                                     "       A second-order exact ensemble square root filter\n    " \
                                     "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                     "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                     ".. [3] Nerger, L.. (2022) \n    " \
                                     "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                     "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['omi_generate_obs'] = "Generation of synthetic observations based on given error statistics and observation operator for diagonal observation error covariance matrix.\n\n    " \
                                 "If non-diagonal observation error covariance matrix has to be used,\n    " \
                                 "the generic :func:`pyPDAF.PDAF.generate_obs` can be used.\n\n    "\
                                 "The generated synthetic observations are based on each member of model forecast.\n    " \
                                 "Therefore, an ensemble of observations can be obtained. In a typical experiment,\n    "\
                                 "one may only need one ensemble member.\n    " \
                                 "The implementation strategy is similar to an assimilation step. This means that, \n    " \
                                 "one can reuse many user-supplied functions for assimilation and observation generation.\n\n    " \
                                 "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_generate_obs`\n    " \
                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                 "This function executes the user-supplied function in the following sequence:\n    " \
                                 "    1. py__collect_state_pdaf\n    " \
                                 "    2. py__prepoststep_state_pdaf\n    " \
                                 "    3. py__init_dim_obs_pdaf\n    " \
                                 "    4. py__obs_op_pda\n    " \
                                 "    5. py__get_obs_f_pdaf\n    " \
                                 "    6. py__prepoststep_state_pdaf\n    " \
                                 "    7. py__distribute_state_pdaf\n    " \
                                 "    8. py__next_observation_pdaf"

docstrings['omi_put_state_3dvar'] = "3DVar DA for a single DA step using diagnoal observation error covariance matrix\n    " \
                                    "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                    "See :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR` for non-diagonal observation error covariance matrix.\n\n    "\
                                    "Compared to :func:`pyPDAF.PDAF.omi_assimilate_3dvar`, this function has no :func:`get_state` call.\n    " \
                                    "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                    "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                    "function call to ensure the sequential DA.\n\n    " \
                                    "When 3DVar is used, the background error covariance matrix\n    "\
                                    "has to be modelled for cotrol variable transformation.\n    " \
                                    "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                    "This function should be called at each model time step.\n\n    " \
                                    "User-supplied functions are executed in the following sequence:\n    " \
                                    "    1. py__collect_state_pdaf\n    " \
                                    "    2. py__prepoststep_state_pdaf\n    " \
                                    "    3. py__init_dim_obs_pdaf\n    " \
                                    "    4. py__obs_op_pdaf\n    " \
                                    "    5. Iterative optimisation:\n    " \
                                    "        1. py__cvt_pdaf\n    " \
                                    "        2. py__obs_op_lin_pdaf\n    " \
                                    "        3. py__obs_op_adj_pdaf\n    " \
                                    "        4. py__cvt_adj_pdaf\n    " \
                                    "        5. core DA algorithm\n    " \
                                    "    6. py__cvt_pdaf"
docstrings['omi_put_state_en3dvar_estkf'] = "3DEnVar for a single DA step using diagnoal observation error covariance matrix\n    " \
                                            "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                            "See :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR` for non-diagonal observation error covariance matrix.\n\n    " \
                                            "Compared to :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`, this function has no :func:`get_state` call.\n    " \
                                            "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                            "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                            "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                            "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                            "function call to ensure the sequential DA.\n\n    " \
                                            "The background error covariance matrix is estimated by an ensemble.\n    " \
                                            "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                            "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                            "This function should be called at each model time step.\n\n    " \
                                            "The user-supplied functions are executed in the following sequence:\n    " \
                                            "    1. py__collect_state_pdaf\n    " \
                                            "    2. py__prepoststep_state_pdaf\n    " \
                                            "    3. py__init_dim_obs_pdaf\n    " \
                                            "    4. py__obs_op_pdaf\n    " \
                                            "    5. the iterative optimisation:\n    " \
                                            "        1. py__cvt_ens_pdaf\n    " \
                                            "        2. py__obs_op_lin_pdaf\n    " \
                                            "        3. py__obs_op_adj_pdaf\n    " \
                                            "        4. py__cvt_adj_ens_pdaf\n    " \
                                            "        5. core 3DEnVar algorithm\n    " \
                                            "    6. py__cvt_ens_pdaf\n    " \
                                            "    7. ESTKF:\n    " \
                                            "        1. py__init_dim_obs_pdaf\n    " \
                                            "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                            "        3. py__obs_op_pdaf (for each ensemble member)\n    " \
                                            "        4. core ESTKF algorithm"
docstrings['omi_put_state_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    "\
                                             "or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`.\n\n    "\
                                             "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                             "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using diagnoal observation error covariance matrix.\n\n    " \
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                             "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "The background error covariance matrix is estimated by ensemble.\n    " \
                                             "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                             "An LESTKF is used to generate ensemble perturbations.\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "The user-supplied function are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf\n    " \
                                             "    5. Starting the iterative optimisation:\n    " \
                                             "        1. py__cvt_ens_pdaf\n    " \
                                             "        2. py__obs_op_lin_pdaf\n    " \
                                             "        3. py__obs_op_adj_pdaf\n    " \
                                             "        4. py__cvt_adj_ens_pdaf\n    " \
                                             "        5. core DA algorithm\n    " \
                                             "    6. py__cvt_ens_pdaf\n    " \
                                             "    7. Perform LESTKF:\n    " \
                                             "        1. py__init_n_domains_p_pdaf\n    " \
                                             "        2. py__init_dim_obs_pdaf\n    " \
                                             "        3. py__obs_op_pdaf\n    "\
                                             "           (for each ensemble member)\n    " \
                                             "        4. loop over each local domain:\n    " \
                                             "            1. py__init_dim_l_pdaf\n    " \
                                             "            2. py__init_dim_obs_l_pdaf\n    " \
                                             "            3. py__g2l_state_pdaf\n    " \
                                             "            4. core DA algorithm\n    " \
                                             "            5. py__l2g_state_pdaf\n" \
                                             "\n    " \
                                             ".. deprecated:: 1.0.0\n\n    " \
                                             "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
                                             "   and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`"
docstrings['omi_put_state_hyb3dvar_estkf'] = "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix\n    "\
                                             "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                             "See :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR` for non-diagonal observation error covariance matrix.\n\n    " \
                                             "Hybrid 3DEnVar for a single DA step where\n    " \
                                             "the background error covariance is hybridised by a static background error covariance,\n    " \
                                             "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                             "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                             "ESTKF in this implementation.\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "The user-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf\n    " \
                                             "    5. the iterative optimisation:\n    " \
                                             "        1. py__cvt_pdaf\n    " \
                                             "        2. py__cvt_ens_pdaf\n    " \
                                             "        3. py__obs_op_lin_pdaf\n    " \
                                             "        4. py__obs_op_adj_pdaf\n    " \
                                             "        5. py__cvt_adj_pdaf\n    " \
                                             "        6. py__cvt_adj_ens_pdaf\n    " \
                                             "        7. core 3DEnVar algorithm\n    " \
                                             "    6. py__cvt_pdaf\n    " \
                                             "    7. py__cvt_ens_pdaf\n    " \
                                             "    8. Perform ESTKF:\n    " \
                                             "        1. py__init_dim_obs_pdaf\n    " \
                                             "        2. py__obs_op_pdaf\n    "\
                                             "           (for ensemble mean)\n    " \
                                             "        3. py__obs_op_pdaf\n    " \
                                             "           (for each ensemble member)\n    " \
                                             "        4. core ESTKF algorithm\n"
docstrings['omi_put_state_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    "\
                                              "or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                              "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                              "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix\n    " \
                                              "without post-processing, distributing analysis, and setting next observation step, where\n    " \
                                              "the background error covariance is hybridised by a static background error covariance,\n    " \
                                              "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                              "Compared to :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                              "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                              "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                              "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                              "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                              "function call to ensure the sequential DA.\n\n    " \
                                              "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                              "LESTKF in this implementation.\n    " \
                                              "This function should be called at each model time step.\n\n    " \
                                              "The user-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf\n    " \
                                              "    5. The iterative optimisation:\n    " \
                                              "        1. py__cvt_pdaf\n    " \
                                              "        2. py__cvt_ens_pdaf\n    " \
                                              "        3. py__obs_op_lin_pdaf\n    " \
                                              "        4. py__obs_op_adj_pdaf\n    " \
                                              "        5. py__cvt_adj_pdaf\n    " \
                                              "        6. py__cvt_adj_ens_pdaf\n    " \
                                              "        7. core DA algorithm\n    " \
                                              "    6. py__cvt_pdaf\n    " \
                                              "    7. py__cvt_ens_pdaf\n    " \
                                              "    8. Perform LESTKF:\n    " \
                                              "        1. py__init_n_domains_p_pdaf\n    " \
                                              "        2. py__init_dim_obs_pdaf\n    " \
                                              "        3. py__obs_op_pdaf\n    " \
                                              "           (for each ensemble member)\n    " \
                                              "        4. loop over each local domain:\n    " \
                                              "            1. py__init_dim_l_pdaf\n    " \
                                              "            2. py__init_dim_obs_l_pdaf\n    " \
                                              "            3. py__g2l_state_pdaf\n    " \
                                              "            4. core DA algorithm\n    " \
                                              "            5. py__l2g_state_pdaf\n" \
                                              "\n    " \
                                              ".. deprecated:: 1.0.0\n\n    " \
                                              "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                              "   and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`"
docstrings['omi_put_state_global'] = "Global filters except for 3DVar for a single DA step using diagnoal observation error covariance matrix\n    " \
                                     "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                     "See :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR`, \n    " \
                                     "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`,\n    " \
                                     "or :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR` for non-diagonal observation error covariance matrix.\n\n    " \
                                     "OMI functions need fewer user-supplied functions and improve DA efficiency.\n\n    " \
                                     "Here, this function call is used for global stochastic EnKF [1]_, E(S)TKF [2]_, \n    " \
                                     "SEEK [2]_, SEIK [2]_, NETF [3]_, and particle filter [4]_.\n    " \
                                     "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                     "Compared to :func:`pyPDAF.PDAF.omi_assimilate_global`, this function has no :func:`get_state` call.\n    " \
                                     "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                     "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                     "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                     "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                     "function call to ensure the sequential DA.\n\n    " \
                                     "The ESTKF is a more efficient equivalent to the ETKF.\n\n    " \
                                     "The function should be called at each model time step.\n\n    " \
                                     "User-supplied functions are executed in the following sequence:\n    " \
                                     "    1. py__collect_state_pdaf\n    " \
                                     "    2. py__prepoststep_state_pdaf\n    " \
                                     "    3. py__init_dim_obs_pdaf\n    " \
                                     "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                     "    5. py__init_obs_pdaf\n    " \
                                     "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
                                     "    7. py__init_obsvar_pdaf (only relevant for adaptive forgetting factor schemes)\n    " \
                                     "    8. py__prodRinvA_pdaf\n    " \
                                     "    9. core DA algorithm\n" \
                                     "\n    " \
                                     ".. deprecated:: 1.0.0\n\n    " \
                                     "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
                                     "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`." \
                                     "\n\n    " \
                                     "References\n    " \
                                     "----------\n    " \
                                     ".. [1] Evensen, G. (1994), \n    "\
                                     "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                     "       using Monte Carlo methods to forecast error statistics,\n    "\
                                     "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572.\n    " \
                                     ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                     "       A unification of ensemble square root Kalman filters. \n    " \
                                     "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                     ".. [3] Tödter, J., and B. Ahrens, 2015:\n    "\
                                     "       A second-order exact ensemble square root filter\n    " \
                                     "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                     "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                     ".. [4] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                                     "       Particle filters for high‐dimensional geoscience applications:\n    "\
                                     "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['omi_put_state_lenkf'] = "Stochastic EnKF (ensemble Kalman filter) with covariance localisation using diagnoal observation error covariance matrix\n    "\
                                    "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                    "See :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR` for non-diagonal observation error covariance matrix.\n\n    "\
                                    "Stochastic EnKF (ensemble Kalman filter) with covariance localisation [1]_\n    " \
                                    "for a single DA step.\n\n    " \
                                    "Compared to :func:`pyPDAF.PDAF.omi_assimilate_lenkf`, this function has no :func:`get_state` call.\n    " \
                                    "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                    "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                    "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                    "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                    "function call to ensure the sequential DA.\n\n    " \
                                    "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                    "This function should be called at each model time step.\n\n    " \
                                    "The user-supplied function is executed in the following sequence:\n    " \
                                    "    1. py__collect_state_pdaf\n    " \
                                    "    2. py__prepoststep_state_pdaf\n    " \
                                    "    3. py__init_dim_obs_pdaf\n    " \
                                    "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                    "    5. py__localize_pdaf\n    " \
                                    "    6. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                    "    7. core DA algorith" \
                                    "\n\n    " \
                                    "References\n    " \
                                    "----------\n    " \
                                    ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                    "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                    "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                    "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['omi_put_state_local'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                    "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`,\n    " \
                                    "or :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`,\n    " \
                                    "or :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`.\n\n    " \
                                    "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                    "Domain local filters for a single DA step using diagnoal observation error covariance matrix\n    " \
                                    "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                    "Here, this function call is used for LE(S)TKF [1]_, \n    " \
                                    "LSEIK [1]_, LNETF [2]_, and LKNETF [3]_.\n    " \
                                    "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                    "Compared to :func:`pyPDAF.PDAF.omi_assimilate_local`, this function has no :func:`get_state` call.\n    " \
                                    "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                    "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                    "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                    "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                    "function call to ensure the sequential DA.\n\n    " \
                                    "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                    "This function should be called at each model time step.\n\n    " \
                                    "User-supplied functions are executed in the following sequence:\n    " \
                                    "    1. py__collect_state_pdaf\n    "\
                                    "    2. py__prepoststep_state_pdaf\n    "\
                                    "    3. py__init_n_domains_p_pdaf\n    "\
                                    "    4. py__init_dim_obs_pdaf\n    "\
                                    "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                    "    6. loop over each local domain:\n    " \
                                    "        1. py__init_dim_l_pdaf\n    "\
                                    "        2. py__init_dim_obs_l_pdaf\n    "\
                                    "        3. py__g2l_state_pdaf\n    "\
                                    "        4. core DA algorithm\n    " \
                                    "        5. py__l2g_state_pdaf\n"\
                                    "\n    " \
                                    ".. deprecated:: 1.0.0\n\n    " \
                                    "   This function is replaced by :func:`pyPDAF.PDAF.omi_put_state`\n    " \
                                    "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`,\n    " \
                                    "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`,\n    " \
                                    "   and :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`." \
                                    "\n\n    " \
                                    "References\n    " \
                                    "----------\n    " \
                                    ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                    "       A unification of ensemble square root Kalman filters. \n    " \
                                    "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                    ".. [2] Tödter, J., and B. Ahrens, 2015:\n    "\
                                    "       A second-order exact ensemble square root filter\n    " \
                                    "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                    "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                    ".. [3] Nerger, L.. (2022) \n    " \
                                    "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                    "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['omi_put_state_generate_obs'] = "Generation of synthetic observations based on given error statistics and observation operator for diagonal observation error covariance matrix\n    " \
                                           "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                           "If non-diagonal observation error covariance matrix has to be used,\n    " \
                                           "the generic :func:`pyPDAF.PDAF.put_generate_obs` can be used.\n\n    "\
                                           "The generated synthetic observations are based on each member of model forecast.\n    " \
                                           "Therefore, an ensemble of observations can be obtained. In a typical experiment,\n    "\
                                           "one may only need one ensemble member.\n\n    " \
                                           "Compared to :func:`pyPDAF.PDAF.omi_generate_obs`, this function has no :func:`get_state` call.\n    " \
                                           "This means that the next DA step will not be assigned by user-supplied functions.\n    " \
                                           "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                           "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                           "function call to ensure the sequential DA.\n\n    " \
                                           "The implementation strategy is similar to an assimilation step. This means that, \n    " \
                                           "one can reuse many user-supplied functions for assimilation and observation generation.\n\n    " \
                                           "This function executes the user-supplied function in the following sequence:\n    " \
                                           "    1. py__collect_state_pdaf\n    " \
                                           "    2. py__prepoststep_state_pdaf\n    " \
                                           "    3. py__init_dim_obs_pdaf\n    " \
                                           "    4. py__obs_op_pda\n    " \
                                           "    5. py__get_obs_f_pdaf"

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

docstrings['omi_assimilate_3dvar_nondiagR'] = "3DVar DA for a single DA step using non-diagnoal observation error covariance matrix.\n\n    "\
                                              "See :func:`pyPDAF.PDAF.omi_assimilate_3dvar` for simpler user-supplied functions\n    " \
                                              "using diagonal observation error covariance matrix.\n\n    "\
                                              "When 3DVar is used, the background error covariance matrix\n    "\
                                              "has to be modelled for cotrol variable transformation.\n    " \
                                              "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                              "This function should be called at each model time step.\n\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                              "User-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf\n    " \
                                              "    5. Iterative optimisation:\n    " \
                                              "        1. py__cvt_pdaf\n    " \
                                              "        2. py__obs_op_lin_pdaf\n    " \
                                              "        3. py__prodRinvA_pdaf\n    " \
                                              "        4. py__obs_op_adj_pdaf\n    " \
                                              "        5. py__cvt_adj_pdaf\n    " \
                                              "        6. core DA algorithm\n    " \
                                              "    6. py__cvt_pdaf\n    " \
                                              "    7. py__prepoststep_state_pdaf\n    " \
                                              "    8. py__distribute_state_pdaf\n    " \
                                              "    9. py__next_observation_pdaf"
docstrings['omi_assimilate_en3dvar_estkf_nondiagR'] = "3DEnVar for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                                      "See :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf` for simpler user-supplied functions\n    " \
                                                      "using diagonal observation error covariance matirx.\n\n    " \
                                                      "Here, the background error covariance matrix is estimated by an ensemble.\n    " \
                                                      "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                      "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                                      "This function should be called at each model time step.\n\n    " \
                                                      "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`\n    " \
                                                      "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                      "The user-supplied functions are executed in the following sequence:\n    " \
                                                      "    1. py__collect_state_pdaf\n    " \
                                                      "    2. py__prepoststep_state_pdaf\n    " \
                                                      "    3. py__init_dim_obs_pdaf\n    " \
                                                      "    4. py__obs_op_pdaf\n    " \
                                                      "    5. the iterative optimisation:\n    " \
                                                      "        1. py__cvt_ens_pdaf\n    " \
                                                      "        2. py__obs_op_lin_pdaf\n    " \
                                                      "        3. py__prodRinvA_pdaf\n    " \
                                                      "        4. py__obs_op_adj_pdaf\n    " \
                                                      "        5. py__cvt_adj_ens_pdaf\n    " \
                                                      "        6. core 3DEnVar algorithm\n    " \
                                                      "    6. py__cvt_ens_pdaf\n    " \
                                                      "    7. ESTKF:\n    " \
                                                      "        1. py__init_dim_obs_pdaf\n    " \
                                                      "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                                      "        3. py__obs_op_pdaf (for each ensemble member)\n    " \
                                                      "        4. py__prodRinvA_pdaf\n    " \
                                                      "        5. core ESTKF algorithm\n    " \
                                                      "    8. py__prepoststep_state_pdaf\n    " \
                                                      "    9. py__distribute_state_pdaf\n    " \
                                                      "    10. py__next_observation_pdaf"
docstrings['omi_assimilate_en3dvar_lestkf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`\n    "\
                                                       "or :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`.\n\n    "\
                                                       "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                       "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using non-diagnoal observation error covariance matrix.\n    " \
                                                       "The background error covariance matrix is estimated by ensemble.\n    " \
                                                       "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                       "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                       "This function should be called at each model time step.\n\n    " \
                                                       "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_en3dvar_lestkf_nondiagR`\n    " \
                                                       "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                       "The user-supplied function are executed in the following sequence:\n    " \
                                                       "    1. py__collect_state_pdaf\n    " \
                                                       "    2. py__prepoststep_state_pdaf\n    " \
                                                       "    3. py__init_dim_obs_pdaf\n    " \
                                                       "    4. py__obs_op_pdaf\n    " \
                                                       "    5. Starting the iterative optimisation:\n    " \
                                                       "        1. py__cvt_ens_pdaf\n    " \
                                                       "        2. py__obs_op_lin_pdaf\n    " \
                                                       "        3. py__prodRinvA_pdaf\n    " \
                                                       "        4. py__obs_op_adj_pdaf\n    " \
                                                       "        5. py__cvt_adj_ens_pdaf\n    " \
                                                       "        6. core DA algorithm\n    " \
                                                       "    6. py__cvt_ens_pdaf\n    " \
                                                       "    7. Perform LESTKF:\n    " \
                                                       "        1. py__init_n_domains_p_pdaf\n    " \
                                                       "        2. py__init_dim_obs_pdaf\n    " \
                                                       "        3. py__obs_op_pdaf\n    "\
                                                       "           (for each ensemble member)\n    " \
                                                       "        4. loop over each local domain:\n    " \
                                                       "            1. py__init_dim_l_pdaf\n    " \
                                                       "            2. py__init_dim_obs_l_pdaf\n    " \
                                                       "            3. py__g2l_state_pdaf\n    " \
                                                       "            4. py__prodRinvA_l_pdaf\n    " \
                                                       "            5. core DA algorithm\n    " \
                                                       "            6. py__l2g_state_pdaf\n    " \
                                                       "    8. py__prepoststep_state_pdaf\n    " \
                                                       "    9. py__distribute_state_pdaf\n    " \
                                                       "    10. py__next_observation_pdaf\n" \
                                                       "\n    " \
                                                       ".. deprecated:: 1.0.0\n\n    " \
                                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    " \
                                                       "   and :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`"
docstrings['omi_assimilate_hyb3dvar_estkf_nondiagR'] = "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix.\n\n    "\
                                                       "See :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf` for simpler user-supplied functions\n    " \
                                                       "using diagonal observation error covariance matrix.\n\n    " \
                                                       "Here the background error covariance is hybridised by a static background error covariance,\n    " \
                                                       "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                                       "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                       "ESTKF in this implementation.\n    " \
                                                       "This function should be called at each model time step.\n\n    " \
                                                       "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`\n    " \
                                                       "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                       "The user-supplied functions are executed in the following sequence:\n    " \
                                                       "    1. py__collect_state_pdaf\n    " \
                                                       "    2. py__prepoststep_state_pdaf\n    " \
                                                       "    3. py__init_dim_obs_pdaf\n    " \
                                                       "    4. py__obs_op_pdaf\n    " \
                                                       "    5. the iterative optimisation:\n    " \
                                                       "        1. py__cvt_pdaf\n    " \
                                                       "        2. py__cvt_ens_pdaf\n    " \
                                                       "        3. py__obs_op_lin_pdaf\n    " \
                                                       "        4. py__prodRinvA_pdaf\n    " \
                                                       "        5. py__obs_op_adj_pdaf\n    " \
                                                       "        6. py__cvt_adj_pdaf\n    " \
                                                       "        7. py__cvt_adj_ens_pdaf\n    " \
                                                       "        8. core 3DEnVar algorithm\n    " \
                                                       "    6. py__cvt_pdaf\n    " \
                                                       "    7. py__cvt_ens_pdaf\n    " \
                                                       "    8. Perform ESTKF:\n    " \
                                                       "        1. py__init_dim_obs_pdaf\n    " \
                                                       "        2. py__obs_op_pdaf\n    "\
                                                       "           (for ensemble mean)\n    " \
                                                       "        3. py__obs_op_pdaf\n    " \
                                                       "           (for each ensemble member)\n    " \
                                                       "        4. py__prodRinvA_pdaf\n    " \
                                                       "        5. core ESTKF algorithm\n    " \
                                                       "    9. py__prepoststep_state_pdaf\n    " \
                                                       "    10. py__distribute_state_pdaf\n    " \
                                                       "    11. py__next_observation_pdaf"
docstrings['omi_assimilate_hyb3dvar_lestkf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`\n    "\
                                                        "or :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`.\n\n    "\
                                                        "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                        "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix where\n    " \
                                                        "the background error covariance is hybridised by a static background error covariance,\n    " \
                                                        "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                                        "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                        "LESTKF in this implementation.\n    " \
                                                        "This function should be called at each model time step.\n\n    " \
                                                        "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_lestkf_nondiagR`\n    " \
                                                        "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                        "The user-supplied functions are executed in the following sequence:\n    " \
                                                        "    1. py__collect_state_pdaf\n    " \
                                                        "    2. py__prepoststep_state_pdaf\n    " \
                                                        "    3. py__init_dim_obs_pdaf\n    " \
                                                        "    4. py__obs_op_pdaf\n    " \
                                                        "    5. The iterative optimisation:\n    " \
                                                        "        1. py__cvt_pdaf\n    " \
                                                        "        2. py__cvt_ens_pdaf\n    " \
                                                        "        3. py__obs_op_lin_pdaf\n    " \
                                                        "        4. py__prodRinvA_pdaf\n    " \
                                                        "        5. py__obs_op_adj_pdaf\n    " \
                                                        "        6. py__cvt_adj_pdaf\n    " \
                                                        "        7. py__cvt_adj_ens_pdaf\n    " \
                                                        "        8. core DA algorithm\n    " \
                                                        "    6. py__cvt_pdaf\n    " \
                                                        "    7. py__cvt_ens_pdaf\n    " \
                                                        "    8. Perform LESTKF:\n    " \
                                                        "        1. py__init_n_domains_p_pdaf\n    " \
                                                        "        2. py__init_dim_obs_pdaf\n    " \
                                                        "        3. py__obs_op_pdaf\n    " \
                                                        "           (for each ensemble member)\n    " \
                                                        "        4. loop over each local domain:\n    " \
                                                        "            1. py__init_dim_l_pdaf\n    " \
                                                        "            2. py__init_dim_obs_l_pdaf\n    " \
                                                        "            3. py__g2l_state_pdaf\n    " \
                                                        "            4. py__prodRinvA_l_pdaf\n    " \
                                                        "            5. core DA algorithm\n    " \
                                                        "            6. py__l2g_state_pdaf\n    " \
                                                        "    9. py__prepoststep_state_pdaf\n    " \
                                                        "    10. py__distribute_state_pdaf\n    " \
                                                        "    11. py__next_observation_pdaf\n" \
                                                        "\n    " \
                                                        ".. deprecated:: 1.0.0\n\n    " \
                                                        "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    " \
                                                        "   and :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`"
docstrings['omi_assimilate_enkf_nondiagR'] = "Stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                             "See :func:`pyPDAF.PDAF.omi_assimilate_global` for simpler user-supplied functions\n    " \
                                             "using diagonal observation error covariance matrix.\n\n    " \
                                             "The stochastic EnKF is proposed by Evensen [1]_ and is a Monte Carlo approximation of the KF.\n\n    " \
                                             "This function should be called at each model time step. \n\n    " \
                                             "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR` " \
                                             "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                             "This function executes the user-supplied functions in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                             "    5. py__add_obs_err_pdaf\n    " \
                                             "    6. py__init_obscovar_pdaf\n    " \
                                             "    7. py__obs_op_pdaf (for each ensemble member)\n    " \
                                             "    8. core DA algorithm\n    " \
                                             "    9. py__prepoststep_state_pdaf\n    " \
                                             "    10. py__distribute_state_pdaf\n    " \
                                             "    11. py__next_observation_pdaf" \
                                             "\n\n    " \
                                             "References\n    " \
                                             "----------\n    " \
                                             ".. [1] Evensen, G. (1994), \n    "\
                                             "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                             "       using Monte Carlo methods to forecast error statistics,\n    "\
                                             "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572."
docstrings['omi_assimilate_global_nondiagR'] = "Global filters except for 3DVar and stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix.\n\n    "\
                                               "See :func:`pyPDAF.PDAF.omi_assimilate_global` for simpler user-supplied functions\n    " \
                                               "using diagonal observation error covariance matrix.\n\n    " \
                                               "Here, this function call is used for global, E(S)TKF [1]_, \n    " \
                                               "SEEK [1]_, SEIK [1]_.\n    " \
                                               "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                               "This function should be called at each model time step. \n\n    " \
                                               "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_global` " \
                                               "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                               "This function executes the user-supplied functions in the following sequence:\n    " \
                                               "    1. py__collect_state_pdaf\n    " \
                                               "    2. py__prepoststep_state_pdaf\n    " \
                                               "    3. py__init_dim_obs_pdaf\n    " \
                                               "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                               "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                               "    6. py__prodRinvA_pdaf\n    " \
                                               "    7. core DA algorithm\n    " \
                                               "    8. py__prepoststep_state_pdaf\n    " \
                                               "    9. py__distribute_state_pdaf\n    " \
                                               "    10. py__next_observation_pdaf" \
                                               "\n\n    " \
                                               "References\n    " \
                                               "----------\n    " \
                                               ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                               "       A unification of ensemble square root Kalman filters. \n    " \
                                               "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['omi_assimilate_nonlin_nondiagR'] = "Global nonlinear filters for a single DA step using non-diagnoal observation error covariance matrix.\n\n    "\
                                               "See :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR` for simpler user-supplied functions\n    " \
                                               "using diagonal observation error covariance matrix.\n\n    " \
                                               "Here, this function call is used for global NETF [1]_, and particle filter [2]_.\n    " \
                                               "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                               "This function should be called at each model time step. \n\n    " \
                                               "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR` " \
                                               "and :func:`pyPDAF.PDAF.get_state`.\n\n    "\
                                               "This function executes the user-supplied functions in the following sequence:\n    " \
                                               "    1. py__collect_state_pdaf\n    " \
                                               "    2. py__prepoststep_state_pdaf\n    " \
                                               "    3. py__init_dim_obs_pdaf\n    " \
                                               "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                               "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                               "    6. py__likelihood_pdaf\n    " \
                                               "    7. core DA algorithm\n    " \
                                               "    8. py__prepoststep_state_pdaf\n    " \
                                               "    9. py__distribute_state_pdaf\n    " \
                                               "    10. py__next_observation_pdaf" \
                                               "\n\n    " \
                                               "References\n    " \
                                               "----------\n    " \
                                               ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                               "       A second-order exact ensemble square root filter\n    " \
                                               "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                               "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                               ".. [2] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                                               "       Particle filters for high‐dimensional geoscience applications:\n    "\
                                               "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['omi_assimilate_lenkf_nondiagR'] = "Covariance localised stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                              "See :func:`pyPDAF.PDAF.omi_assimilate_lenkf` for simpler user-supplied functions\n    " \
                                              "using diagnoal observation error covariance matrix.\n\n    "\
                                              "This stochastic EnKF is implemented based on [1]_\n\n    " \
                                              "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                              "This function should be called at each model time step.\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                              "The user-supplied function is executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                              "    5. py__localize_pdaf\n    " \
                                              "    6. py__add_obs_err_pdaf\n    " \
                                              "    7. py__init_obscovar_pdaf\n    " \
                                              "    8. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                              "    9. core DA algorith\n    " \
                                              "    10. py__prepoststep_state_pdaf\n    " \
                                              "    11. py__distribute_state_pdaf\n    " \
                                              "    12. py__next_observation_pdaf" \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                              "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                              "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                              "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['omi_assimilate_local_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`\n    "\
                                              "or :func:`pyPDAF.PDAF.localomi_assimilate`.\n\n    " \
                                              "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                              "Domain local filters for a single DA step using non-diagnoal observation error covariance matrix.\n    " \
                                              "Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_\n    " \
                                              "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                              "This function should be called at each model time step.\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_local_nondiagR`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                              "User-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    "\
                                              "    2. py__prepoststep_state_pdaf\n    "\
                                              "    3. py__init_n_domains_p_pdaf\n    "\
                                              "    4. py__init_dim_obs_pdaf\n    "\
                                              "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                              "    6. loop over each local domain:\n    " \
                                              "        1. py__init_dim_l_pdaf\n    "\
                                              "        2. py__init_dim_obs_l_pdaf\n    "\
                                              "        3. py__g2l_state_pdaf\n    "\
                                              "        4. py__init_obs_l_pdaf\n    "\
                                              "        5. py__prodRinvA_l_pdaf\n    " \
                                              "        6. core DA algorithm\n    " \
                                              "        7. py__l2g_state_pdaf\n    "\
                                              "    7. py__prepoststep_state_pdaf\n    "\
                                              "    8. py__distribute_state_pdaf\n    "\
                                              "    9. py__next_observation_pdaf\n" \
                                              "\n    " \
                                              ".. deprecated:: 1.0.0\n\n    " \
                                              "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                              "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`." \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                              "       A unification of ensemble square root Kalman filters. \n    " \
                                              "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['omi_assimilate_lnetf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`\n    "\
                                              "or :func:`pyPDAF.PDAF.localomi_assimilate`.\n\n    " \
                                              "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                              "LNETF [1]_ for a single DA step using non-diagnoal observation error covariance matrix.\n    " \
                                              "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                              "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                              "This function should be called at each model time step.\n    " \
                                              "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_lnetf_nondiagR`\n    " \
                                              "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                              "User-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    "\
                                              "    2. py__prepoststep_state_pdaf\n    "\
                                              "    3. py__init_n_domains_p_pdaf\n    "\
                                              "    4. py__init_dim_obs_pdaf\n    "\
                                              "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                              "    6. loop over each local domain:\n    " \
                                              "        1. py__init_dim_l_pdaf\n    "\
                                              "        2. py__init_dim_obs_l_pdaf\n    "\
                                              "        3. py__g2l_state_pdaf\n    "\
                                              "        4. py__likelihood_l_pdaf\n    " \
                                              "        5. core DA algorithm\n    " \
                                              "        6. py__l2g_state_pdaf\n    "\
                                              "    7. py__prepoststep_state_pdaf\n    "\
                                              "    8. py__distribute_state_pdaf\n    "\
                                              "    9. py__next_observation_pdaf\n" \
                                              "\n    " \
                                              ".. deprecated:: 1.0.0\n\n    " \
                                              "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                              "   and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`." \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                              "       A second-order exact ensemble square root filter\n    " \
                                              "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                              "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['omi_assimilate_lknetf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`\n    "\
                                               "or :func:`pyPDAF.PDAF.localomi_assimilate`.\n\n    " \
                                               "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                               "LKNETF [1]_ for a single DA step using non-diagnoal observation error covariance matrix.\n    " \
                                               "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                               "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                               "This function should be called at each model time step.\n    " \
                                               "The function is a combination of :func:`pyPDAF.PDAF.omi_put_state_lknetf_nondiagR`\n    " \
                                               "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                               "User-supplied functions are executed in the following sequence:\n    " \
                                               "    1. py__collect_state_pdaf\n    "\
                                               "    2. py__prepoststep_state_pdaf\n    "\
                                               "    3. py__init_n_domains_p_pdaf\n    "\
                                               "    4. py__init_dim_obs_pdaf\n    "\
                                               "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                               "    6. loop over each local domain:\n    " \
                                               "        1. py__init_dim_l_pdaf\n    "\
                                               "        2. py__init_dim_obs_l_pdaf\n    "\
                                               "        3. py__g2l_state_pdaf\n    "\
                                               "        4. py__prodRinvA_pdaf\n    " \
                                               "        5. py__likelihood_l_pdaf\n    " \
                                               "        6. core DA algorithm\n    " \
                                               "        7. py__l2g_state_pdaf\n    "\
                                               "        8. py__obs_op_pdaf\n    " \
                                               "           (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                               "        9. py__likelihood_hyb_l_pda\n    " \
                                               "        10. py__prodRinvA_hyb_l_pdaf\n    " \
                                               "    7. py__prepoststep_state_pdaf\n    "\
                                               "    8. py__distribute_state_pdaf\n    "\
                                               "    9. py__next_observation_pdaf\n" \
                                               "\n    " \
                                               ".. deprecated:: 1.0.0\n\n    " \
                                               "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                               "   and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`." \
                                               "\n\n    " \
                                               "References\n    " \
                                               "----------\n    " \
                                               ".. [1] Nerger, L.. (2022) \n    " \
                                               "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                               "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

docstrings['omi_put_state_3dvar_nondiagR'] = "3DVar DA for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                             "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                             "See :func:`pyPDAF.PDAF.omi_put_state_3dvar` for simpler user-supplied functions using\n    " \
                                             "diagonal observation error covariance matrix.\n\n    "\
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "When 3DVar is used, the background error covariance matrix\n    "\
                                             "has to be modelled for cotrol variable transformation.\n    " \
                                             "This is a deterministic filtering scheme so no ensemble and parallelisation is needed.\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "User-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf\n    " \
                                             "    5. Iterative optimisation:\n    " \
                                             "        1. py__cvt_pdaf\n    " \
                                             "        2. py__obs_op_lin_pdaf\n    " \
                                             "        3. py__prodRinvA_pdaf\n    " \
                                             "        4. py__obs_op_adj_pdaf\n    " \
                                             "        5. py__cvt_adj_pdaf\n    " \
                                             "        6. core DA algorithm\n    " \
                                             "    6. py__cvt_pdaf"
docstrings['omi_put_state_en3dvar_estkf_nondiagR'] = "3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                     "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                     "See :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf` for simpler user-supplied functions\n    " \
                                                     "using diagonal observation error covariance matrix.\n\n    " \
                                                     "Compared to :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                     "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                     "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                     "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                     "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                     "function call to ensure the sequential DA.\n\n    " \
                                                     "The background error covariance matrix is estimated by an ensemble.\n    " \
                                                     "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                     "An ESTKF is used along with 3DEnVar to generate ensemble perturbations.\n    " \
                                                     "This function should be called at each model time step.\n\n    " \
                                                     "The user-supplied functions are executed in the following sequence:\n    " \
                                                     "    1. py__collect_state_pdaf\n    " \
                                                     "    2. py__prepoststep_state_pdaf\n    " \
                                                     "    3. py__init_dim_obs_pdaf\n    " \
                                                     "    4. py__obs_op_pdaf\n    " \
                                                     "    5. the iterative optimisation:\n    " \
                                                     "        1. py__cvt_ens_pdaf\n    " \
                                                     "        2. py__obs_op_lin_pdaf\n    " \
                                                     "        3. py__prodRinvA_pdaf\n    " \
                                                     "        4. py__obs_op_adj_pdaf\n    " \
                                                     "        5. py__cvt_adj_ens_pdaf\n    " \
                                                     "        6. core 3DEnVar algorithm\n    " \
                                                     "    6. py__cvt_ens_pdaf\n    " \
                                                     "    7. ESTKF:\n    " \
                                                     "        1. py__init_dim_obs_pdaf\n    " \
                                                     "        2. py__obs_op_pdaf (for ensemble mean)\n    " \
                                                     "        3. py__obs_op_pdaf (for each ensemble member)\n    " \
                                                     "        4. py__prodRinvA_pdaf\n    " \
                                                     "        5. core ESTKF algorithm"
docstrings['omi_put_state_en3dvar_lestkf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`\n    "\
                                                      "or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`.\n\n    "\
                                                      "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                      "3DEnVar for a single DA step without post-processing, distributing analysis, and setting next observation step,\n     "\
                                                      "where the ensemble anomaly is generated by LESTKF using non-diagnoal observation error covariance matrix.\n\n    " \
                                                      "Compared to :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_lestkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                      "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                      "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                      "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                      "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                      "function call to ensure the sequential DA.\n\n    " \
                                                      "The background error covariance matrix is estimated by ensemble.\n    " \
                                                      "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                      "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                      "This function should be called at each model time step.\n\n    " \
                                                      "The user-supplied function are executed in the following sequence:\n    " \
                                                      "    1. py__collect_state_pdaf\n    " \
                                                      "    2. py__prepoststep_state_pdaf\n    " \
                                                      "    3. py__init_dim_obs_pdaf\n    " \
                                                      "    4. py__obs_op_pdaf\n    " \
                                                      "    5. Starting the iterative optimisation:\n    " \
                                                      "        1. py__cvt_ens_pdaf\n    " \
                                                      "        2. py__obs_op_lin_pdaf\n    " \
                                                      "        3. py__prodRinvA_pdaf\n    " \
                                                      "        4. py__obs_op_adj_pdaf\n    " \
                                                      "        5. py__cvt_adj_ens_pdaf\n    " \
                                                      "        6. core DA algorithm\n    " \
                                                      "    6. py__cvt_ens_pdaf\n    " \
                                                      "    7. Perform LESTKF:\n    " \
                                                      "        1. py__init_n_domains_p_pdaf\n    " \
                                                      "        2. py__init_dim_obs_pdaf\n    " \
                                                      "        3. py__obs_op_pdaf\n    "\
                                                      "           (for each ensemble member)\n    " \
                                                      "        4. loop over each local domain:\n    " \
                                                      "            1. py__init_dim_l_pdaf\n    " \
                                                      "            2. py__init_dim_obs_l_pdaf\n    " \
                                                      "            3. py__g2l_state_pdaf\n    " \
                                                      "            4. py__prodRinvA_l_pdaf\n    " \
                                                      "            5. core DA algorithm\n    " \
                                                      "            6. py__l2g_state_pdaf\n" \
                                                      "\n    " \
                                                      ".. deprecated:: 1.0.0\n\n    " \
                                                      "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`\n    " \
                                                      "   and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`"
docstrings['omi_put_state_hyb3dvar_estkf_nondiagR'] = "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    "\
                                                      "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                      "See :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf` for simpler user-supplied functions\n    " \
                                                      "using diagonal observation error covariance matrix.\n\n    " \
                                                      "Here the background error covariance is hybridised by a static background error covariance,\n    " \
                                                      "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                                      "Compared to :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                      "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                      "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                      "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                      "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                      "function call to ensure the sequential DA.\n\n    " \
                                                      "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                      "ESTKF in this implementation.\n    " \
                                                      "This function should be called at each model time step.\n\n    " \
                                                      "The user-supplied functions are executed in the following sequence:\n    " \
                                                      "    1. py__collect_state_pdaf\n    " \
                                                      "    2. py__prepoststep_state_pdaf\n    " \
                                                      "    3. py__init_dim_obs_pdaf\n    " \
                                                      "    4. py__obs_op_pdaf\n    " \
                                                      "    5. the iterative optimisation:\n    " \
                                                      "        1. py__cvt_pdaf\n    " \
                                                      "        2. py__cvt_ens_pdaf\n    " \
                                                      "        3. py__obs_op_lin_pdaf\n    " \
                                                      "        4. py__prodRinvA_pdaf\n    " \
                                                      "        5. py__obs_op_adj_pdaf\n    " \
                                                      "        6. py__cvt_adj_pdaf\n    " \
                                                      "        7. py__cvt_adj_ens_pdaf\n    " \
                                                      "        8. core 3DEnVar algorithm\n    " \
                                                      "    6. py__cvt_pdaf\n    " \
                                                      "    7. py__cvt_ens_pdaf\n    " \
                                                      "    8. Perform ESTKF:\n    " \
                                                      "        1. py__init_dim_obs_pdaf\n    " \
                                                      "        2. py__obs_op_pdaf\n    "\
                                                      "           (for ensemble mean)\n    " \
                                                      "        3. py__obs_op_pdaf\n    " \
                                                      "           (for each ensemble member)\n    " \
                                                      "        4. py__prodRinvA_pdaf\n    " \
                                                      "        5. core ESTKF algorithm"
docstrings['omi_put_state_hyb3dvar_lestkf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`\n    "\
                                                       "or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`.\n\n    "\
                                                       "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                       "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                       "without post-processing, distributing analysis, and setting next observation step, where\n    " \
                                                       "the background error covariance is hybridised by a static background error covariance,\n    " \
                                                       "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                                       "Compared to :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_lestkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                       "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                       "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                       "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                       "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                       "function call to ensure the sequential DA.\n\n    " \
                                                       "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                       "LESTKF in this implementation.\n    " \
                                                       "This function should be called at each model time step.\n\n    " \
                                                       "The user-supplied functions are executed in the following sequence:\n    " \
                                                       "    1. py__collect_state_pdaf\n    " \
                                                       "    2. py__prepoststep_state_pdaf\n    " \
                                                       "    3. py__init_dim_obs_pdaf\n    " \
                                                       "    4. py__obs_op_pdaf\n    " \
                                                       "    5. The iterative optimisation:\n    " \
                                                       "        1. py__cvt_pdaf\n    " \
                                                       "        2. py__cvt_ens_pdaf\n    " \
                                                       "        3. py__obs_op_lin_pdaf\n    " \
                                                       "        4. py__prodRinvA_pdaf\n    " \
                                                       "        5. py__obs_op_adj_pdaf\n    " \
                                                       "        6. py__cvt_adj_pdaf\n    " \
                                                       "        7. py__cvt_adj_ens_pdaf\n    " \
                                                       "        8. core DA algorithm\n    " \
                                                       "    6. py__cvt_pdaf\n    " \
                                                       "    7. py__cvt_ens_pdaf\n    " \
                                                       "    8. Perform LESTKF:\n    " \
                                                       "        1. py__init_n_domains_p_pdaf\n    " \
                                                       "        2. py__init_dim_obs_pdaf\n    " \
                                                       "        3. py__obs_op_pdaf\n    " \
                                                       "           (for each ensemble member)\n    " \
                                                       "        4. loop over each local domain:\n    " \
                                                       "            1. py__init_dim_l_pdaf\n    " \
                                                       "            2. py__init_dim_obs_l_pdaf\n    " \
                                                       "            3. py__g2l_state_pdaf\n    " \
                                                       "            4. py__prodRinvA_l_pdaf\n    " \
                                                       "            5. core DA algorithm\n    " \
                                                       "            6. py__l2g_state_pdaf\n" \
                                                       "\n    " \
                                                       ".. deprecated:: 1.0.0\n\n    " \
                                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                                       "   and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`"
docstrings['omi_put_state_enkf_nondiagR'] = "Stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                            "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                            "See :func:`pyPDAF.PDAF.omi_put_state_global` for simpler user-supplied functions\n    " \
                                            "using diagonal observation error covariance matrix.\n\n    " \
                                            "The stochastic EnKF is implemented based on [1]_.\n\n    " \
                                            "Compared to :func:`pyPDAF.PDAF.omi_assimilate_enkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                            "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                            "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                            "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                            "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                            "function call to ensure the sequential DA.\n\n    " \
                                            "This function should be called at each model time step. \n\n    " \
                                            "This function executes the user-supplied functions in the following sequence:\n    " \
                                            "    1. py__collect_state_pdaf\n    " \
                                            "    2. py__prepoststep_state_pdaf\n    " \
                                            "    3. py__init_dim_obs_pdaf\n    " \
                                            "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                            "    5. py__add_obs_err_pdaf\n    " \
                                            "    6. py__init_obscovar_pdaf\n    " \
                                            "    7. py__obs_op_pdaf (for each ensemble member)\n    " \
                                            "    8. core DA algorithm" \
                                            "\n\n    " \
                                            "References\n    " \
                                            "----------\n    " \
                                            ".. [1] Evensen, G. (1994), \n    "\
                                            "       Sequential data assimilation with a nonlinear quasi-geostrophic model\n    "\
                                            "       using Monte Carlo methods to forecast error statistics,\n    "\
                                            "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572."
docstrings['omi_put_state_global_nondiagR'] = "Global filters except for 3DVar and stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix\n    "\
                                              "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                              "See :func:`pyPDAF.PDAF.omi_put_state_global` for simpler user-supplied functions\n    " \
                                              "using diagonal observation error covariance matrix.\n\n    " \
                                              "Here, this function call is used for global, E(S)TKF [1]_, \n    " \
                                              "SEEK [1]_, SEIK [1]_.\n    " \
                                              "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                              "Compared to :func:`pyPDAF.PDAF.omi_assimilate_global_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                              "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                              "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                              "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                              "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                              "function call to ensure the sequential DA.\n\n    " \
                                              "This function should be called at each model time step. \n\n    " \
                                              "This function executes the user-supplied functions in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                              "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                              "    6. py__prodRinvA_pdaf\n    " \
                                              "    7. core DA algorithm" \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                              "       A unification of ensemble square root Kalman filters. \n    " \
                                              "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['omi_put_state_nonlin_nondiagR'] = "Global nonlinear filters for a single DA step using non-diagnoal observation error covariance matrix\n    "\
                                              "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                              "See :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR` for simpler user-supplied functions\n    " \
                                              "using diagonal observation error covariance matrix.\n\n    " \
                                              "Here, this function call is used for global NETF [1]_, and particle filter [2]_.\n    " \
                                              "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                              "Compared to :func:`pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                              "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                              "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                              "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                              "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                              "function call to ensure the sequential DA.\n\n    " \
                                              "This function should be called at each model time step. \n\n    " \
                                              "This function executes the user-supplied functions in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    " \
                                              "    2. py__prepoststep_state_pdaf\n    " \
                                              "    3. py__init_dim_obs_pdaf\n    " \
                                              "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
                                              "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                              "    6. py__likelihood_pdaf\n    " \
                                              "    7. core DA algorithm" \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                              "       A second-order exact ensemble square root filter\n    " \
                                              "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                              "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                              ".. [2] Van Leeuwen, P. J., Künsch, H. R., Nerger, L., Potthast, R., & Reich, S. (2019).\n    "\
                                              "       Particle filters for high‐dimensional geoscience applications:\n    "\
                                              "       A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365."
docstrings['omi_put_state_lenkf_nondiagR'] = "Covariance localised stochastic EnKF for a single DA step using non-diagnoal observation error covariance matrix\n    "\
                                             "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                             "See :func:`pyPDAF.PDAF.omi_put_state_lenkf` for simpler user-supplied functions\n    " \
                                             "using diagnoal observation error covariance matrix.\n\n    "\
                                             "This function is implemented based on [1]_.\n\n    " \
                                             "This is the only scheme for covariance localisation in PDAF.\n\n    " \
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                             "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "The user-supplied function is executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    " \
                                             "    2. py__prepoststep_state_pdaf\n    " \
                                             "    3. py__init_dim_obs_pdaf\n    " \
                                             "    4. py__obs_op_pdaf (for each ensemble member)\n    " \
                                             "    5. py__localize_pdaf\n    " \
                                             "    6. py__add_obs_err_pdaf\n    " \
                                             "    7. py__init_obscovar_pdaf\n    " \
                                             "    8. py__obs_op_pdaf (repeated to reduce storage)\n    " \
                                             "    9. core DA algorithm" \
                                             "\n\n    " \
                                             "References\n    " \
                                             "----------\n    " \
                                             ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
                                             "       Data Assimilation Using an Ensemble Kalman Filter Technique.\n    "\
                                             "       Mon. Wea. Rev., 126, 796–811,\n    "\
                                             "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['omi_put_state_local_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`\n    "\
                                             "or :func:`pyPDAF.PDAF.localomi_put_state`.\n\n    " \
                                             "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                             "Domain local filters for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                             "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                             "Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_\n    " \
                                             "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_local_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                             "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "User-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    "\
                                             "    2. py__prepoststep_state_pdaf\n    "\
                                             "    3. py__init_n_domains_p_pdaf\n    "\
                                             "    4. py__init_dim_obs_pdaf\n    "\
                                             "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                             "    6. loop over each local domain:\n    " \
                                             "        1. py__init_dim_l_pdaf\n    "\
                                             "        2. py__init_dim_obs_l_pdaf\n    "\
                                             "        3. py__g2l_state_pdaf\n    "\
                                             "        4. py__init_obs_l_pdaf\n    "\
                                             "        5. py__prodRinvA_l_pdaf\n    " \
                                             "        6. core DA algorithm\n    " \
                                             "        7. py__l2g_state_pdaf\n" \
                                             "\n    " \
                                             ".. deprecated:: 1.0.0\n\n    " \
                                             "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                             "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`." \
                                             "\n\n    " \
                                             "References\n    " \
                                             "----------\n    " \
                                             ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                             "       A unification of ensemble square root Kalman filters. \n    " \
                                             "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['omi_put_state_lnetf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`\n    "\
                                             "or :func:`pyPDAF.PDAF.localomi_put_state`.\n\n    " \
                                             "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                             "LNETF [1]_ for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                             "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                             "See :func:`pyPDAF.PDAF.localomi_put_state` for using diagnoal observation error covariance matrix.\n    " \
                                             "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                             "Compared to :func:`pyPDAF.PDAF.omi_assimilate_lnetf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                             "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                             "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                             "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                             "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                             "function call to ensure the sequential DA.\n\n    " \
                                             "This function should be called at each model time step.\n\n    " \
                                             "User-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    "\
                                             "    2. py__prepoststep_state_pdaf\n    "\
                                             "    3. py__init_n_domains_p_pdaf\n    "\
                                             "    4. py__init_dim_obs_pdaf\n    "\
                                             "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                             "    6. loop over each local domain:\n    " \
                                             "        1. py__init_dim_l_pdaf\n    "\
                                             "        2. py__init_dim_obs_l_pdaf\n    "\
                                             "        3. py__g2l_state_pdaf\n    "\
                                             "        4. py__likelihood_l_pdaf\n    " \
                                             "        5. core DA algorithm\n    " \
                                             "        6. py__l2g_state_pdaf\n" \
                                             "\n    " \
                                             ".. deprecated:: 1.0.0\n\n    " \
                                             "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                             "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`." \
                                             "\n\n    " \
                                             "References\n    " \
                                             "----------\n    " \
                                             ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                             "       A second-order exact ensemble square root filter\n    " \
                                             "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                             "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['omi_put_state_lknetf_nondiagR'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`\n    "\
                                              "or :func:`pyPDAF.PDAF.localomi_put_state`.\n\n    " \
                                              "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                              "LKNETF [1]_ for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                              "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                              "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                              "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                              "Compared to :func:`pyPDAF.PDAF.omi_assimilate_lknetf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                              "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                              "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                              "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                              "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                              "function call to ensure the sequential DA.\n\n    " \
                                              "This function should be called at each model time step.\n\n    " \
                                              "User-supplied functions are executed in the following sequence:\n    " \
                                              "    1. py__collect_state_pdaf\n    "\
                                              "    2. py__prepoststep_state_pdaf\n    "\
                                              "    3. py__init_n_domains_p_pdaf\n    "\
                                              "    4. py__init_dim_obs_pdaf\n    "\
                                              "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                              "    6. loop over each local domain:\n    " \
                                              "        1. py__init_dim_l_pdaf\n    "\
                                              "        2. py__init_dim_obs_l_pdaf\n    "\
                                              "        3. py__g2l_state_pdaf\n    "\
                                              "        4. py__prodRinvA_pdaf\n    " \
                                              "        5. py__likelihood_l_pdaf\n    " \
                                              "        6. core DA algorithm\n    " \
                                              "        7. py__l2g_state_pdaf\n    "\
                                              "        8. py__obs_op_pdaf\n    " \
                                              "           (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                              "        9. py__likelihood_hyb_l_pda\n    " \
                                              "        10. py__prodRinvA_hyb_l_pdaf\n" \
                                              "\n    " \
                                              ".. deprecated:: 1.0.0\n\n    " \
                                              "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                              "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`." \
                                              "\n\n    " \
                                              "References\n    " \
                                              "----------\n    " \
                                              ".. [1] Nerger, L.. (2022) \n    " \
                                              "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter.\n    " \
                                              "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

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

docstrings['local_assimilate_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    "\
                                                "or :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`.\n\n    "\
                                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF.\n    " \
                                                "The background error covariance matrix is estimated by ensemble.\n    " \
                                                "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                "This function should be called at each model time step.\n\n    " \
                                                "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_en3dvar_lestkf`\n    " \
                                                "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                "The user-supplied function are executed in the following sequence:\n    " \
                                                "    1. py__collect_state_pdaf\n    " \
                                                "    2. py__prepoststep_state_pdaf\n    " \
                                                "    3. py__init_dim_obs_pdaf\n    " \
                                                "    4. py__obs_op_pdaf\n    " \
                                                "    5. py__init_obs_pdaf\n    " \
                                                "    6. Starting the iterative optimisation:\n    " \
                                                "        1. py__cvt_ens_pdaf\n    " \
                                                "        2. py__obs_op_lin_pdaf\n    " \
                                                "        3. py__prodRinvA_pdaf\n    " \
                                                "        4. py__obs_op_adj_pdaf\n    " \
                                                "        5. py__cvt_adj_ens_pdaf\n    " \
                                                "        6. core DA algorithm\n    " \
                                                "    7. py__cvt_ens_pdaf\n    " \
                                                "    8. Perform LESTKF:\n    " \
                                                "        1. py__init_n_domains_p_pdaf\n    " \
                                                "        2. py__init_dim_obs_pdaf\n    " \
                                                "        3. py__obs_op_pdaf\n    "\
                                                "           (for each ensemble member)\n    " \
                                                "        4. py__init_obs_pdaf\n    "\
                                                "           (if global adaptive forgetting factor is used\n    "\
                                                "           `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                                "        5. py__init_obsvar_pdaf\n    "\
                                                "           (if global adaptive forgetting factor is used)\n    " \
                                                "        6. loop over each local domain:\n    " \
                                                "            1. py__init_dim_l_pdaf\n    " \
                                                "            2. py__init_dim_obs_l_pdaf\n    " \
                                                "            3. py__g2l_obs_pdaf\n    "\
                                                "               (localise mean ensemble in observation space)\n    " \
                                                "            4. py__init_obs_l_pdaf\n    "\
                                                "            5. py__g2l_obs_pdaf\n    " \
                                                "               (localise each ensemble member in observation space)\n    " \
                                                "            6. py__init_obsvar_l_pdaf\n    " \
                                                "               (only called if local adaptive forgetting factor\n    " \
                                                "               `type_forget=2` is used)\n    "\
                                                "            7. py__prodRinvA_l_pdaf\n    " \
                                                "            8. core DA algorithm\n    " \
                                                "    9. py__prepoststep_state_pdaf\n    " \
                                                "    10. py__distribute_state_pdaf\n    " \
                                                "    11. py__next_observation_pdaf\n" \
                                                "\n    " \
                                                ".. deprecated:: 1.0.0\n\n    " \
                                                "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`\n    " \
                                                "   and :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`"
docstrings['local_assimilate_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    "\
                                                 "or :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                                 "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                 "Hybrid 3DEnVar for a single DA step where\n    " \
                                                 "the background error covariance is hybridised by a static background error covariance,\n    " \
                                                 "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                                 "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                 "LESTKF in this implementation.\n    " \
                                                 "This function should be called at each model time step.\n\n    " \
                                                 "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_hyb3dvar_lestkf`\n    " \
                                                 "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                 "The user-supplied functions are executed in the following sequence:\n    " \
                                                 "    1. py__collect_state_pdaf\n    " \
                                                 "    2. py__prepoststep_state_pdaf\n    " \
                                                 "    3. py__init_dim_obs_pdaf\n    " \
                                                 "    4. py__obs_op_pdaf\n    " \
                                                 "    5. py__init_obs_pdaf\n    " \
                                                 "    6. The iterative optimisation:\n    " \
                                                 "        1. py__cvt_pdaf\n    " \
                                                 "        2. py__cvt_ens_pdaf\n    " \
                                                 "        3. py__obs_op_lin_pdaf\n    " \
                                                 "        4. py__prodRinvA_pdaf\n    " \
                                                 "        5. py__obs_op_adj_pdaf\n    " \
                                                 "        6. py__cvt_adj_pdaf\n    " \
                                                 "        7. py__cvt_adj_ens_pdaf\n    " \
                                                 "        8. core DA algorithm\n    " \
                                                 "    7. py__cvt_pdaf\n    " \
                                                 "    8. py__cvt_ens_pdaf\n    " \
                                                 "    9. Perform LESTKF:\n    " \
                                                 "        1. py__init_n_domains_p_pdaf\n    " \
                                                 "        2. py__init_dim_obs_pdaf\n    " \
                                                 "        3. py__obs_op_pdaf\n    " \
                                                 "           (for each ensemble member)\n    " \
                                                 "        4. py__init_obs_pdaf\n    " \
                                                 "           (if global adaptive forgetting factor `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                                 "        5. py__init_obsvar_pdaf\n    " \
                                                 "           (if global adaptive forgetting factor is used)\n    " \
                                                 "        6. loop over each local domain:\n    " \
                                                 "            1. py__init_dim_l_pdaf\n    " \
                                                 "            2. py__init_dim_obs_l_pdaf\n    " \
                                                 "            3. py__g2l_obs_pdaf\n    "\
                                                 "               (localise mean ensemble in observation space)\n    " \
                                                 "            4. py__init_obs_l_pdaf\n    "\
                                                 "            5. py__g2l_obs_pdaf\n    " \
                                                 "               (localise each ensemble member in observation space)\n    " \
                                                 "            6. py__init_obsvar_l_pdaf\n    " \
                                                 "               (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                                 "            7. py__prodRinvA_l_pdaf\n    " \
                                                 "            8. core DA algorithm\n    " \
                                                 "    10. py__prepoststep_state_pdaf\n    " \
                                                 "    11. py__distribute_state_pdaf\n    " \
                                                 "    12. py__next_observation_pdaf\n" \
                                                 "\n    " \
                                                 ".. deprecated:: 1.0.0\n\n    " \
                                                 "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`\n    " \
                                                 "   and :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`"
docstrings['local_assimilate_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                        "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                        "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                        "Local ESTKF (error space transform " \
                                        "Kalman filter) [1]_ for a single DA step without OMI.\n    " \
                                        "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                        "This function should be called at each model time step.\n    " \
                                        "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_lestkf`\n    " \
                                        "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                        "User-supplied functions are executed in the following sequence:\n    " \
                                        "    1. py__collect_state_pdaf\n    "\
                                        "    2. py__prepoststep_state_pdaf\n    "\
                                        "    3. py__init_n_domains_p_pdaf\n    "\
                                        "    4. py__init_dim_obs_pdaf\n    "\
                                        "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                        "    6. py__init_obs_pdaf\n    " \
                                        "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                        "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                        "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                        "    8. loop over each local domain:\n    " \
                                        "        1. py__init_dim_l_pdaf\n    "\
                                        "        2. py__init_dim_obs_l_pdaf\n    "\
                                        "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                        "        4. py__init_obs_l_pdaf\n    " \
                                        "        5. py__g2l_obs_pdaf\n    "\
                                        "           (localise each ensemble member in observation space)\n    "\
                                        "        6. py__init_obsvar_l_pdaf\n    " \
                                        "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                        "        7. py__prodRinvA_l_pdaf\n    "\
                                        "        8. core DA algorithm\n    " \
                                        "    9. py__prepoststep_state_pdaf\n    "\
                                        "    10. py__distribute_state_pdaf\n    "\
                                        "    11. py__next_observation_pdaf\n" \
                                        "\n    " \
                                        ".. deprecated:: 1.0.0\n\n    " \
                                        "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                        "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                        "\n\n    " \
                                        "References\n    " \
                                        "----------\n    " \
                                        ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                        "       A unification of ensemble square root Kalman filters. \n    " \
                                        "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['local_assimilate_letkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                       "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                       "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "Local ensemble transform Kalman filter (LETKF) [1]_ " \
                                       "for a single DA step without OMI. Implementation is based on [2]_.\n    " \
                                       "Note that the LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                       "This function should be called at each model time step.\n    " \
                                       "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_letkf`\n    " \
                                       "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                       "User-supplied functions are executed in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    "\
                                       "    2. py__prepoststep_state_pdaf\n    "\
                                       "    3. py__init_n_domains_p_pdaf\n    "\
                                       "    4. py__init_dim_obs_pdaf\n    "\
                                       "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                       "    6. py__init_obs_pdaf\n    " \
                                       "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                       "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                       "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                       "    8. loop over each local domain:\n    " \
                                       "        1. py__init_dim_l_pdaf\n    "\
                                       "        2. py__init_dim_obs_l_pdaf\n    "\
                                       "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                       "        4. py__init_obs_l_pdaf\n    " \
                                       "        5. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    "\
                                       "        6. py__init_obsvar_l_pdaf\n    " \
                                       "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                       "        7. py__prodRinvA_l_pdaf\n    "\
                                       "        8. core DA algorithm\n    " \
                                       "    9. py__prepoststep_state_pdaf\n    "\
                                       "    10. py__distribute_state_pdaf\n    "\
                                       "    11. py__next_observation_pdaf\n" \
                                       "\n    " \
                                       ".. deprecated:: 1.0.0\n\n    " \
                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                       "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007).\n    "\
                                       "       Efficient data assimilation for spatiotemporal chaos:\n    "\
                                       "       A local ensemble transform Kalman filter. \n    "\
                                       "       Physica D: Nonlinear Phenomena, 230(1-2), 112-126.\n    " \
                                       ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                       "       A unification of ensemble square root Kalman filters. \n    " \
                                       "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['local_assimilate_lseik'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                       "or :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`.\n\n    "\
                                       "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "Local singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n    " \
                                       "This function should be called at each model time step.\n\n    " \
                                       "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_lseik` " \
                                       "and :func:`pyPDAF.PDAF.get_state`\n\n    "\
                                       "This function  executes the user-supplied functions in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    " \
                                       "    2. py__prepoststep_state_pdaf\n    " \
                                       "    3. py__init_n_domains_p_pdaf\n    " \
                                       "    4. py__init_dim_obs_pdaf\n    " \
                                       "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                       "    6. py__init_obs_pdaf\n    "\
                                       "       (if global adaptive forgetting factor `type_forget=1`\n    " \
                                       "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                       "    7. py__init_obsvar_pdaf\n    "\
                                       "       (if global adaptive forgetting factor is used)\n    " \
                                       "    8. loop over each local domain:\n    " \
                                       "        1. py__init_dim_l_pdaf\n    " \
                                       "        2. py__init_dim_obs_l_pdaf\n    " \
                                       "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                       "        4. py__init_obs_l_pdaf\n    "\
                                       "        5. py__g2l_obs_pdaf\n    "\
                                       "           (localise each ensemble member in observation space)\n    " \
                                       "        6. py__init_obsvar_l_pdaf\n    "\
                                       "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                       "        7. py__prodRinvA_l_pdaf\n    " \
                                       "        8. core DA algorithm\n    " \
                                       "    9. py__prepoststep_state_pdaf\n    " \
                                       "    10. py__distribute_state_pdaf\n    " \
                                       "    11. py__next_observation_pdaf\n" \
                                       "\n    " \
                                       ".. deprecated:: 1.0.0\n\n    " \
                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                       "   and :func:`pyPDAF.PDAF.localomi_assimilate_nondiagR`" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                       "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                       "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['local_assimilate_lnetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                       "or :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`.\n\n    "\
                                       "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "Local Nonlinear Ensemble Transform Filter (LNETF) [1]_ for a single DA step.\n    " \
                                       "The nonlinear filter computes the distribution up to\n    " \
                                       "the second moment similar to Kalman filters but it uses a nonlinear weighting similar to\n    " \
                                       "particle filters. This leads to an equal weights assumption for the prior ensemble at each step.\n    " \
                                       "This function should be called at each model time step.\n\n    " \
                                       "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_lnetf`\n    " \
                                       "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                       "This function executes the user-supplied function in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    " \
                                       "    2. py__prepoststep_state_pdaf\n    " \
                                       "    3. py__init_n_domains_p_pdaf\n    " \
                                       "    4. py__init_dim_obs_pdaf\n    " \
                                       "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                       "    6. loop over each local domain:\n    " \
                                       "        1. py__init_dim_l_pdaf\n    " \
                                       "        2. py__init_dim_obs_l_pdaf\n    " \
                                       "        3. py__init_obs_l_pdaf\n    "\
                                       "        4. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                       "        5. py__likelihood_l_pdaf\n    " \
                                       "        6. core DA algorithm\n    " \
                                       "    7. py__prepoststep_state_pdaf\n    " \
                                       "    8. py__distribute_state_pdaf\n    " \
                                       "    9. py__next_observation_pdaf\n" \
                                       "\n    " \
                                       ".. deprecated:: 1.0.0\n\n    " \
                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                       "   and :func:`pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR`" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                       "       A second-order exact ensemble square root filter\n    " \
                                       "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                       "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['local_assimilate_lknetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_assimilate`\n    "\
                                        "or :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`.\n\n    "\
                                        "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                        "A hybridised LETKF and LNETF [1]_ for a single DA step.\n    " \
                                        "The LNETF computes the distribution up to\n    " \
                                        "the second moment similar to Kalman filters but using a nonlinear weighting similar to\n    " \
                                        "particle filters. This leads to an equal weights assumption for the prior ensemble.\n    " \
                                        "The hybridisation with LETKF is expected to lead to improved performance for\n    " \
                                        "quasi-Gaussian problems.\n    " \
                                        "The function should be called at each model step.\n\n    " \
                                        "The function is a combination of :func:`pyPDAF.PDAF.local_put_state_lknetf`\n    " \
                                        "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                        "This function executes the user-supplied function in the following sequence:\n    " \
                                        "    1. py__collect_state_pdaf\n    " \
                                        "    2. py__prepoststep_state_pdaf\n    " \
                                        "    3. py__init_n_domains_p_pdaf\n    " \
                                        "    4. py__init_dim_obs_pdaf\n    " \
                                        "    5. py__obs_op_pdaf\n    "\
                                        "       (for each ensemble member)\n    " \
                                        "    6. py__init_obs_pdaf\n    " \
                                        "       (if global adaptive forgetting factor `type_forget=1`\n    "\
                                        "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                        "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                        "    8. loop over each local domain:\n    " \
                                        "        1. py__init_dim_l_pdaf\n    " \
                                        "        2. py__init_dim_obs_l_pdaf\n    " \
                                        "        3. py__g2l_obs_pdaf\n    "\
                                        "           (localise each ensemble member in observation space)\n    " \
                                        "        4. py__init_obs_l_pdaf\n    "\
                                        "        5. py__init_obsvar_l_pdaf\n    "\
                                        "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                        "        6. py__prodRinvA_pdaf\n    " \
                                        "        7. py__likelihood_l_pdaf\n    " \
                                        "        8. core DA algorithm\n    " \
                                        "    9. py__obs_op_pdaf\n    " \
                                        "       (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                        "    10. py__likelihood_hyb_l_pda\n    " \
                                        "    11. py__init_obsvar_l_pdaf\n    " \
                                        "        (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                        "    12. py__prodRinvA_hyb_l_pdaf\n    " \
                                        "    13. py__prepoststep_state_pdaf\n    " \
                                        "    14. py__distribute_state_pdaf\n    " \
                                        "    15. py__next_observation_pdaf\n" \
                                        "\n    " \
                                        ".. deprecated:: 1.0.0\n\n    " \
                                        "   This function is replaced by :func:`pyPDAF.PDAF.localomi_assimilate`\n    " \
                                        "   and :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`" \
                                        "\n\n    " \
                                        "References\n    " \
                                        "----------\n    " \
                                        ".. [1] Nerger, L.. (2022) \n    " \
                                        "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                        "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

docstrings['local_put_state_en3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    "\
                                               "or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`.\n\n    "\
                                               "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                               "3DEnVar for a single DA step without post-processing, distributing analysis, and setting next observation step,\n     "\
                                               "where the ensemble anomaly is generated by LESTKF.\n\n    " \
                                               "Compared to :func:`pyPDAF.PDAF.local_assimilate_en3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                               "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                               "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                               "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                               "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                               "function call to ensure the sequential DA.\n\n    " \
                                               "The background error covariance matrix is estimated by ensemble.\n    " \
                                               "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                               "An LESTKF is used to generate ensemble perturbations.\n    " \
                                               "This function should be called at each model time step.\n\n    " \
                                               "The user-supplied function are executed in the following sequence:\n    " \
                                               "    1. py__collect_state_pdaf\n    " \
                                               "    2. py__prepoststep_state_pdaf\n    " \
                                               "    3. py__init_dim_obs_pdaf\n    " \
                                               "    4. py__obs_op_pdaf\n    " \
                                               "    5. py__init_obs_pdaf\n    " \
                                               "    6. Starting the iterative optimisation:\n    " \
                                               "        1. py__cvt_ens_pdaf\n    " \
                                               "        2. py__obs_op_lin_pdaf\n    " \
                                               "        3. py__prodRinvA_pdaf\n    " \
                                               "        4. py__obs_op_adj_pdaf\n    " \
                                               "        5. py__cvt_adj_ens_pdaf\n    " \
                                               "        6. core DA algorithm\n    " \
                                               "    7. py__cvt_ens_pdaf\n    " \
                                               "    8. Perform LESTKF:\n    " \
                                               "        1. py__init_n_domains_p_pdaf\n    " \
                                               "        2. py__init_dim_obs_pdaf\n    " \
                                               "        3. py__obs_op_pdaf\n    "\
                                               "           (for each ensemble member)\n    " \
                                               "        4. py__init_obs_pdaf\n    "\
                                               "           (if global adaptive forgetting factor is used\n    "\
                                               "           `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                               "        5. py__init_obsvar_pdaf\n    "\
                                               "           (if global adaptive forgetting factor is used)\n    " \
                                               "        6. loop over each local domain:\n    " \
                                               "            1. py__init_dim_l_pdaf\n    " \
                                               "            2. py__init_dim_obs_l_pdaf\n    " \
                                               "            3. py__g2l_obs_pdaf\n    "\
                                               "               (localise mean ensemble in observation space)\n    " \
                                               "            4. py__init_obs_l_pdaf\n    "\
                                               "            5. py__g2l_obs_pdaf\n    " \
                                               "               (localise each ensemble member in observation space)\n    " \
                                               "            6. py__init_obsvar_l_pdaf\n    " \
                                               "               (only called if local adaptive forgetting factor\n    " \
                                               "               `type_forget=2` is used)\n    "\
                                               "            7. py__prodRinvA_l_pdaf\n    " \
                                               "            8. core DA algorithm\n" \
                                               "\n    " \
                                               ".. deprecated:: 1.0.0\n\n    " \
                                               "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
                                               "   and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`"
docstrings['local_put_state_hyb3dvar_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    "\
                                                "or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`.\n\n    "\
                                                "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                                "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                "without post-processing, distributing analysis, and setting next observation step, where\n    " \
                                                "the background error covariance is hybridised by a static background error covariance,\n    " \
                                                "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                                "Compared to :func:`pyPDAF.PDAF.local_assimilate_hyb3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                                "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                "function call to ensure the sequential DA.\n\n    " \
                                                "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                "LESTKF in this implementation.\n    " \
                                                "This function should be called at each model time step.\n\n    " \
                                                "The user-supplied functions are executed in the following sequence:\n    " \
                                                "    1. py__collect_state_pdaf\n    " \
                                                "    2. py__prepoststep_state_pdaf\n    " \
                                                "    3. py__init_dim_obs_pdaf\n    " \
                                                "    4. py__obs_op_pdaf\n    " \
                                                "    5. py__init_obs_pdaf\n    " \
                                                "    6. The iterative optimisation:\n    " \
                                                "        1. py__cvt_pdaf\n    " \
                                                "        2. py__cvt_ens_pdaf\n    " \
                                                "        3. py__obs_op_lin_pdaf\n    " \
                                                "        4. py__prodRinvA_pdaf\n    " \
                                                "        5. py__obs_op_adj_pdaf\n    " \
                                                "        6. py__cvt_adj_pdaf\n    " \
                                                "        7. py__cvt_adj_ens_pdaf\n    " \
                                                "        8. core DA algorithm\n    " \
                                                "    7. py__cvt_pdaf\n    " \
                                                "    8. py__cvt_ens_pdaf\n    " \
                                                "    9. Perform LESTKF:\n    " \
                                                "        1. py__init_n_domains_p_pdaf\n    " \
                                                "        2. py__init_dim_obs_pdaf\n    " \
                                                "        3. py__obs_op_pdaf\n    " \
                                                "           (for each ensemble member)\n    " \
                                                "        4. py__init_obs_pdaf\n    " \
                                                "           (if global adaptive forgetting factor `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
                                                "        5. py__init_obsvar_pdaf\n    " \
                                                "           (if global adaptive forgetting factor is used)\n    " \
                                                "        6. loop over each local domain:\n    " \
                                                "            1. py__init_dim_l_pdaf\n    " \
                                                "            2. py__init_dim_obs_l_pdaf\n    " \
                                                "            3. py__g2l_obs_pdaf\n    "\
                                                "               (localise mean ensemble in observation space)\n    " \
                                                "            4. py__init_obs_l_pdaf\n    "\
                                                "            5. py__g2l_obs_pdaf\n    " \
                                                "               (localise each ensemble member in observation space)\n    " \
                                                "            6. py__init_obsvar_l_pdaf\n    " \
                                                "               (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                                "            7. py__prodRinvA_l_pdaf\n    " \
                                                "            8. core DA algorithm\n" \
                                                "\n    " \
                                                ".. deprecated:: 1.0.0\n\n    " \
                                                "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                                "   and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`"
docstrings['local_put_state_lestkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                       "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                       "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "Local ESTKF (error space transform " \
                                       "Kalman filter) [1]_ for a single DA step without OMI.\n\n    " \
                                       "Compared to :func:`pyPDAF.PDAF.local_assimilate_lestkf`, this function has no :func:`get_state` call.\n    " \
                                       "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                       "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                       "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                       "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                       "function call to ensure the sequential DA.\n\n    " \
                                       "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                       "This function should be called at each model time step.\n\n    " \
                                       "User-supplied functions are executed in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    "\
                                       "    2. py__prepoststep_state_pdaf\n    "\
                                       "    3. py__init_n_domains_p_pdaf\n    "\
                                       "    4. py__init_dim_obs_pdaf\n    "\
                                       "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                       "    6. py__init_obs_pdaf\n    " \
                                       "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                       "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                       "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                       "    8. loop over each local domain:\n    " \
                                       "        1. py__init_dim_l_pdaf\n    "\
                                       "        2. py__init_dim_obs_l_pdaf\n    "\
                                       "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                       "        4. py__init_obs_l_pdaf\n    " \
                                       "        5. py__g2l_obs_pdaf\n    "\
                                       "           (localise each ensemble member in observation space)\n    "\
                                       "        6. py__init_obsvar_l_pdaf\n    " \
                                       "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                       "        7. py__prodRinvA_l_pdaf\n    "\
                                       "        8. core DA algorithm\n" \
                                       "\n    " \
                                       ".. deprecated:: 1.0.0\n\n    " \
                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                       "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                       "       A unification of ensemble square root Kalman filters. \n    " \
                                       "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['local_put_state_letkf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                      "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                      "PDAFlocal-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                      "Local ensemble transform Kalman filter (LETKF) [1]_\n    " \
                                      "for a single DA step without OMI. Implementation is based on [2]_.\n\n    " \
                                      "Compared to :func:`pyPDAF.PDAF.local_assimilate_letkf`, this function has no :func:`get_state` call.\n    " \
                                      "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                      "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                      "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                      "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                      "function call to ensure the sequential DA.\n\n    " \
                                      "Note that the LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                      "This function should be called at each model time step.\n\n    " \
                                      "User-supplied functions are executed in the following sequence:\n    " \
                                      "    1. py__collect_state_pdaf\n    "\
                                      "    2. py__prepoststep_state_pdaf\n    "\
                                      "    3. py__init_n_domains_p_pdaf\n    "\
                                      "    4. py__init_dim_obs_pdaf\n    "\
                                      "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                      "    6. py__init_obs_pdaf\n    " \
                                      "       (if global adaptive forgetting factor `type_forget=1` is used\n    " \
                                      "       in :func:`pyPDAF.PDAF.init`)\n    "\
                                      "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    "\
                                      "    8. loop over each local domain:\n    " \
                                      "        1. py__init_dim_l_pdaf\n    "\
                                      "        2. py__init_dim_obs_l_pdaf\n    "\
                                      "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    "\
                                      "        4. py__init_obs_l_pdaf\n    " \
                                      "        5. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    "\
                                      "        6. py__init_obsvar_l_pdaf\n    " \
                                      "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    " \
                                      "        7. py__prodRinvA_l_pdaf\n    "\
                                      "        8. core DA algorithm\n" \
                                      "\n    " \
                                      ".. deprecated:: 1.0.0\n\n    " \
                                      "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                      "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                      "\n\n    " \
                                      "References\n    " \
                                      "----------\n    " \
                                      ".. [1] Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007).\n    "\
                                      "       Efficient data assimilation for spatiotemporal chaos:\n    "\
                                      "       A local ensemble transform Kalman filter. \n    "\
                                      "       Physica D: Nonlinear Phenomena, 230(1-2), 112-126.\n    " \
                                      ".. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                      "       A unification of ensemble square root Kalman filters. \n    " \
                                      "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['local_put_state_lseik'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                      "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
                                      "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                      "Local singular evolutive interpolated Kalman filter [1]_ for a single DA step.\n\n    " \
                                      "Compared to :func:`pyPDAF.PDAF.local_assimilate_lseik`, this function has no :func:`get_state` call.\n    " \
                                      "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                      "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                      "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                      "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                      "function call to ensure the sequential DA.\n\n    " \
                                      "This function should be called at each model time step.\n\n    " \
                                      "This function  executes the user-supplied functions in the following sequence:\n    " \
                                      "    1. py__collect_state_pdaf\n    " \
                                      "    2. py__prepoststep_state_pdaf\n    " \
                                      "    3. py__init_n_domains_p_pdaf\n    " \
                                      "    4. py__init_dim_obs_pdaf\n    " \
                                      "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                      "    6. py__init_obs_pdaf\n    "\
                                      "       (if global adaptive forgetting factor `type_forget=1`\n    " \
                                      "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                      "    7. py__init_obsvar_pdaf\n    "\
                                      "       (if global adaptive forgetting factor is used)\n    " \
                                      "    8. loop over each local domain:\n    " \
                                      "        1. py__init_dim_l_pdaf\n    " \
                                      "        2. py__init_dim_obs_l_pdaf\n    " \
                                      "        3. py__g2l_obs_pdaf (localise mean ensemble in observation space)\n    " \
                                      "        4. py__init_obs_l_pdaf\n    "\
                                      "        5. py__g2l_obs_pdaf\n    "\
                                      "           (localise each ensemble member in observation space)\n    " \
                                      "        6. py__init_obsvar_l_pdaf\n    "\
                                      "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                      "        7. py__prodRinvA_l_pdaf\n    " \
                                      "        8. core DA algorithm\n    " \
                                      "\n    " \
                                      ".. deprecated:: 1.0.0\n\n    " \
                                      "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                      "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
                                      "\n\n    " \
                                      "References\n    " \
                                      "----------\n    " \
                                      ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
                                      "       A singular evolutive extended Kalman filter for data assimilation\n    "\
                                      "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['local_put_state_lnetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                      "or :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`.\n\n    "\
                                      "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                      "Local Nonlinear Ensemble Transform Filter (LNETF) [1]_ for a single DA step.\n\n    " \
                                      "Compared to :func:`pyPDAF.PDAF.local_assimilate_lnetf`, this function has no :func:`get_state` call.\n    " \
                                      "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                      "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                      "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                      "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                      "function call to ensure the sequential DA.\n\n    " \
                                      "The nonlinear filter computes the distribution up to\n    " \
                                      "the second moment similar to Kalman filters but it uses a nonlinear weighting similar to\n    " \
                                      "particle filters. This leads to an equal weights assumption for the prior ensemble at each step.\n    " \
                                      "This function should be called at each model time step.\n\n    " \
                                      "This function executes the user-supplied function in the following sequence:\n    " \
                                      "    1. py__collect_state_pdaf\n    " \
                                      "    2. py__prepoststep_state_pdaf\n    " \
                                      "    3. py__init_n_domains_p_pdaf\n    " \
                                      "    4. py__init_dim_obs_pdaf\n    " \
                                      "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
                                      "    6. loop over each local domain:\n    " \
                                      "        1. py__init_dim_l_pdaf\n    " \
                                      "        2. py__init_dim_obs_l_pdaf\n    " \
                                      "        3. py__init_obs_l_pdaf\n    "\
                                      "        4. py__g2l_obs_pdaf (localise each ensemble member in observation space)\n    " \
                                      "        5. py__likelihood_l_pdaf\n    " \
                                      "        6. core DA algorithm\n" \
                                      "\n    " \
                                      ".. deprecated:: 1.0.0\n\n    " \
                                      "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                      "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`" \
                                      "\n\n    " \
                                      "References\n    " \
                                      "----------\n    " \
                                      ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                      "       A second-order exact ensemble square root filter\n    " \
                                      "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                      "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['local_put_state_lknetf'] = "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
                                       "or :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`.\n\n    "\
                                       "PDAF-OMI modules require fewer user-supplied functions and improved efficiency.\n\n    " \
                                       "A hybridised LETKF and LNETF [1]_ for a single DA step.\n\n    " \
                                       "Compared to :func:`pyPDAF.PDAF.local_assimilate_lknetf`, this function has no :func:`get_state` call.\n    " \
                                       "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                       "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                       "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                       "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                       "function call to ensure the sequential DA.\n\n    " \
                                       "The LNETF computes the distribution up to\n    " \
                                       "the second moment similar to Kalman filters but using a nonlinear weighting similar to\n    " \
                                       "particle filters. This leads to an equal weights assumption for the prior ensemble.\n    " \
                                       "The hybridisation with LETKF is expected to lead to improved performance for\n    " \
                                       "quasi-Gaussian problems.\n    " \
                                       "The function should be called at each model step.\n\n    " \
                                       "This function executes the user-supplied function in the following sequence:\n    " \
                                       "    1. py__collect_state_pdaf\n    " \
                                       "    2. py__prepoststep_state_pdaf\n    " \
                                       "    3. py__init_n_domains_p_pdaf\n    " \
                                       "    4. py__init_dim_obs_pdaf\n    " \
                                       "    5. py__obs_op_pdaf\n    "\
                                       "       (for each ensemble member)\n    " \
                                       "    6. py__init_obs_pdaf\n    " \
                                       "       (if global adaptive forgetting factor `type_forget=1`\n    "\
                                       "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
                                       "    7. py__init_obsvar_pdaf (if global adaptive forgetting factor is used)\n    " \
                                       "    8. loop over each local domain:\n    " \
                                       "        1. py__init_dim_l_pdaf\n    " \
                                       "        2. py__init_dim_obs_l_pdaf\n    " \
                                       "        3. py__g2l_obs_pdaf\n    "\
                                       "           (localise each ensemble member in observation space)\n    " \
                                       "        4. py__init_obs_l_pdaf\n    "\
                                       "        5. py__init_obsvar_l_pdaf\n    "\
                                       "           (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                       "        6. py__prodRinvA_pdaf\n    " \
                                       "        7. py__likelihood_l_pdaf\n    " \
                                       "        8. core DA algorithm\n    " \
                                       "    9. py__obs_op_pdaf\n    " \
                                       "       (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                       "    10. py__likelihood_hyb_l_pda\n    " \
                                       "    11. py__init_obsvar_l_pdaf\n    " \
                                       "        (only called if local adaptive forgetting factor `type_forget=2` is used)\n    "\
                                       "    12. py__prodRinvA_hyb_l_pdaf\n" \
                                       "\n    " \
                                       ".. deprecated:: 1.0.0\n\n    " \
                                       "   This function is replaced by :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
                                       "   and :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`" \
                                       "\n\n    " \
                                       "References\n    " \
                                       "----------\n    " \
                                       ".. [1] Nerger, L.. (2022) \n" \
                                       "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n" \
                                       "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

docstrings['localomi_assimilate'] = "Domain local filters for a single DA step using diagnoal observation error covariance matrix.\n\n    " \
                                    "Here, this function call is used for LE(S)TKF [1]_, \n    " \
                                    "LSEIK [1]_, LNETF [2]_, and LKNETF [3]_.\n    " \
                                    "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                    "This function should be called at each model time step.\n    " \
                                    "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_local`\n    " \
                                    "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                    "User-supplied functions are executed in the following sequence:\n    " \
                                    "    1. py__collect_state_pdaf\n    "\
                                    "    2. py__prepoststep_state_pdaf\n    "\
                                    "    3. py__init_n_domains_p_pdaf\n    "\
                                    "    4. py__init_dim_obs_pdaf\n    "\
                                    "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                    "    6. loop over each local domain:\n    " \
                                    "        1. py__init_dim_l_pdaf\n    "\
                                    "        2. py__init_dim_obs_l_pdaf\n    "\
                                    "        3. core DA algorithm\n    " \
                                    "    7. py__prepoststep_state_pdaf\n    "\
                                    "    8. py__distribute_state_pdaf\n    "\
                                    "    9. py__next_observation_pdaf\n" \
                                    "\n    " \
                                    "References\n    " \
                                    "----------\n    " \
                                    ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                    "       A unification of ensemble square root Kalman filters. \n    " \
                                    "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                    ".. [2] Tödter, J., and B. Ahrens, 2015:\n    "\
                                    "       A second-order exact ensemble square root filter\n    " \
                                    "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                    "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                    ".. [3] Nerger, L.. (2022) \n    " \
                                    "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                    "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['localomi_assimilate_en3dvar_lestkf'] = "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using diagnoal observation error covariance matrix.\n\n    " \
                                                   "The background error covariance matrix is estimated by ensemble.\n    " \
                                                   "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                   "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                   "This function should be called at each model time step.\n\n    " \
                                                   "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
                                                   "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                   "The user-supplied function are executed in the following sequence:\n    " \
                                                   "    1. py__collect_state_pdaf\n    " \
                                                   "    2. py__prepoststep_state_pdaf\n    " \
                                                   "    3. py__init_dim_obs_pdaf\n    " \
                                                   "    4. py__obs_op_pdaf\n    " \
                                                   "    5. Starting the iterative optimisation:\n    " \
                                                   "        1. py__cvt_ens_pdaf\n    " \
                                                   "        2. py__obs_op_lin_pdaf\n    " \
                                                   "        3. py__obs_op_adj_pdaf\n    " \
                                                   "        4. py__cvt_adj_ens_pdaf\n    " \
                                                   "        5. core DA algorithm\n    " \
                                                   "    6. py__cvt_ens_pdaf\n    " \
                                                   "    7. Perform LESTKF:\n    " \
                                                   "        1. py__init_n_domains_p_pdaf\n    " \
                                                   "        2. py__init_dim_obs_pdaf\n    " \
                                                   "        3. py__obs_op_pdaf\n    "\
                                                   "           (for each ensemble member)\n    " \
                                                   "        4. loop over each local domain:\n    " \
                                                   "            1. py__init_dim_l_pdaf\n    " \
                                                   "            2. py__init_dim_obs_l_pdaf\n    " \
                                                   "            3. core DA algorithm\n    " \
                                                   "    8. py__prepoststep_state_pdaf\n    " \
                                                   "    9. py__distribute_state_pdaf\n    " \
                                                   "    10. py__next_observation_pdaf"
docstrings['localomi_assimilate_hyb3dvar_lestkf'] = "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix.\n\n    " \
                                                    "Here, the background error covariance is hybridised by a static background error covariance,\n    " \
                                                    "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                                    "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                    "LESTKF in this implementation.\n    " \
                                                    "This function should be called at each model time step.\n\n    " \
                                                    "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
                                                    "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                    "The user-supplied functions are executed in the following sequence:\n    " \
                                                    "    1. py__collect_state_pdaf\n    " \
                                                    "    2. py__prepoststep_state_pdaf\n    " \
                                                    "    3. py__init_dim_obs_pdaf\n    " \
                                                    "    4. py__obs_op_pdaf\n    " \
                                                    "    5. The iterative optimisation:\n    " \
                                                    "        1. py__cvt_pdaf\n    " \
                                                    "        2. py__cvt_ens_pdaf\n    " \
                                                    "        3. py__obs_op_lin_pdaf\n    " \
                                                    "        4. py__obs_op_adj_pdaf\n    " \
                                                    "        5. py__cvt_adj_pdaf\n    " \
                                                    "        6. py__cvt_adj_ens_pdaf\n    " \
                                                    "        7. core DA algorithm\n    " \
                                                    "    6. py__cvt_pdaf\n    " \
                                                    "    7. py__cvt_ens_pdaf\n    " \
                                                    "    8. Perform LESTKF:\n    " \
                                                    "        1. py__init_n_domains_p_pdaf\n    " \
                                                    "        2. py__init_dim_obs_pdaf\n    " \
                                                    "        3. py__obs_op_pdaf\n    " \
                                                    "           (for each ensemble member)\n    " \
                                                    "        4. loop over each local domain:\n    " \
                                                    "            1. py__init_dim_l_pdaf\n    " \
                                                    "            2. py__init_dim_obs_l_pdaf\n    " \
                                                    "            3. core DA algorithm\n    " \
                                                    "    9. py__prepoststep_state_pdaf\n    " \
                                                    "    10. py__distribute_state_pdaf\n    " \
                                                    "    11. py__next_observation_pdaf"
docstrings['localomi_assimilate_en3dvar_lestkf_nondiagR'] = "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using non-diagnoal observation error covariance matrix.\n\n    " \
                                                            "Here, the background error covariance matrix is estimated by ensemble.\n    " \
                                                            "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                            "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                            "This function should be called at each model time step.\n\n    " \
                                                            "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`\n    " \
                                                            "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                            "The user-supplied function are executed in the following sequence:\n    " \
                                                            "    1. py__collect_state_pdaf\n    " \
                                                            "    2. py__prepoststep_state_pdaf\n    " \
                                                            "    3. py__init_dim_obs_pdaf\n    " \
                                                            "    4. py__obs_op_pdaf\n    " \
                                                            "    5. Starting the iterative optimisation:\n    " \
                                                            "        1. py__cvt_ens_pdaf\n    " \
                                                            "        2. py__obs_op_lin_pdaf\n    " \
                                                            "        3. py__prodRinvA_pdaf\n    " \
                                                            "        4. py__obs_op_adj_pdaf\n    " \
                                                            "        5. py__cvt_adj_ens_pdaf\n    " \
                                                            "        6. core DA algorithm\n    " \
                                                            "    6. py__cvt_ens_pdaf\n    " \
                                                            "    7. Perform LESTKF:\n    " \
                                                            "        1. py__init_n_domains_p_pdaf\n    " \
                                                            "        2. py__init_dim_obs_pdaf\n    " \
                                                            "        3. py__obs_op_pdaf\n    "\
                                                            "           (for each ensemble member)\n    " \
                                                            "        4. loop over each local domain:\n    " \
                                                            "            1. py__init_dim_l_pdaf\n    " \
                                                            "            2. py__init_dim_obs_l_pdaf\n    " \
                                                            "            3. py__prodRinvA_l_pdaf\n    " \
                                                            "            4. core DA algorithm\n    " \
                                                            "    8. py__prepoststep_state_pdaf\n    " \
                                                            "    9. py__distribute_state_pdaf\n    " \
                                                            "    10. py__next_observation_pdaf"
docstrings['localomi_assimilate_hyb3dvar_lestkf_nondiagR'] = "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix.\n\n    " \
                                                             "Here, the background error covariance is hybridised by a static background error covariance,\n    " \
                                                             "and a flow-dependent background error covariance estimated from ensemble.\n    " \
                                                             "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                             "LESTKF in this implementation.\n    " \
                                                             "This function should be called at each model time step.\n\n    " \
                                                             "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`\n    " \
                                                             "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                             "The user-supplied functions are executed in the following sequence:\n    " \
                                                             "    1. py__collect_state_pdaf\n    " \
                                                             "    2. py__prepoststep_state_pdaf\n    " \
                                                             "    3. py__init_dim_obs_pdaf\n    " \
                                                             "    4. py__obs_op_pdaf\n    " \
                                                             "    5. The iterative optimisation:\n    " \
                                                             "        1. py__cvt_pdaf\n    " \
                                                             "        2. py__cvt_ens_pdaf\n    " \
                                                             "        3. py__obs_op_lin_pdaf\n    " \
                                                             "        4. py__prodRinvA_pdaf\n    " \
                                                             "        5. py__obs_op_adj_pdaf\n    " \
                                                             "        6. py__cvt_adj_pdaf\n    " \
                                                             "        7. py__cvt_adj_ens_pdaf\n    " \
                                                             "        8. core DA algorithm\n    " \
                                                             "    6. py__cvt_pdaf\n    " \
                                                             "    7. py__cvt_ens_pdaf\n    " \
                                                             "    8. Perform LESTKF:\n    " \
                                                             "        1. py__init_n_domains_p_pdaf\n    " \
                                                             "        2. py__init_dim_obs_pdaf\n    " \
                                                             "        3. py__obs_op_pdaf\n    " \
                                                             "           (for each ensemble member)\n    " \
                                                             "        4. loop over each local domain:\n    " \
                                                             "            1. py__init_dim_l_pdaf\n    " \
                                                             "            2. py__init_dim_obs_l_pdaf\n    " \
                                                             "            3. py__prodRinvA_l_pdaf\n    " \
                                                             "            4. core DA algorithm\n    " \
                                                             "    9. py__prepoststep_state_pdaf\n    " \
                                                             "    10. py__distribute_state_pdaf\n    " \
                                                             "    11. py__next_observation_pdaf"
docstrings['localomi_assimilate_nondiagR'] = "Domain local filters for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                             "Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_\n    " \
                                             "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                             "This function should be called at each model time step.\n    " \
                                             "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_local_nondiagR`\n    " \
                                             "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                             "User-supplied functions are executed in the following sequence:\n    " \
                                             "    1. py__collect_state_pdaf\n    "\
                                             "    2. py__prepoststep_state_pdaf\n    "\
                                             "    3. py__init_n_domains_p_pdaf\n    "\
                                             "    4. py__init_dim_obs_pdaf\n    "\
                                             "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                             "    6. loop over each local domain:\n    " \
                                             "        1. py__init_dim_l_pdaf\n    "\
                                             "        2. py__init_dim_obs_l_pdaf\n    "\
                                             "        3. py__init_obs_l_pdaf\n    "\
                                             "        4. py__prodRinvA_l_pdaf\n    " \
                                             "        5. core DA algorithm\n    " \
                                             "    7. py__prepoststep_state_pdaf\n    "\
                                             "    8. py__distribute_state_pdaf\n    "\
                                             "    9. py__next_observation_pdaf\n" \
                                             "\n    " \
                                             "References\n    " \
                                             "----------\n    " \
                                             ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                             "       A unification of ensemble square root Kalman filters. \n    " \
                                             "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['localomi_assimilate_lnetf_nondiagR'] = "LNETF for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                                   "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                                   "The non-linear filter is proposed in [1]_.\n    " \
                                                   "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                                   "This function should be called at each model time step.\n    " \
                                                   "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`\n    " \
                                                   "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                   "User-supplied functions are executed in the following sequence:\n    " \
                                                   "    1. py__collect_state_pdaf\n    "\
                                                   "    2. py__prepoststep_state_pdaf\n    "\
                                                   "    3. py__init_n_domains_p_pdaf\n    "\
                                                   "    4. py__init_dim_obs_pdaf\n    "\
                                                   "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                                   "    6. loop over each local domain:\n    " \
                                                   "        1. py__init_dim_l_pdaf\n    "\
                                                   "        2. py__init_dim_obs_l_pdaf\n    "\
                                                   "        3. py__likelihood_l_pdaf\n    " \
                                                   "        4. core DA algorithm\n    " \
                                                   "    7. py__prepoststep_state_pdaf\n    "\
                                                   "    8. py__distribute_state_pdaf\n    "\
                                                   "    9. py__next_observation_pdaf\n" \
                                                   "\n    " \
                                                   "References\n    " \
                                                   "----------\n    " \
                                                   ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                                   "       A second-order exact ensemble square root filter\n    " \
                                                   "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                                   "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['localomi_assimilate_lknetf_nondiagR'] = "LKNETF for a single DA step using non-diagnoal observation error covariance matrix.\n\n    " \
                                                    "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                                    "The non-linear filter is proposed in [1]_.\n    " \
                                                    "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                                    "This function should be called at each model time step.\n    " \
                                                    "The function is a combination of :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`\n    " \
                                                    "and :func:`pyPDAF.PDAF.get_state`.\n\n    " \
                                                    "User-supplied functions are executed in the following sequence:\n    " \
                                                    "    1. py__collect_state_pdaf\n    "\
                                                    "    2. py__prepoststep_state_pdaf\n    "\
                                                    "    3. py__init_n_domains_p_pdaf\n    "\
                                                    "    4. py__init_dim_obs_pdaf\n    "\
                                                    "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                                    "    6. loop over each local domain:\n    " \
                                                    "        1. py__init_dim_l_pdaf\n    "\
                                                    "        2. py__init_dim_obs_l_pdaf\n    "\
                                                    "        3. py__prodRinvA_pdaf\n    " \
                                                    "        4. py__likelihood_l_pdaf\n    " \
                                                    "        5. core DA algorithm\n    " \
                                                    "        6. py__obs_op_pdaf\n    " \
                                                    "           (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                                    "        7. py__likelihood_hyb_l_pda\n    " \
                                                    "        8. py__prodRinvA_hyb_l_pdaf\n    " \
                                                    "    7. py__prepoststep_state_pdaf\n    "\
                                                    "    8. py__distribute_state_pdaf\n    "\
                                                    "    9. py__next_observation_pdaf\n" \
                                                    "\n    " \
                                                    "References\n    " \
                                                    "----------\n    " \
                                                    ".. [1] Nerger, L.. (2022) \n    " \
                                                    "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                                    "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

docstrings['localomi_put_state'] = "Domain local filters for a single DA step using diagnoal observation error covariance matrix\n    " \
                                   "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                   "Here, this function call is used for LE(S)TKF [1]_, \n    " \
                                   "LSEIK [1]_, LNETF [2]_, and LKNETF [3]_.\n    " \
                                   "The filter type is set in :func:`pyPDAF.PDAF.init`.\n    " \
                                   "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_local`, this function has no :func:`get_state` call.\n    " \
                                   "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                   "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                   "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                   "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                   "function call to ensure the sequential DA.\n\n    " \
                                   "The LESTKF is a more efficient equivalent to the LETKF.\n\n    " \
                                   "This function should be called at each model time step.\n\n    " \
                                   "User-supplied functions are executed in the following sequence:\n    " \
                                   "    1. py__collect_state_pdaf\n    "\
                                   "    2. py__prepoststep_state_pdaf\n    "\
                                   "    3. py__init_n_domains_p_pdaf\n    "\
                                   "    4. py__init_dim_obs_pdaf\n    "\
                                   "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                   "    6. loop over each local domain:\n    " \
                                   "        1. py__init_dim_l_pdaf\n    "\
                                   "        2. py__init_dim_obs_l_pdaf\n    "\
                                   "        3. core DA algorithm\n" \
                                   "\n    " \
                                   "References\n    " \
                                   "----------\n    " \
                                   ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                   "       A unification of ensemble square root Kalman filters. \n    " \
                                   "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1\n    " \
                                   ".. [2] Tödter, J., and B. Ahrens, 2015:\n    "\
                                   "       A second-order exact ensemble square root filter\n    " \
                                   "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                   "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.\n    " \
                                   ".. [3] Nerger, L.. (2022) \n    " \
                                   "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter. \n    " \
                                   "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['localomi_put_state_en3dvar_lestkf'] = "3DEnVar for a single DA step where the ensemble anomaly is generated by LESTKF using diagnoal observation error covariance matrix\n    " \
                                                  "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                  "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                                  "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                  "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                  "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                  "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                  "function call to ensure the sequential DA.\n\n    " \
                                                  "The background error covariance matrix is estimated by ensemble.\n    " \
                                                  "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                  "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                  "This function should be called at each model time step.\n\n    " \
                                                  "The user-supplied function are executed in the following sequence:\n    " \
                                                  "    1. py__collect_state_pdaf\n    " \
                                                  "    2. py__prepoststep_state_pdaf\n    " \
                                                  "    3. py__init_dim_obs_pdaf\n    " \
                                                  "    4. py__obs_op_pdaf\n    " \
                                                  "    5. Starting the iterative optimisation:\n    " \
                                                  "        1. py__cvt_ens_pdaf\n    " \
                                                  "        2. py__obs_op_lin_pdaf\n    " \
                                                  "        3. py__obs_op_adj_pdaf\n    " \
                                                  "        4. py__cvt_adj_ens_pdaf\n    " \
                                                  "        5. core DA algorithm\n    " \
                                                  "    6. py__cvt_ens_pdaf\n    " \
                                                  "    7. Perform LESTKF:\n    " \
                                                  "        1. py__init_n_domains_p_pdaf\n    " \
                                                  "        2. py__init_dim_obs_pdaf\n    " \
                                                  "        3. py__obs_op_pdaf\n    "\
                                                  "           (for each ensemble member)\n    " \
                                                  "        4. loop over each local domain:\n    " \
                                                  "            1. py__init_dim_l_pdaf\n    " \
                                                  "            2. py__init_dim_obs_l_pdaf\n    " \
                                                  "            3. core DA algorithm"
docstrings['localomi_put_state_hyb3dvar_lestkf'] = "Hybrid 3DEnVar for a single DA step using diagnoal observation error covariance matrix\n    " \
                                                   "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                   "Here, the background error covariance is hybridised by a static background error covariance,\n    " \
                                                   "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                                   "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`, this function has no :func:`get_state` call.\n    " \
                                                   "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                   "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                   "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                   "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                   "function call to ensure the sequential DA.\n\n    " \
                                                   "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                   "LESTKF in this implementation.\n    " \
                                                   "This function should be called at each model time step.\n\n    " \
                                                   "The user-supplied functions are executed in the following sequence:\n    " \
                                                   "    1. py__collect_state_pdaf\n    " \
                                                   "    2. py__prepoststep_state_pdaf\n    " \
                                                   "    3. py__init_dim_obs_pdaf\n    " \
                                                   "    4. py__obs_op_pdaf\n    " \
                                                   "    5. The iterative optimisation:\n    " \
                                                   "        1. py__cvt_pdaf\n    " \
                                                   "        2. py__cvt_ens_pdaf\n    " \
                                                   "        3. py__obs_op_lin_pdaf\n    " \
                                                   "        4. py__obs_op_adj_pdaf\n    " \
                                                   "        5. py__cvt_adj_pdaf\n    " \
                                                   "        6. py__cvt_adj_ens_pdaf\n    " \
                                                   "        7. core DA algorithm\n    " \
                                                   "    6. py__cvt_pdaf\n    " \
                                                   "    7. py__cvt_ens_pdaf\n    " \
                                                   "    8. Perform LESTKF:\n    " \
                                                   "        1. py__init_n_domains_p_pdaf\n    " \
                                                   "        2. py__init_dim_obs_pdaf\n    " \
                                                   "        3. py__obs_op_pdaf\n    " \
                                                   "           (for each ensemble member)\n    " \
                                                   "        4. loop over each local domain:\n    " \
                                                   "            1. py__init_dim_l_pdaf\n    " \
                                                   "            2. py__init_dim_obs_l_pdaf\n    " \
                                                   "            3. core DA algorithm"
docstrings['localomi_put_state_en3dvar_lestkf_nondiagR'] = "3DEnVar for a single DA step without post-processing, distributing analysis, and setting next observation step.\n\n     "\
                                                           "Here, the ensemble anomaly is generated by LESTKF using non-diagnoal observation error covariance matrix.\n\n    " \
                                                           "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                           "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                           "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                           "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                           "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                           "function call to ensure the sequential DA.\n\n    " \
                                                           "The background error covariance matrix is estimated by ensemble.\n    " \
                                                           "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
                                                           "An LESTKF is used to generate ensemble perturbations.\n    " \
                                                           "This function should be called at each model time step.\n\n    " \
                                                           "The user-supplied function are executed in the following sequence:\n    " \
                                                           "    1. py__collect_state_pdaf\n    " \
                                                           "    2. py__prepoststep_state_pdaf\n    " \
                                                           "    3. py__init_dim_obs_pdaf\n    " \
                                                           "    4. py__obs_op_pdaf\n    " \
                                                           "    5. Starting the iterative optimisation:\n    " \
                                                           "        1. py__cvt_ens_pdaf\n    " \
                                                           "        2. py__obs_op_lin_pdaf\n    " \
                                                           "        3. py__prodRinvA_pdaf\n    " \
                                                           "        4. py__obs_op_adj_pdaf\n    " \
                                                           "        5. py__cvt_adj_ens_pdaf\n    " \
                                                           "        6. core DA algorithm\n    " \
                                                           "    6. py__cvt_ens_pdaf\n    " \
                                                           "    7. Perform LESTKF:\n    " \
                                                           "        1. py__init_n_domains_p_pdaf\n    " \
                                                           "        2. py__init_dim_obs_pdaf\n    " \
                                                           "        3. py__obs_op_pdaf\n    "\
                                                           "           (for each ensemble member)\n    " \
                                                           "        4. loop over each local domain:\n    " \
                                                           "            1. py__init_dim_l_pdaf\n    " \
                                                           "            2. py__init_dim_obs_l_pdaf\n    " \
                                                           "            3. py__prodRinvA_l_pdaf\n    " \
                                                           "            4. core DA algorithm"
docstrings['localomi_put_state_hyb3dvar_lestkf_nondiagR'] = "Hybrid 3DEnVar for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                            "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                            "Here, the background error covariance is hybridised by a static background error covariance,\n    " \
                                                            "and a flow-dependent background error covariance estimated from ensemble.\n\n    " \
                                                            "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                            "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                            "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                            "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                            "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                            "function call to ensure the sequential DA.\n\n    " \
                                                            "The 3DVar generates an ensemble mean and the ensemble perturbation is generated by\n    " \
                                                            "LESTKF in this implementation.\n    " \
                                                            "This function should be called at each model time step.\n\n    " \
                                                            "The user-supplied functions are executed in the following sequence:\n    " \
                                                            "    1. py__collect_state_pdaf\n    " \
                                                            "    2. py__prepoststep_state_pdaf\n    " \
                                                            "    3. py__init_dim_obs_pdaf\n    " \
                                                            "    4. py__obs_op_pdaf\n    " \
                                                            "    5. The iterative optimisation:\n    " \
                                                            "        1. py__cvt_pdaf\n    " \
                                                            "        2. py__cvt_ens_pdaf\n    " \
                                                            "        3. py__obs_op_lin_pdaf\n    " \
                                                            "        4. py__prodRinvA_pdaf\n    " \
                                                            "        5. py__obs_op_adj_pdaf\n    " \
                                                            "        6. py__cvt_adj_pdaf\n    " \
                                                            "        7. py__cvt_adj_ens_pdaf\n    " \
                                                            "        8. core DA algorithm\n    " \
                                                            "    6. py__cvt_pdaf\n    " \
                                                            "    7. py__cvt_ens_pdaf\n    " \
                                                            "    8. Perform LESTKF:\n    " \
                                                            "        1. py__init_n_domains_p_pdaf\n    " \
                                                            "        2. py__init_dim_obs_pdaf\n    " \
                                                            "        3. py__obs_op_pdaf\n    " \
                                                            "           (for each ensemble member)\n    " \
                                                            "        4. loop over each local domain:\n    " \
                                                            "            1. py__init_dim_l_pdaf\n    " \
                                                            "            2. py__init_dim_obs_l_pdaf\n    " \
                                                            "            3. py__prodRinvA_l_pdaf\n    " \
                                                            "            4. core DA algorithm"
docstrings['localomi_put_state_nondiagR'] = "Domain local filters for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                            "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                            "Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_\n    " \
                                            "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                            "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_local_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                            "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                            "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                            "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                            "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                            "function call to ensure the sequential DA.\n\n    " \
                                            "This function should be called at each model time step.\n\n    " \
                                            "User-supplied functions are executed in the following sequence:\n    " \
                                            "    1. py__collect_state_pdaf\n    "\
                                            "    2. py__prepoststep_state_pdaf\n    "\
                                            "    3. py__init_n_domains_p_pdaf\n    "\
                                            "    4. py__init_dim_obs_pdaf\n    "\
                                            "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                            "    6. loop over each local domain:\n    " \
                                            "        1. py__init_dim_l_pdaf\n    "\
                                            "        2. py__init_dim_obs_l_pdaf\n    "\
                                            "        3. py__init_obs_l_pdaf\n    "\
                                            "        4. py__prodRinvA_l_pdaf\n    " \
                                            "        5. core DA algorithm\n    " \
                                            "\n    " \
                                            "References\n    " \
                                            "----------\n    " \
                                            ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
                                            "       A unification of ensemble square root Kalman filters. \n    " \
                                            "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['localomi_put_state_lnetf_nondiagR'] = "LNETF for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                  "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                  "See :func:`pyPDAF.PDAF.localomi_put_state` for using diagnoal observation error covariance matrix.\n    " \
                                                  "The non-linear filter is proposed in [1]_.\n    " \
                                                  "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                                  "Compared to :func:`pyPDAF.PDAF.omi_assimilate_lnetf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                  "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                  "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                  "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                  "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                  "function call to ensure the sequential DA.\n\n    " \
                                                  "This function should be called at each model time step.\n\n    " \
                                                  "User-supplied functions are executed in the following sequence:\n    " \
                                                  "    1. py__collect_state_pdaf\n    "\
                                                  "    2. py__prepoststep_state_pdaf\n    "\
                                                  "    3. py__init_n_domains_p_pdaf\n    "\
                                                  "    4. py__init_dim_obs_pdaf\n    "\
                                                  "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                                  "    6. loop over each local domain:\n    " \
                                                  "        1. py__init_dim_l_pdaf\n    "\
                                                  "        2. py__init_dim_obs_l_pdaf\n    "\
                                                  "        3. py__likelihood_l_pdaf\n    " \
                                                  "        4. core DA algorithm\n    " \
                                                  "\n    " \
                                                  "References\n    " \
                                                  "----------\n    " \
                                                  ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
                                                  "       A second-order exact ensemble square root filter\n    " \
                                                  "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
                                                  "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['localomi_put_state_lknetf_nondiagR'] = "LKNETF for a single DA step using non-diagnoal observation error covariance matrix\n    " \
                                                   "without post-processing, distributing analysis, and setting next observation step.\n\n    " \
                                                   "See :func:`pyPDAF.PDAF.localomi_assimilate` for using diagnoal observation error covariance matrix.\n    " \
                                                   "The non-linear filter is proposed in [1]_.\n    " \
                                                   "The filter type is set in :func:`pyPDAF.PDAF.init`.\n\n    " \
                                                   "Compared to :func:`pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR`, this function has no :func:`get_state` call.\n    " \
                                                   "This means that the analysis is not post-processed, and distributed to the model forecast\n    " \
                                                   "by user-supplied functions. The next DA step will not be assigned by user-supplied functions as well.\n    " \
                                                   "This function is typically used when there are not enough CPUs to run the ensemble in parallel,\n    "\
                                                   "and some ensemble members have to be run serially. The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
                                                   "function call to ensure the sequential DA.\n\n    " \
                                                   "This function should be called at each model time step.\n\n    " \
                                                   "User-supplied functions are executed in the following sequence:\n    " \
                                                   "    1. py__collect_state_pdaf\n    "\
                                                   "    2. py__prepoststep_state_pdaf\n    "\
                                                   "    3. py__init_n_domains_p_pdaf\n    "\
                                                   "    4. py__init_dim_obs_pdaf\n    "\
                                                   "    5. py__obs_op_pdaf (for each ensemble member)\n    "\
                                                   "    6. loop over each local domain:\n    " \
                                                   "        1. py__init_dim_l_pdaf\n    "\
                                                   "        2. py__init_dim_obs_l_pdaf\n    "\
                                                   "        3. py__prodRinvA_pdaf\n    " \
                                                   "        4. py__likelihood_l_pdaf\n    " \
                                                   "        5. core DA algorithm\n    " \
                                                   "        6. py__obs_op_pdaf\n    " \
                                                   "           (only called with `HKN` and `HNK` options called for each ensemble member)\n    " \
                                                   "        7. py__likelihood_hyb_l_pda\n    " \
                                                   "        8. py__prodRinvA_hyb_l_pdaf\n" \
                                                   "\n    " \
                                                   "References\n    " \
                                                   "----------\n    " \
                                                   ".. [1] Nerger, L.. (2022) \n    " \
                                                   "       Data assimilation for nonlinear systems with a hybrid nonlinear Kalman ensemble transform filter.\n    " \
                                                   "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"

