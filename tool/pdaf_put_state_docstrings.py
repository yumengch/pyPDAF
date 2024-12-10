"""docstrings for PDAF_put_state_xxx functions

These functions are mostly deprecated.
"""

docstrings = {}

docstrings['put_state_3dvar'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_3dvar`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied\n    " \
    "functions and improved efficiency.\n\n    " \
    "3DVar DA for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_3dvar`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will not\n    " \
    "be assigned by user-supplied functions as well.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "When 3DVar is used, the background error covariance matrix\n    "\
    "has to be modelled for cotrol variable transformation.\n    " \
    "This is a deterministic filtering scheme so no ensemble\n    " \
    "and parallelisation is needed.\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_3dvar`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_3dvar_nondiagR`"
docstrings['put_state_en3dvar_estkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`." \
    "\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "3DEnVar for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_estkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will not be\n    " \
    "assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The background error covariance matrix is\n    " \
    "estimated by an ensemble.\n    " \
    "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
    "An ESTKF is used along with 3DEnVar to\n    " \
    "generate ensemble perturbations.\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR`"
docstrings['put_state_en3dvar_lestkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "3DEnVar for a single DA step without post-processing,\n    " \
    "distributing analysis, and setting next observation step,\n    "\
    "where the ensemble anomaly is generated by LESTKF.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_en3dvar_lestkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The background error covariance matrix is estimated by ensemble.\n    " \
    "The 3DEnVar only calculates the analysis of the ensemble mean.\n    " \
    "An LESTKF is used to generate ensemble perturbations.\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR`"
docstrings['put_state_hyb3dvar_estkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Hybrid 3DEnVar for a single DA step where\n    " \
    "the background error covariance is hybridised by\n    " \
    "a static background error covariance,\n    " \
    "and a flow-dependent background error covariance\n    " \
    "estimated from ensemble.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_estkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The 3DVar generates an ensemble mean and\n    " \
    "the ensemble perturbation is generated by\n    " \
    "ESTKF in this implementation.\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "           (only relevant for adaptive\n    " \
    "           forgetting factor schemes)\n    " \
    "        6. py__prodRinvA_pdaf\n    " \
    "        7. core ESTKF algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR`"
docstrings['put_state_hyb3dvar_lestkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Hybrid 3DEnVar for a single DA step using\n    " \
    "non-diagnoal observation error covariance matrix\n    " \
    "without post-processing, distributing analysis,\n    " \
    "and setting next observation step, where\n    " \
    "the background error covariance is hybridised by\n    " \
    "a static background error covariance,\n    " \
    "and a flow-dependent background error covariance\n    " \
    "estimated from ensemble.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_hyb3dvar_lestkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The 3DVar generates an ensemble mean and\n    " \
    "the ensemble perturbation is generated by\n    " \
    "LESTKF in this implementation.\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "           (if global adaptive forgetting factor\n    " \
    "           `type_forget=1` in :func:`pyPDAF.PDAF.init`)\n    " \
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
    "               (localise each ensemble member\n    " \
    "               in observation space)\n    " \
    "            7. py__init_obsvar_l_pdaf\n    " \
    "               (only called if local adaptive forgetting\n    " \
    "               factor `type_forget=2` is used)\n    "\
    "            8. py__prodRinvA_l_pdaf\n    " \
    "            9. core DA algorithm\n    " \
    "            10. py__l2g_state_pdaf\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR`"
docstrings['put_state_enkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Stochastic EnKF (ensemble " \
    "Kalman filter) [1]_ for a single DA step without OMI.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_enkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This function should be called at each model time step. \n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_enkf_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Evensen, G. (1994), \n    "\
    "       Sequential data assimilation with a\n    " \
    "       nonlinear quasi-geostrophic model\n    "\
    "       using Monte Carlo methods to forecast error statistics,\n    "\
    "       J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572."
docstrings['put_state_estkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`\n    " \
    "instead of this function.\n\n    " \
    "OMI functions need fewer user-supplied functions\n    " \
    "and improve DA efficiency.\n\n    " \
    "This function calls ESTKF\n    " \
    "(error space transform Kalman filter) [1]_.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_estkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
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
    "    7. py__init_obsvar_pdaf (only relevant for\n    " \
    "       adaptive forgetting factor schemes)\n    " \
    "    8. py__prodRinvA_pdaf\n    " \
    "    9. core DA algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`." \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
    "       A unification of ensemble square root Kalman filters. \n    " \
    "       Monthly Weather Review, 140, 2335-2345. doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_etkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
    "PDAFlocal-OMI modules require fewer user-supplied\n    " \
    "functions and improved efficiency.\n\n    " \
    "Using ETKF (ensemble transform " \
    "Kalman filter) [1]_ for a single DA step without OMI.\n    " \
    "The implementation is baed on [2]_.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_etkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_dim_obs_pdaf\n    " \
    "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
    "    5. py__init_obs_pdaf\n    " \
    "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
    "    7. py__init_obsvar_pdaf (only relevant for\n    " \
    "       adaptive forgetting factor schemes)\n    " \
    "    8. py__prodRinvA_pdaf\n    " \
    "    9. core DA algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
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
docstrings['put_state_seek'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "This function will use\n    " \
    "singular evolutive extended Kalman filter [1]_ for\n    " \
    "a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_seek`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This is a deterministic Kalman filter.\n    " \
    "The function should be called at each model step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
    "       A singular evolutive extended Kalman filter for data assimilation\n    "\
    "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_seik'] =  \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "This function will use\n    " \
    "singular evolutive interpolated Kalman filter [1]_ for\n    " \
    "a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_seik`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions.\n    " \
    "The next DA step will not be assigned by user-supplied\n    " \
    "functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The function should be called at each model step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_dim_obs_pdaf\n    " \
    "    4. py__obs_op_pdaf (for ensemble mean)\n    " \
    "    5. py__init_obs_pdaf\n    " \
    "    6. py__obs_op_pdaf (for each ensemble member)\n    " \
    "    7. py__init_obsvar_pdaf (only relevant for\n    " \
    "       adaptive forgetting factor schemes)\n    " \
    "    8. py__prodRinvA_pdaf\n    " \
    "    9. core DA algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
    "       A singular evolutive extended Kalman filter\n    " \
    "       for data assimilation\n    "\
    "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_netf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "This function will use\n    " \
    "Nonlinear Ensemble Transform Filter (NETF) [1]_ \n    " \
    "for a single DA step.\n\n    "\
    "Compared to :func:`pyPDAF.PDAF.assimilate_netf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The nonlinear filter computes the distribution up to\n    " \
    "the second moment similar to KF but using\n    " \
    "a nonlinear weighting similar to\n    " \
    "particle filter. This leads to an equal weights\n    " \
    "assumption for prior ensemble.\n    " \
    "The function should be called at each model step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_dim_obs_pdaf\n    " \
    "    4. py__init_obs_pdaf\n    " \
    "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
    "    6. py__likelihood_pdaf\n    " \
    "    7. core DA algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
    "       A second-order exact ensemble square root filter\n    " \
    "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
    "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['put_state_pf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_global`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "This function will use particle filter for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_pf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This is a fully nonlinear filter, and may require\n    " \
    "a high number of ensemble members.\n    " \
    "A review of particle filter can be found at [1]_.\n    " \
    "The function should be called at each model step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_dim_obs_pdaf\n    " \
    "    4. py__init_obs_pdaf\n    " \
    "    5. py__obs_op_pdaf (for each ensemble member)\n    " \
    "    6. py__likelihood_pdaf\n    " \
    "    7. core DA algorithm\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_global`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_nonlin_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Van Leeuwen, P. J., Künsch, H. R., Nerger, L.,\n    " \
    "       Potthast, R., & Reich, S. (2019).\n    "\
    "       Particle filters for high‐dimensional\n    " \
    "       geoscience applications:\n    "\
    "       A review. \n    " \
    "       Quarterly Journal of the Royal Meteorological Society,\n    " \
    "       145(723), 2335-2365."
docstrings['put_state_lenkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_put_state_lenkf`\n    "\
    "or :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Stochastic EnKF (ensemble Kalman filter)\n    " \
    "with covariance localisation [1]_\n    " \
    "for a single DA step without OMI.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_lenkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This is the only scheme for covariance localisation in PDAF." \
    "\n\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.omi_put_state_lenkf`\n    " \
    "   and :func:`pyPDAF.PDAF.omi_put_state_lenkf_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Houtekamer, P. L., and H. L. Mitchell (1998): \n    " \
    "       Data Assimilation Using an Ensemble Kalman\n    " \
    "       Filter Technique.\n    "\
    "       Mon. Wea. Rev., 126, 796–811,\n    "\
    "       doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2."
docstrings['put_state_lestkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
    "PDAFlocal-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Local ESTKF (error space transform " \
    "Kalman filter) [1]_ for a single DA step without OMI.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_lestkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
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
    "    7. py__init_obsvar_pdaf (if global adaptive\n    " \
    "       forgetting factor is used)\n    "\
    "    8. loop over each local domain:\n    " \
    "        1. py__init_dim_l_pdaf\n    "\
    "        2. py__init_dim_obs_l_pdaf\n    "\
    "        3. py__g2l_state_pdaf\n    "\
    "        4. py__g2l_obs_pdaf (localise mean ensemble\n    " \
    "           in observation space)\n    "\
    "        5. py__init_obs_l_pdaf\n    " \
    "        6. py__g2l_obs_pdaf\n    "\
    "           (localise each ensemble member in observation space)\n    "\
    "        7. py__init_obsvar_l_pdaf\n    " \
    "           (only called if local adaptive forgetting factor\n    " \
    "           `type_forget=2` is used)\n    " \
    "        8. py__prodRinvA_l_pdaf\n    "\
    "        9. core DA algorithm\n    " \
    "        10. py__l2g_state_pdaf\n"\
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012). \n    " \
    "       A unification of ensemble square root Kalman filters. \n    " \
    "       Monthly Weather Review, 140, 2335-2345.\n    " \
    "       doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_letkf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
    "PDAFlocal-OMI modules require fewer user-supplied\n    " \
    "functions and improved efficiency.\n\n    " \
    "Local ensemble transform Kalman filter (LETKF) [1]_\n    " \
    "for a single DA step without OMI.\n    " \
    "Implementation is based on [2]_.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_letkf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
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
    "       (if global adaptive forgetting factor\n    " \
    "       `type_forget=1` is used\n    " \
    "       in :func:`pyPDAF.PDAF.init`)\n    "\
    "    7. py__init_obsvar_pdaf (if global adaptive forgetting\n    " \
    "       factor is used)\n    "\
    "    8. loop over each local domain:\n    " \
    "        1. py__init_dim_l_pdaf\n    "\
    "        2. py__init_dim_obs_l_pdaf\n    "\
    "        3. py__g2l_state_pdaf\n    "\
    "        4. py__g2l_obs_pdaf (localise mean ensemble\n    " \
    "           in observation space)\n    "\
    "        5. py__init_obs_l_pdaf\n    " \
    "        6. py__g2l_obs_pdaf (localise each ensemble member\n    " \
    "           in observation space)\n    "\
    "        7. py__init_obsvar_l_pdaf\n    " \
    "           (only called if local adaptive forgetting factor\n    " \
    "           `type_forget=2` is used)\n    " \
    "        8. py__prodRinvA_l_pdaf\n    "\
    "        9. core DA algorithm\n    " \
    "        10. py__l2g_state_pdaf\n"\
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
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
    "       Monthly Weather Review, 140, 2335-2345.\n    " \
    "       doi:10.1175/MWR-D-11-00102.1"
docstrings['put_state_lseik'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Local singular evolutive interpolated Kalman filter [1]_\n    " \
    "for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_lseik`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "        4. py__g2l_obs_pdaf (localise mean ensemble\n    " \
    "           in observation space)\n    " \
    "        5. py__init_obs_l_pdaf\n    "\
    "        6. py__g2l_obs_pdaf\n    "\
    "           (localise each ensemble member in observation space)\n    " \
    "        7. py__init_obsvar_l_pdaf\n    "\
    "           (only called if local adaptive forgetting factor\n    " \
    "           `type_forget=2` is used)\n    "\
    "        8. py__prodRinvA_l_pdaf\n    " \
    "        9. core DA algorithm\n    " \
    "        10. py__l2g_state_pdaf\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Pham, D. T., Verron, J., & Roubaud, M. C. (1998).\n    "\
    "       A singular evolutive extended Kalman filter for data assimilation\n    "\
    "       in oceanography. Journal of Marine systems, 16(3-4), 323-340."
docstrings['put_state_lnetf'] = \
    "It is recommended to use :func:`pyPDAF.PDAF.localomi_put_state`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "Local Nonlinear Ensemble Transform Filter (LNETF) [1]_\n    " \
    "for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_lnetf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The nonlinear filter computes the distribution up to\n    " \
    "the second moment similar to Kalman filters\n    " \
    "but it uses a nonlinear weighting similar to\n    " \
    "particle filters. This leads to an equal weights\n    " \
    "assumption for the prior ensemble at each step.\n    " \
    "This function should be called at each model time step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
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
    "        5. py__g2l_obs_pdaf (localise each ensemble member\n    " \
    "           in observation space)\n    " \
    "        6. py__likelihood_l_pdaf\n    " \
    "        7. core DA algorithm\n    " \
    "        8. py__l2g_state_pdaf\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Tödter, J., and B. Ahrens, 2015:\n    "\
    "       A second-order exact ensemble square root filter\n    " \
    "       for nonlinear data assimilation. Mon. Wea. Rev.,\n    " \
    "       143, 1347–1367, doi:10.1175/MWR-D-14-00108.1."
docstrings['put_state_lknetf'] = \
    "It is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.localomi_put_state`\n    "\
    "or :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`.\n\n    "\
    "PDAF-OMI modules require fewer user-supplied functions\n    " \
    "and improved efficiency.\n\n    " \
    "A hybridised LETKF and LNETF [1]_ for a single DA step.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_lknetf`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the analysis is not post-processed,\n    " \
    "and distributed to the model forecast\n    " \
    "by user-supplied functions. The next DA step will\n    " \
    "not be assigned by user-supplied functions as well.\n    " \
    "This function is typically used when there are\n    " \
    "not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The LNETF computes the distribution up to\n    " \
    "the second moment similar to Kalman filters but\n    " \
    "using a nonlinear weighting similar to\n    " \
    "particle filters. This leads to an equal weights\n    " \
    "assumption for the prior ensemble.\n    " \
    "The hybridisation with LETKF is expected to lead to\n    " \
    "improved performance for\n    " \
    "quasi-Gaussian problems.\n    " \
    "The function should be called at each model step.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_n_domains_p_pdaf\n    " \
    "    4. py__init_dim_obs_pdaf\n    " \
    "    5. py__obs_op_pdaf\n    "\
    "       (for each ensemble member)\n    " \
    "    6. py__init_obs_pdaf\n    " \
    "       (if global adaptive forgetting factor `type_forget=1`\n    "\
    "       is used in :func:`pyPDAF.PDAF.init`)\n    " \
    "    7. py__init_obsvar_pdaf (if global adaptive\n    " \
    "       forgetting factor is used)\n    " \
    "    8. loop over each local domain:\n    " \
    "        1. py__init_dim_l_pdaf\n    " \
    "        2. py__init_dim_obs_l_pdaf\n    " \
    "        3. py__g2l_state_pdaf\n    " \
    "        4. py__g2l_obs_pdaf\n    "\
    "           (localise each ensemble member in observation space)\n    " \
    "        5. py__init_obs_l_pdaf\n    "\
    "        6. py__init_obsvar_l_pdaf\n    "\
    "           (only called if local adaptive forgetting factor\n    " \
    "           `type_forget=2` is used)\n    "\
    "        7. py__prodRinvA_pdaf\n    " \
    "        8. py__likelihood_l_pdaf\n    " \
    "        9. core DA algorithm\n    " \
    "        10. py__l2g_state_pdaf\n    " \
    "    9. py__obs_op_pdaf\n    " \
    "       (only called with `HKN` and `HNK` options called\n    " \
    "       for each ensemble member)\n    " \
    "    10. py__likelihood_hyb_l_pda\n    " \
    "    11. py__init_obsvar_l_pdaf\n    " \
    "        (only called if local adaptive forgetting factor\n    " \
    "        `type_forget=2` is used)\n    "\
    "    12. py__prodRinvA_hyb_l_pdaf\n" \
    "\n    " \
    ".. deprecated:: 1.0.0\n\n    " \
    "   This function is replaced by\n    " \
    "   :func:`pyPDAF.PDAF.localomi_put_state`\n    " \
    "   and :func:`pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR`" \
    "\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Nerger, L.. (2022) \n    " \
    "       Data assimilation for nonlinear systems with\n    " \
    "       a hybrid nonlinear Kalman ensemble transform filter.\n    " \
    "       Q J R Meteorol Soc, 620–640. doi:10.1002/qj.4221"
docstrings['put_state_prepost'] = \
    "It is used to preprocess and postprocess of the ensemble.\n\n    " \
    "No DA is performed in this function.\n    " \
    "Compared to :func:`pyPDAF.PDAF.assimilate_prepost`,\n    " \
    "this function does not set assimilation flag, \n    " \
    "and does not distribute the processed ensemble to the model field.\n    " \
    "This function also does not set the next assimilation step as\n    " \
    ":func:`pyPDAF.PDAF.assimilate_prepost`\n    " \
    "because it does not call :func:`pyPDAF.PDAF.get_state`.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf (preprocess, step < 0)"
docstrings['put_state_generate_obs'] = \
    "Generation of synthetic observations\n    " \
    "based on given error statistics and observation operator\n    " \
    "without post-processing, distributing analysis,\n    " \
    "and setting next observation step.\n\n    " \
    "When diagonal observation error covariance matrix is used,\n    " \
    "it is recommended to use\n    " \
    ":func:`pyPDAF.PDAF.omi_generate_obs` functionalities\n    "\
    "for fewer user-supplied functions and improved efficiency.\n\n    " \
    "The generated synthetic observations are\n    " \
    "based on each member of model forecast.\n    " \
    "Therefore, an ensemble of observations can be obtained.\n    " \
    "In a typical experiment,\n    "\
    "one may only need one ensemble member.\n\n    " \
    "Compared to :func:`pyPDAF.PDAF.generate_obs`,\n    " \
    "this function has no :func:`get_state` call.\n    " \
    "This means that the next DA step will\n    " \
    "not be assigned by user-supplied functions.\n    " \
    "This function is typically used when there\n    " \
    "are not enough CPUs to run the ensemble in parallel,\n    "\
    "and some ensemble members have to be run serially.\n    " \
    "The :func:`pyPDAF.PDAF.get_state` function follows this\n    "\
    "function call to ensure the sequential DA.\n\n    " \
    "The implementation strategy is similar to\n    " \
    "an assimilation step. This means that, \n    " \
    "one can reuse many user-supplied functions for\n    " \
    "assimilation and observation generation.\n\n    " \
    "User-supplied functions are executed in the following sequence:\n    " \
    "    1. py__collect_state_pdaf\n    " \
    "    2. py__prepoststep_state_pdaf\n    " \
    "    3. py__init_dim_obs_pdaf\n    " \
    "    4. py__obs_op_pda\n    " \
    "    5. py__init_obserr_f_pdaf\n    " \
    "    6. py__get_obs_f_pdaf"
