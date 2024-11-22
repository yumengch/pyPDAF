API
===
.. contents::
   :local:
   :depth: 2

Initialisation and finalisation
-------------------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.init
   pyPDAF.PDAF.omi_init
   pyPDAF.PDAF.deallocate
   pyPDAF.PDAF.print_info

`Fully parallel` DA algorithms
------------------------------

Sequential DA
^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_assimilate_global
   pyPDAF.PDAF.omi_assimilate_lenkf
   pyPDAF.PDAF.localomi_assimilate

non-diagnoal observation matrix
"""""""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_assimilate_global_nondiagR
   pyPDAF.PDAF.omi_assimilate_enkf_nondiagR
   pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR
   pyPDAF.PDAF.localomi_assimilate_nondiagR
   pyPDAF.PDAF.omi_assimilate_nonlin_nondiagR
   pyPDAF.PDAF.localomi_assimilate_lnetf_nondiagR
   pyPDAF.PDAF.localomi_assimilate_lknetf_nondiagR

Variational DA
^^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_assimilate_3dvar
   pyPDAF.PDAF.omi_assimilate_en3dvar_estkf
   pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf
   pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf
   pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf

non-diagnoal observation matrix
"""""""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR
   pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR
   pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR
   pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR
   pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR


`Flexible` DA algorithms
------------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:


Sequential DA
^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_put_state_global
   pyPDAF.PDAF.omi_put_state_lenkf
   pyPDAF.PDAF.localomi_put_state


diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_put_state_global_nondiagR
   pyPDAF.PDAF.omi_put_state_enkf_nondiagR
   pyPDAF.PDAF.omi_put_state_lenkf_nondiagR
   pyPDAF.PDAF.localomi_put_state_nondiagR
   pyPDAF.PDAF.omi_put_state_nonlin_nondiagR
   pyPDAF.PDAF.localomi_put_state_lnetf_nondiagR
   pyPDAF.PDAF.localomi_put_state_lknetf_nondiagR

Variational DA
^^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_put_state_3dvar
   pyPDAF.PDAF.omi_put_state_en3dvar_estkf
   pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf
   pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf
   pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf

non-diagnoal observation matrix
"""""""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_put_state_3dvar_nondiagR
   pyPDAF.PDAF.omi_put_state_en3dvar_estkf_nondiagR
   pyPDAF.PDAF.omi_put_state_hyb3dvar_estkf_nondiagR
   pyPDAF.PDAF.localomi_put_state_en3dvar_lestkf_nondiagR
   pyPDAF.PDAF.localomi_put_state_hyb3dvar_lestkf_nondiagR



OMI functions
-------------

obs_f setter functions
^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_set_doassim
   pyPDAF.PDAF.omi_set_disttype
   pyPDAF.PDAF.omi_set_ncoord
   pyPDAF.PDAF.omi_set_id_obs_p
   pyPDAF.PDAF.omi_set_icoeff_p
   pyPDAF.PDAF.omi_set_domainsize
   pyPDAF.PDAF.omi_set_obs_err_type
   pyPDAF.PDAF.omi_set_use_global_obs
   pyPDAF.PDAF.omi_set_inno_omit
   pyPDAF.PDAF.omi_set_inno_omit_ivar

Observation operators
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_obs_op_gridpoint
   pyPDAF.PDAF.omi_obs_op_gridavg
   pyPDAF.PDAF.omi_obs_op_interp_lin
   pyPDAF.PDAF.omi_obs_op_adj_gridavg
   pyPDAF.PDAF.omi_obs_op_adj_gridpoint
   pyPDAF.PDAF.omi_obs_op_adj_interp_lin

Interpolations
^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_get_interp_coeff_tri
   pyPDAF.PDAF.omi_get_interp_coeff_lin1D
   pyPDAF.PDAF.omi_get_interp_coeff_lin

Localisation
^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_obsstats_l
   pyPDAF.PDAF.omi_weights_l
   pyPDAF.PDAF.omi_weights_l_sgnl
   pyPDAF.PDAF.omi_init_dim_obs_l_iso
   pyPDAF.PDAF.omi_init_dim_obs_l_noniso
   pyPDAF.PDAF.omi_init_dim_obs_l_noniso_locweights
   pyPDAF.PDAF.omi_localize_covar_iso
   pyPDAF.PDAF.omi_localize_covar_noniso
   pyPDAF.PDAF.omi_localize_covar_noniso_locweights

Others
^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_set_domain_limits
   pyPDAF.PDAF.omi_set_debug_flag
   pyPDAF.PDAF.omi_gather_obs
   pyPDAF.PDAF.omi_gather_obsstate
   pyPDAF.PDAF.omi_gather_obsdims
   pyPDAF.PDAF.omi_obsstats
   pyPDAF.PDAF.omit_obs_omi
   pyPDAF.PDAF.omi_deallocate_obs
   pyPDAF.PDAF.omi_check_error

Local module functions
----------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.local_set_indices
   pyPDAF.PDAF.local_set_increment_weights
   pyPDAF.PDAF.local_clear_increment_weights

Utilities
---------

Synthetic experiments
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_generate_obs
   pyPDAF.PDAF.generate_obs

   pyPDAF.PDAF.omi_put_state_generate_obs
   pyPDAF.PDAF.put_state_generate_obs

Statistical diagnostics
^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.diag_effsample
   pyPDAF.PDAF.diag_ensstats
   pyPDAF.PDAF.diag_histogram
   pyPDAF.PDAF.diag_CRPS
   pyPDAF.PDAF.diag_CRPS_nompi

Ensemble generation
^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.eofcovar
   pyPDAF.PDAF.SampleEns

Diagnostics
^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.get_assim_flag
   pyPDAF.PDAF.get_ensstats
   pyPDAF.PDAF.get_localfilter
   pyPDAF.PDAF.get_memberid
   pyPDAF.PDAF.get_obsmemberid
   pyPDAF.PDAF.get_smootherens
   pyPDAF.PDAF.set_debug_flag
   pyPDAF.PDAF.print_local_obsstats
   pyPDAF.PDAF.print_domain_stats

Advanced manipulation
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.set_ens_pointer
   pyPDAF.PDAF.set_smootherens
   pyPDAF.PDAF.set_memberid
   pyPDAF.PDAF.set_comm_pdaf
   pyPDAF.PDAF.set_offline_mode
   pyPDAF.PDAF.reset_forget

   pyPDAF.PDAF.init_local_obsstats
   pyPDAF.PDAF.incr_local_obsstats
   pyPDAF.PDAF.force_analysis

Others
^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.get_state
   pyPDAF.PDAF.assimilate_prepost
   pyPDAF.PDAF.prepost
   pyPDAF.PDAF.put_state_prepost

   pyPDAF.PDAF.gather_dim_obs_f
   pyPDAF.PDAF.gather_obs_f
   pyPDAF.PDAF.gather_obs_f2

   pyPDAF.PDAF.incremental
   pyPDAF.PDAF.add_increment
   pyPDAF.PDAF.local_weight
   pyPDAF.PDAF.local_weights

   pyPDAF.PDAF.gather_obs_f2_flex
   pyPDAF.PDAF.gather_obs_f_flex

Internal matrix operations
--------------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.seik_TtimesA
   pyPDAF.PDAF.etkf_Tleft
   pyPDAF.PDAF.estkf_OmegaA
   pyPDAF.PDAF.enkf_omega
   pyPDAF.PDAF.seik_omega


Internal callback functions
---------------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.omi_init_obs_f_cb
   pyPDAF.PDAF.omi_init_obsvar_cb
   pyPDAF.PDAF.omi_g2l_obs_cb
   pyPDAF.PDAF.omi_init_obs_l_cb
   pyPDAF.PDAF.omi_init_obsvar_l_cb
   pyPDAF.PDAF.omi_prodRinvA_l_cb
   pyPDAF.PDAF.omi_likelihood_l_cb
   pyPDAF.PDAF.omi_prodRinvA_cb
   pyPDAF.PDAF.omi_likelihood_cb
   pyPDAF.PDAF.omi_add_obs_error_cb
   pyPDAF.PDAF.omi_init_obscovar_cb
   pyPDAF.PDAF.omi_init_obserr_f_cb
   pyPDAF.PDAF.omi_prodRinvA_hyb_l_cb
   pyPDAF.PDAF.omi_likelihood_hyb_l_cb
   pyPDAF.PDAF.omi_omit_by_inno_l_cb
   pyPDAF.PDAF.omi_omit_by_inno_cb
   pyPDAF.PDAF.local_g2l_cb
   pyPDAF.PDAF.local_l2g_cb
