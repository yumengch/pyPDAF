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

   pyPDAF.init
   pyPDAF.set_parallel
   pyPDAF.init_forecast
   pyPDAF.PDAFomi.init
   pyPDAF.PDAFomi.init_local
   pyPDAF.deallocate

DA algorithms
------------------------------

Sequential DA
^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate
   pyPDAF.assim_offline

non-diagnoal observation matrix
"""""""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate_local_nondiagr
   pyPDAF.assimilate_global_nondiagr
   pyPDAF.assimilate_lnetf_nondiagr
   pyPDAF.assimilate_lknetf_nondiagr
   pyPDAF.assimilate_enkf_nondiagr
   pyPDAF.assimilate_nonlin_nondiagr

   pyPDAF.assim_offline_local_nondiagr
   pyPDAF.assim_offline_global_nondiagr
   pyPDAF.assim_offline_lnetf_nondiagr
   pyPDAF.assim_offline_lknetf_nondiagr
   pyPDAF.assim_offline_enkf_nondiagr
   pyPDAF.assim_offline_lenkf_nondiagr
   pyPDAF.assim_offline_nonlin_nondiagr

Variational DA
^^^^^^^^^^^^^^

diagnoal observation matrix
"""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate_3dvar_all
   pyPDAF.assim_offline_3dvar_all

non-diagnoal observation matrix
"""""""""""""""""""""""""""""""
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate_3dvar_nondiagr
   pyPDAF.assimilate_en3dvar_estkf_nondiagr
   pyPDAF.assimilate_en3dvar_lestkf_nondiagr
   pyPDAF.assimilate_hyb3dvar_estkf_nondiagr
   pyPDAF.assimilate_hyb3dvar_lestkf_nondiagr

   pyPDAF.assim_offline_3dvar_nondiagr
   pyPDAF.assim_offline_en3dvar_estkf_nondiagr
   pyPDAF.assim_offline_en3dvar_lestkf_nondiagr
   pyPDAF.assim_offline_hyb3dvar_estkf_nondiagr
   pyPDAF.assim_offline_hyb3dvar_lestkf_nondiagr


OMI functions
-------------

setter functions
^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.set_doassim
   pyPDAF.PDAFomi.set_disttype
   pyPDAF.PDAFomi.set_ncoord
   pyPDAF.PDAFomi.set_obs_err_type
   pyPDAF.PDAFomi.set_use_global_obs
   pyPDAF.PDAFomi.set_inno_omit
   pyPDAF.PDAFomi.set_inno_omit_ivar
   pyPDAF.PDAFomi.set_id_obs_p
   pyPDAF.PDAFomi.set_icoeff_p
   pyPDAF.PDAFomi.set_domainsize
   pyPDAF.PDAFomi.set_name
   pyPDAF.PDAFomi.gather_obs


Observation operators
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.obs_op_gridpoint
   pyPDAF.PDAFomi.obs_op_gridavg
   pyPDAF.PDAFomi.obs_op_extern
   pyPDAF.PDAFomi.obs_op_interp_lin
   pyPDAF.PDAFomi.obs_op_adj_gridavg
   pyPDAF.PDAFomi.obs_op_adj_gridpoint
   pyPDAF.PDAFomi.obs_op_adj_interp_lin
   pyPDAF.PDAFomi.gather_obsstate

Interpolations
^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.get_interp_coeff_tri
   pyPDAF.PDAFomi.get_interp_coeff_lin1d
   pyPDAF.PDAFomi.get_interp_coeff_lin

Localisation
^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.init_dim_obs_l_iso
   pyPDAF.PDAFomi.init_dim_obs_l_noniso
   pyPDAF.PDAFomi.init_dim_obs_l_noniso_locweights
   pyPDAF.PDAFomi.observation_localization_weights
   pyPDAF.PDAFomi.set_domain_limits
   pyPDAF.PDAFomi.get_domain_limits_unstr
   pyPDAF.PDAFomi.set_localize_covar_iso
   pyPDAF.PDAFomi.set_localize_covar_noniso
   pyPDAF.PDAFomi.set_localize_covar_noniso_locweights

Custom local observation initialisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.set_localization
   pyPDAF.PDAFomi.set_localization_noniso
   pyPDAF.PDAFomi.set_dim_obs_l
   pyPDAF.PDAFomi.store_obs_l_index
   pyPDAF.PDAFomi.store_obs_l_index_vdist

Diagnostics
^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.check_error
   pyPDAF.PDAFomi.set_debug_flag
   pyPDAF.PDAFomi.set_obs_diag
   pyPDAF.PDAFomi.check_error
   pyPDAF.PDAFomi.diag_dimobs
   pyPDAF.PDAFomi.diag_get_hx
   pyPDAF.PDAFomi.diag_get_hxmean
   pyPDAF.PDAFomi.diag_get_ivar
   pyPDAF.PDAFomi.diag_get_obs
   pyPDAF.PDAFomi.diag_nobstypes
   pyPDAF.PDAFomi.diag_obs_rmsd
   pyPDAF.PDAFomi.diag_stats


Localisation functions
----------------------
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFlocal.set_indices
   pyPDAF.PDAFlocal.set_increment_weights
   pyPDAF.PDAFlocal.clear_increment_weights
   pyPDAF.PDAF.correlation_function
   pyPDAF.PDAF.local_weight
   pyPDAF.PDAF.local_weights

Utilities
---------

PDAF state and setup information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.get_assim_flag
   pyPDAF.PDAF.get_localfilter
   pyPDAF.PDAF.get_local_type
   pyPDAF.PDAF.get_memberid
   pyPDAF.PDAF.get_obsmemberid
   pyPDAF.PDAF.get_smoother_ens
   pyPDAF.PDAF.print_filter_types
   pyPDAF.PDAF.print_da_types
   pyPDAF.PDAF.print_info


Observation MPI handling
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.gather_dim_obs_f
   pyPDAF.PDAF.gather_obs_f
   pyPDAF.PDAF.gather_obs_f2
   pyPDAF.PDAF.gather_obs_f_flex
   pyPDAF.PDAF.gather_obs_f2_flex


Synthetic experiments
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.generate_obs
   pyPDAF.generate_obs_offline


Incremental analysis update
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.iau_init
   pyPDAF.PDAF.iau_reset
   pyPDAF.PDAF.iau_set_pointer

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
   pyPDAF.PDAF.sample_ens

PDAF debug options
^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.set_debug_flag

Advanced manipulation
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.set_iparam
   pyPDAF.PDAF.set_rparam
   pyPDAF.PDAF.set_comm_pdaf

   pyPDAF.PDAF.set_ens_pointer
   pyPDAF.PDAF.set_memberid
   pyPDAF.PDAF.set_offline_mode
   pyPDAF.PDAF.set_seedset
   pyPDAF.PDAF.set_smoother_ens

   pyPDAF.PDAF.force_analysis
   pyPDAF.PDAF.reset_forget

