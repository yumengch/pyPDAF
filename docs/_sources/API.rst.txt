API Reference
=============

This page lists the functions exported by the public ``pyPDAF``
``__init__.py`` files. The top-level namespace contains the functions most
users need for a PDAF3 workflow. The subpackages expose lower-level PDAF,
PDAFomi, localisation, diagnostic, and compatibility routines.

The function names follow PDAF closely. See :doc:`naming_convention` for the
meaning of suffixes such as ``_p``, ``_l``, and ``_f``.

.. contents::
   :local:
   :depth: 2

Top-Level Workflow Functions
----------------------------

Initialisation and finalisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.init
   pyPDAF.init_parallel
   pyPDAF.init_forecast
   pyPDAF.set_parallel
   pyPDAF.deallocate

Assimilation drivers
^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate
   pyPDAF.assim_offline
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

Variational and hybrid drivers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.assimilate_3dvar_all
   pyPDAF.assim_offline_3dvar_all
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

Observation generation and hooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.generate_obs
   pyPDAF.generate_obs_offline
   pyPDAF.prepost
   pyPDAF.prepost_offline

Top-level utilities
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.get_fcst_info
   pyPDAF.print_version
   pyPDAF.print_filter_types
   pyPDAF.configinfo_filters
   pyPDAF.options_filters
   pyPDAF.flush_fortran_stdout
   pyPDAF.global_except_hook

``pyPDAF.PDAF`` Utilities and Legacy Interface
----------------------------------------------

Core setup, execution, and information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.init
   pyPDAF.PDAF.init_forecast
   pyPDAF.PDAF.deallocate
   pyPDAF.PDAF.finalize
   pyPDAF.PDAF.abort
   pyPDAF.PDAF.get_state
   pyPDAF.PDAF.get_fcst_info
   pyPDAF.PDAF.force_analysis
   pyPDAF.PDAF.print_version
   pyPDAF.PDAF.print_filter_types
   pyPDAF.PDAF.print_da_types
   pyPDAF.PDAF.print_info
   pyPDAF.PDAF.configinfo_filters
   pyPDAF.PDAF.options_filters
   pyPDAF.PDAF.flush_fortran_stdout

Observation gathering
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.gather_dim_obs_f
   pyPDAF.PDAF.gather_obs_f
   pyPDAF.PDAF.gather_obs_f2
   pyPDAF.PDAF.gather_obs_f_flex
   pyPDAF.PDAF.gather_obs_f2_flex

Localisation helpers
^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.correlation_function
   pyPDAF.PDAF.local_weight
   pyPDAF.PDAF.local_weights

Random numbers and ensemble generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.generate_rndvec
   pyPDAF.PDAF.eofcovar
   pyPDAF.PDAF.sample_ens
   pyPDAF.PDAF.reset_forget

Setters
^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.set_comm_pdaf
   pyPDAF.PDAF.set_debug_flag
   pyPDAF.PDAF.set_ens_pointer
   pyPDAF.PDAF.set_iparam
   pyPDAF.PDAF.set_memberid
   pyPDAF.PDAF.set_offline_mode
   pyPDAF.PDAF.set_rparam
   pyPDAF.PDAF.genobs_set_rparam
   pyPDAF.PDAF.set_seedset
   pyPDAF.PDAF.set_seed
   pyPDAF.PDAF.set_seedvec
   pyPDAF.PDAF.set_smoother_ens

Getters and state flags
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.get_assim_flag
   pyPDAF.PDAF.get_localfilter
   pyPDAF.PDAF.get_local_type
   pyPDAF.PDAF.get_memberid
   pyPDAF.PDAF.get_obsmemberid
   pyPDAF.PDAF.get_seed
   pyPDAF.PDAF.get_seedvec
   pyPDAF.PDAF.get_rndcount
   pyPDAF.PDAF.reset_fcst_flag
   pyPDAF.PDAF.get_smoother_ens

Incremental analysis update
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.iau_init
   pyPDAF.PDAF.iau_reset
   pyPDAF.PDAF.iau_set_pointer

Diagnostics
^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF.diag_ensmean
   pyPDAF.PDAF.diag_stddev_nompi
   pyPDAF.PDAF.diag_stddev
   pyPDAF.PDAF.diag_variance_nompi
   pyPDAF.PDAF.diag_variance
   pyPDAF.PDAF.diag_rmsd_nompi
   pyPDAF.PDAF.diag_rmsd
   pyPDAF.PDAF.diag_crps_mpi
   pyPDAF.PDAF.diag_crps_nompi
   pyPDAF.PDAF.diag_crps
   pyPDAF.PDAF.diag_effsample
   pyPDAF.PDAF.diag_ensstats
   pyPDAF.PDAF.diag_compute_moments
   pyPDAF.PDAF.diag_histogram
   pyPDAF.PDAF.diag_reliability_budget
   pyPDAF.PDAF.diag_diffstats

``pyPDAF.PDAF3`` Assimilation Interface
---------------------------------------

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAF3.init
   pyPDAF.PDAF3.init_parallel
   pyPDAF.PDAF3.init_forecast
   pyPDAF.PDAF3.set_parallel
   pyPDAF.PDAF3.assimilate
   pyPDAF.PDAF3.assim_offline
   pyPDAF.PDAF3.assimilate_3dvar_all
   pyPDAF.PDAF3.assim_offline_3dvar_all
   pyPDAF.PDAF3.assimilate_local_nondiagr
   pyPDAF.PDAF3.assimilate_global_nondiagr
   pyPDAF.PDAF3.assimilate_lnetf_nondiagr
   pyPDAF.PDAF3.assimilate_lknetf_nondiagr
   pyPDAF.PDAF3.assimilate_enkf_nondiagr
   pyPDAF.PDAF3.assimilate_nonlin_nondiagr
   pyPDAF.PDAF3.assimilate_3dvar_nondiagr
   pyPDAF.PDAF3.assimilate_en3dvar_estkf_nondiagr
   pyPDAF.PDAF3.assimilate_en3dvar_lestkf_nondiagr
   pyPDAF.PDAF3.assimilate_hyb3dvar_estkf_nondiagr
   pyPDAF.PDAF3.assimilate_hyb3dvar_lestkf_nondiagr
   pyPDAF.PDAF3.assim_offline_local_nondiagr
   pyPDAF.PDAF3.assim_offline_global_nondiagr
   pyPDAF.PDAF3.assim_offline_lnetf_nondiagr
   pyPDAF.PDAF3.assim_offline_lknetf_nondiagr
   pyPDAF.PDAF3.assim_offline_enkf_nondiagr
   pyPDAF.PDAF3.assim_offline_lenkf_nondiagr
   pyPDAF.PDAF3.assim_offline_nonlin_nondiagr
   pyPDAF.PDAF3.assim_offline_3dvar_nondiagr
   pyPDAF.PDAF3.assim_offline_en3dvar_estkf_nondiagr
   pyPDAF.PDAF3.assim_offline_en3dvar_lestkf_nondiagr
   pyPDAF.PDAF3.assim_offline_hyb3dvar_estkf_nondiagr
   pyPDAF.PDAF3.assim_offline_hyb3dvar_lestkf_nondiagr
   pyPDAF.PDAF3.generate_obs
   pyPDAF.PDAF3.generate_obs_offline
   pyPDAF.PDAF3.prepost
   pyPDAF.PDAF3.prepost_offline

``pyPDAF.PDAFomi`` Observation Module Interface
-----------------------------------------------

Setup and error handling
^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.init
   pyPDAF.PDAFomi.check_error
   pyPDAF.PDAFomi.set_debug_flag

Observation storage and metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.gather_obs
   pyPDAF.PDAFomi.gather_obsstate
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
   pyPDAF.PDAFomi.set_searchtype

Interpolation coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.get_interp_coeff_tri
   pyPDAF.PDAFomi.get_interp_coeff_lin1d
   pyPDAF.PDAFomi.get_interp_coeff_lin
   pyPDAF.PDAFomi.get_interp_coeff_tri_vec
   pyPDAF.PDAFomi.get_interp_coeff_lin1d_vec
   pyPDAF.PDAFomi.get_interp_coeff_lin_vec

Observation operators
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.obs_op_gridpoint
   pyPDAF.PDAFomi.obs_op_gridavg
   pyPDAF.PDAFomi.obs_op_extern
   pyPDAF.PDAFomi.obs_op_interp_lin
   pyPDAF.PDAFomi.obs_op_adj_gridpoint
   pyPDAF.PDAFomi.obs_op_adj_gridavg
   pyPDAF.PDAFomi.obs_op_adj_interp_lin

Local observations and localisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.init_dim_obs_l_iso
   pyPDAF.PDAFomi.init_dim_obs_l_noniso
   pyPDAF.PDAFomi.init_dim_obs_l_noniso_locweights
   pyPDAF.PDAFomi.observation_localization_weights
   pyPDAF.PDAFomi.set_dim_obs_l
   pyPDAF.PDAFomi.set_localization
   pyPDAF.PDAFomi.set_localization_noniso
   pyPDAF.PDAFomi.set_localize_covar_iso
   pyPDAF.PDAFomi.set_localize_covar_noniso
   pyPDAF.PDAFomi.set_localize_covar_noniso_locweights
   pyPDAF.PDAFomi.set_domain_limits
   pyPDAF.PDAFomi.get_domain_limits_unstr
   pyPDAF.PDAFomi.store_obs_l_index
   pyPDAF.PDAFomi.store_obs_l_index_vdist

Observation diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFomi.set_obs_diag
   pyPDAF.PDAFomi.diag_dimobs
   pyPDAF.PDAFomi.diag_get_hx
   pyPDAF.PDAFomi.diag_get_hxmean
   pyPDAF.PDAFomi.diag_get_ivar
   pyPDAF.PDAFomi.diag_get_obs
   pyPDAF.PDAFomi.diag_nobstypes
   pyPDAF.PDAFomi.diag_obs_rmsd
   pyPDAF.PDAFomi.diag_stats
   pyPDAF.PDAFomi.diag_rmsd
   pyPDAF.PDAFomi.diag_diffstats
   pyPDAF.PDAFomi.diag_crps

``pyPDAF.PDAFlocal`` Localisation Helpers
-----------------------------------------

.. autosummary::
   :toctree: _autosummary
   :recursive:

   pyPDAF.PDAFlocal.set_indices
   pyPDAF.PDAFlocal.set_increment_weights
   pyPDAF.PDAFlocal.clear_increment_weights

``pyPDAF.PDAFlocalomi``
-----------------------

``pyPDAF.PDAFlocalomi`` is exported as a namespace package, but it currently
does not expose functions from its ``__init__.py``.
