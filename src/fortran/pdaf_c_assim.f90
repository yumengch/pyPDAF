MODULE pdaf_c_assim
use iso_c_binding, only: c_int, c_double, c_bool
use PDAF
! use PDAF_analysis_utils
use pdaf_c_cb_interface
use pdaf_c_f_interface

implicit none

contains
   SUBROUTINE c__PDAF_get_state(steps, time, doexit, u_next_observation,  &
      u_distribute_state, u_prepoststep, outflag) bind(c)
      ! Flag and number of time steps
      INTEGER(c_int), INTENT(inout) :: steps
      ! current model time1
      REAL(c_double), INTENT(out) :: time
      ! Whether to exit from forecasts
      INTEGER(c_int), INTENT(inout) :: doexit
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: u_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      next_observation_pdaf_c_ptr => u_next_observation
      distribute_state_pdaf_c_ptr => u_distribute_state
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_get_state(steps, time, doexit, f__next_observation_pdaf,  &
         f__distribute_state_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_get_state

   SUBROUTINE c__PDAF_assimilate_estkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_prodrinva,  &
      u_init_obsvar, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_estkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, &
         f__prepoststep_pdaf, f__prodrinva_pdaf,  &
         f__init_obsvar_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_estkf

   SUBROUTINE c__PDAF_assim_offline_estkf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_prodrinva, u_init_obsvar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar

      call PDAF_assim_offline_estkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__prodrinva_pdaf, f__init_obsvar_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_estkf

   SUBROUTINE c__PDAF_assimilate_3dvar(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva, u_cvt, u_cvt_adj,  &
      u_obs_op_lin, u_obs_op_adj, u_prepoststep, u_next_observation,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_3dvar(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, &
         f__prodrinva_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf, &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_3dvar

   SUBROUTINE c__PDAF_assim_offline_3dvar(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj, u_prepoststep,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_3dvar(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prodrinva_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_3dvar

   SUBROUTINE c__PDAF_assimilate_en3dvar_lestkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f,  &
      u_obs_op_f, u_init_obs_f, u_init_obs_l, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_prepoststep,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: u_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obs_f_pdaf_c_ptr => u_init_obs_f
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_en3dvar_lestkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, &
         f__prodrinva_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, &
         f__init_dim_obs_f_pdaf,  &
         f__obs_op_f_pdaf, f__init_obs_f_pdaf, f__init_obs_l_pdaf, &
         f__prodrinva_l_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, &
         f__g2l_state_pdaf,  &
         f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, &
         f__init_obsvar_l_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAF_assim_offline_en3dvar_lestkf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: u_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obs_f_pdaf_c_ptr => u_init_obs_f
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_en3dvar_lestkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, &
         f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__init_dim_obs_f_pdaf, f__obs_op_f_pdaf, &
         f__init_obs_f_pdaf,  &
         f__init_obs_l_pdaf, f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, &
         f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, &
         f__g2l_obs_pdaf, f__init_obsvar_pdaf,  &
         f__init_obsvar_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_en3dvar_lestkf

   SUBROUTINE c__PDAF_assimilate_ensrf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obsvars,  &
      u_localize_covar_serial, u_prepoststep, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Initialize vector of observation error variances
      procedure(c__init_obsvars_pdaf) :: u_init_obsvars
      ! Apply localization for single-observation vectors
      procedure(c__localize_covar_serial_pdaf) :: u_localize_covar_serial
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obsvars_pdaf_c_ptr => u_init_obsvars
      localize_covar_serial_pdaf_c_ptr => u_localize_covar_serial
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_ensrf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obsvars_pdaf,  &
         f__localize_covar_serial_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_ensrf

   SUBROUTINE c__PDAF_assim_offline_ensrf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_init_obsvars, u_localize_covar_serial, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Initialize vector of observation error variances
      procedure(c__init_obsvars_pdaf) :: u_init_obsvars
      ! Apply localization for single-observation vectors
      procedure(c__localize_covar_serial_pdaf) :: u_localize_covar_serial
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obsvars_pdaf_c_ptr => u_init_obsvars
      localize_covar_serial_pdaf_c_ptr => u_localize_covar_serial
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_ensrf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obsvars_pdaf, f__localize_covar_serial_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_ensrf

   SUBROUTINE c__PDAF_assimilate_lknetf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
      u_prodrinva_l, u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar,  &
      u_init_obsvar_l, u_likelihood_l, u_likelihood_hyb_l, u_next_observation,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: u_prodrinva_hyb_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Compute likelihood
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l
      ! Compute likelihood with hybrid weight
      procedure(c__likelihood_hyb_l_pdaf) :: u_likelihood_hyb_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      prodrinva_hyb_l_pdaf_c_ptr => u_prodrinva_hyb_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      likelihood_l_pdaf_c_ptr => u_likelihood_l
      likelihood_hyb_l_pdaf_c_ptr => u_likelihood_hyb_l
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_lknetf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__prodrinva_l_pdaf, f__prodrinva_hyb_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf,  &
         f__init_obsvar_l_pdaf, f__likelihood_l_pdaf, f__likelihood_hyb_l_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_lknetf

   SUBROUTINE c__PDAF_assim_offline_lknetf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_likelihood_l, u_likelihood_hyb_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: u_prodrinva_hyb_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Compute likelihood
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l
      ! Compute likelihood with hybrid weight
      procedure(c__likelihood_hyb_l_pdaf) :: u_likelihood_hyb_l

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      prodrinva_hyb_l_pdaf_c_ptr => u_prodrinva_hyb_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      likelihood_l_pdaf_c_ptr => u_likelihood_l
      likelihood_hyb_l_pdaf_c_ptr => u_likelihood_hyb_l

      call PDAF_assim_offline_lknetf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obs_l_pdaf, f__prepoststep_pdaf, f__prodrinva_l_pdaf, f__prodrinva_hyb_l_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf,  &
         f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, f__init_obsvar_l_pdaf,  &
         f__likelihood_l_pdaf, f__likelihood_hyb_l_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_lknetf

   SUBROUTINE c__PDAF_assimilate_hyb3dvar_estkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
      u_init_obsvar, u_prepoststep, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_hyb3dvar_estkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_obsvar_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_hyb3dvar_estkf

   SUBROUTINE c__PDAF_assim_offline_hyb3dvar_estkf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens,  &
      u_obs_op_lin, u_obs_op_adj, u_init_obsvar, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_hyb3dvar_estkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_obsvar_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_hyb3dvar_estkf

   SUBROUTINE c__PDAF_assimilate_hyb3dvar_lestkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
      u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_prepoststep,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: u_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obs_f_pdaf_c_ptr => u_init_obs_f
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_hyb3dvar_lestkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_dim_obs_f_pdaf, f__obs_op_f_pdaf, f__init_obs_f_pdaf, f__init_obs_l_pdaf,  &
         f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, f__init_obsvar_l_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAF_assim_offline_hyb3dvar_lestkf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_cvt, u_cvt_adj,  &
      u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f, u_init_obs_f,  &
      u_init_obs_l, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar,  &
      u_init_obsvar_l, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: u_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      cvt_pdaf_c_ptr => u_cvt
      cvt_adj_pdaf_c_ptr => u_cvt_adj
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obs_f_pdaf_c_ptr => u_init_obs_f
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_hyb3dvar_lestkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_dim_obs_f_pdaf, f__obs_op_f_pdaf,  &
         f__init_obs_f_pdaf, f__init_obs_l_pdaf, f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf,  &
         f__init_obsvar_pdaf, f__init_obsvar_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_hyb3dvar_lestkf

   SUBROUTINE c__PDAF_assimilate_lestkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_lestkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, f__init_obsvar_l_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_lestkf

   SUBROUTINE c__PDAF_assim_offline_lestkf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l

      call PDAF_assim_offline_lestkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obs_l_pdaf, f__prepoststep_pdaf, f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf,  &
         f__init_obsvar_pdaf, f__init_obsvar_l_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_lestkf

   SUBROUTINE c__PDAF_assimilate_enkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_add_obs_error,  &
      u_init_obs_covar, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: u_add_obs_error
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      add_obs_err_pdaf_c_ptr => u_add_obs_error
      init_obs_covar_pdaf_c_ptr => u_init_obs_covar
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_enkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__add_obs_err_pdaf,  &
         f__init_obs_covar_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_enkf

   SUBROUTINE c__PDAF_assim_offline_enkf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_add_obs_err, u_init_obs_covar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      add_obs_err_pdaf_c_ptr => u_add_obs_err
      init_obs_covar_pdaf_c_ptr => u_init_obs_covar

      call PDAF_assim_offline_enkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__add_obs_err_pdaf, f__init_obs_covar_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_enkf

   SUBROUTINE c__PDAF_assimilate_letkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_letkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, f__init_obsvar_l_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_letkf

   SUBROUTINE c__PDAF_assim_offline_letkf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_init_obs_l, u_prepoststep, u_prodrinva_l, u_init_n_domains_p,  &
      u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
      u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l

      call PDAF_assim_offline_letkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obs_l_pdaf, f__prepoststep_pdaf, f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf,  &
         f__init_obsvar_pdaf, f__init_obsvar_l_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_letkf

   SUBROUTINE c__PDAF_assimilate_seik(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_prodrinva,  &
      u_init_obsvar, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_seik(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__prodrinva_pdaf,  &
         f__init_obsvar_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_seik

   SUBROUTINE c__PDAF_assim_offline_seik(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_prodrinva, u_init_obsvar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar

      call PDAF_assim_offline_seik(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__prodrinva_pdaf, f__init_obsvar_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_seik

   SUBROUTINE c__PDAF_assimilate_lnetf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
      u_likelihood_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_l_pdaf_c_ptr => u_likelihood_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_lnetf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__likelihood_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_lnetf

   SUBROUTINE c__PDAF_assim_offline_lnetf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_init_obs_l, u_prepoststep, u_likelihood_l, u_init_n_domains_p,  &
      u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_l_pdaf_c_ptr => u_likelihood_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs

      call PDAF_assim_offline_lnetf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obs_l_pdaf, f__prepoststep_pdaf, f__likelihood_l_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf,  &
         outflag)

   END SUBROUTINE c__PDAF_assim_offline_lnetf

   SUBROUTINE c__PDAF_assimilate_prepost(u_collect_state, u_distribute_state,  &
      u_prepoststep, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_prepost(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_prepost

   SUBROUTINE c__PDAF_assimilate_en3dvar_estkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, u_init_obsvar,  &
      u_prepoststep, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_en3dvar_estkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_obsvar_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_en3dvar_estkf

   SUBROUTINE c__PDAF_assim_offline_en3dvar_estkf(u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, u_init_obsvar, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prodrinva_pdaf_c_ptr => u_prodrinva
      cvt_ens_pdaf_c_ptr => u_cvt_ens
      cvt_adj_ens_pdaf_c_ptr => u_cvt_adj_ens
      obs_op_lin_pdaf_c_ptr => u_obs_op_lin
      obs_op_adj_pdaf_c_ptr => u_obs_op_adj
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_assim_offline_en3dvar_estkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_obs_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__init_obsvar_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_en3dvar_estkf

   SUBROUTINE c__PDAF_assimilate_netf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_likelihood,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_pdaf_c_ptr => u_likelihood
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_netf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__likelihood_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_netf

   SUBROUTINE c__PDAF_assim_offline_netf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_likelihood, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_pdaf_c_ptr => u_likelihood

      call PDAF_assim_offline_netf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__likelihood_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_netf

   SUBROUTINE c__PDAF_assimilate_pf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_likelihood,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_pdaf_c_ptr => u_likelihood
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_pf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__likelihood_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_pf

   SUBROUTINE c__PDAF_assim_offline_pf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_likelihood, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      likelihood_pdaf_c_ptr => u_likelihood

      call PDAF_assim_offline_pf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__likelihood_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_pf

   SUBROUTINE c__PDAF_assimilate_lenkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_localize,  &
      u_add_obs_error, u_init_obs_covar, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: u_localize
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: u_add_obs_error
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      localize_covar_pdaf_c_ptr => u_localize
      add_obs_err_pdaf_c_ptr => u_add_obs_error
      init_obs_covar_pdaf_c_ptr => u_init_obs_covar
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_lenkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__localize_covar_pdaf,  &
         f__add_obs_err_pdaf, f__init_obs_covar_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_lenkf

   SUBROUTINE c__PDAF_assim_offline_lenkf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_localize, u_add_obs_err, u_init_obs_covar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: u_localize
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      localize_covar_pdaf_c_ptr => u_localize
      add_obs_err_pdaf_c_ptr => u_add_obs_err
      init_obs_covar_pdaf_c_ptr => u_init_obs_covar

      call PDAF_assim_offline_lenkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__localize_covar_pdaf, f__add_obs_err_pdaf, f__init_obs_covar_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_lenkf

   SUBROUTINE c__PDAF_assimilate_etkf(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prepoststep, u_prodrinva,  &
      u_init_obsvar, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_etkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__prepoststep_pdaf, f__prodrinva_pdaf,  &
         f__init_obsvar_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_etkf

   SUBROUTINE c__PDAF_assim_offline_etkf(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prepoststep, u_prodrinva, u_init_obsvar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_pdaf_c_ptr => u_prodrinva
      init_obsvar_pdaf_c_ptr => u_init_obsvar

      call PDAF_assim_offline_etkf(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__prepoststep_pdaf, f__prodrinva_pdaf, f__init_obsvar_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_etkf

   SUBROUTINE c__PDAF_assimilate_lseik(u_collect_state, u_distribute_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_assimilate_lseik(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf, f__init_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf, f__init_obsvar_pdaf, f__init_obsvar_l_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_assimilate_lseik

   SUBROUTINE c__PDAF_assim_offline_lseik(u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_init_obs_l, u_prepoststep, u_prodrinva_l, u_init_n_domains_p,  &
      u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
      u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l

      init_dim_obs_pdaf_c_ptr => u_init_dim_obs
      obs_op_pdaf_c_ptr => u_obs_op
      init_obs_pdaf_c_ptr => u_init_obs
      init_obs_l_pdaf_c_ptr => u_init_obs_l
      prepoststep_pdaf_c_ptr => u_prepoststep
      prodrinva_l_pdaf_c_ptr => u_prodrinva_l
      init_n_domains_p_pdaf_c_ptr => u_init_n_domains_p
      init_dim_l_pdaf_c_ptr => u_init_dim_l
      init_dim_obs_l_pdaf_c_ptr => u_init_dim_obs_l
      g2l_state_pdaf_c_ptr => u_g2l_state
      l2g_state_pdaf_c_ptr => u_l2g_state
      g2l_obs_pdaf_c_ptr => u_g2l_obs
      init_obsvar_pdaf_c_ptr => u_init_obsvar
      init_obsvar_l_pdaf_c_ptr => u_init_obsvar_l

      call PDAF_assim_offline_lseik(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_obs_pdaf,  &
         f__init_obs_l_pdaf, f__prepoststep_pdaf, f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__g2l_obs_pdaf,  &
         f__init_obsvar_pdaf, f__init_obsvar_l_pdaf, outflag)

   END SUBROUTINE c__PDAF_assim_offline_lseik

   SUBROUTINE c__PDAF_generate_obs(u_collect_state, u_distribute_state,  &
      u_init_dim_obs_f, u_obs_op_f, u_init_obserr_f, u_get_obs_f,  &
      u_prepoststep, u_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize vector of observation error standard deviations
      procedure(c__init_obserr_f_pdaf) :: u_init_obserr_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: u_get_obs_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state
      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obserr_f_pdaf_c_ptr => u_init_obserr_f
      get_obs_f_pdaf_c_ptr => u_get_obs_f
      prepoststep_pdaf_c_ptr => u_prepoststep
      next_observation_pdaf_c_ptr => u_next_observation

      call PDAF_generate_obs(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_f_pdaf, f__obs_op_f_pdaf, f__init_obserr_f_pdaf, f__get_obs_f_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF_generate_obs

   SUBROUTINE c__PDAF_generate_obs_offline(u_init_dim_obs_f, u_obs_op_f,  &
      u_init_obserr_f, u_get_obs_f, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Initialize vector of observation error standard deviations
      procedure(c__init_obserr_f_pdaf) :: u_init_obserr_f
      ! Provide observation vector
      procedure(c__get_obs_f_pdaf) :: u_get_obs_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      init_dim_obs_f_pdaf_c_ptr => u_init_dim_obs_f
      obs_op_f_pdaf_c_ptr => u_obs_op_f
      init_obserr_f_pdaf_c_ptr => u_init_obserr_f
      get_obs_f_pdaf_c_ptr => u_get_obs_f
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF_generate_obs_offline(f__init_dim_obs_f_pdaf, f__obs_op_f_pdaf,  &
         f__init_obserr_f_pdaf, f__get_obs_f_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_generate_obs_offline
END MODULE pdaf_c_assim
