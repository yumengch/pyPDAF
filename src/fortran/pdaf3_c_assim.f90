MODULE pdaf3_c_assim
use iso_c_binding, only: c_double, c_int, c_bool
use PDAF
use pdaf_c_cb_interface
use pdaf_c_f_interface

implicit none

contains
   SUBROUTINE c__PDAF_get_fcst_info(steps, time, doexit) bind(c)
      ! Flag and number of time steps
      INTEGER(c_int), INTENT(inout) :: steps
      ! current model time
      REAL(c_double), INTENT(inout) :: time
      ! Whether to exit from forecasts
      INTEGER(c_int), INTENT(inout) :: doexit


      call PDAF_get_fcst_info(steps, time, doexit)

   END SUBROUTINE c__PDAF_get_fcst_info

   SUBROUTINE c__PDAF3_assimilate_3dvar_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_3dvar_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_3dvar_nondiagR

   SUBROUTINE c__PDAF3_assimilate_en3dvar_estkf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_en3dvar_estkf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_en3dvar_estkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      prodrinva_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide product R^-1 A and apply localizations
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_en3dvar_lestkf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__prodrinva_l_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_en3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_hyb3dvar_estkf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_hyb3dvar_estkf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_hyb3dvar_estkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prodrinva_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide product R^-1 A and apply localizations
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_hyb3dvar_lestkf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prodrinva_l_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_hyb3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_3dvar_all(init_dim_obs_pdaf, obs_op_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_3dvar_all(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_3dvar_all

   SUBROUTINE c__PDAF3_assim_offline_3dvar(init_dim_obs_pdaf, obs_op_pdaf,  &
      cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_3dvar(f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_pdaf,  &
         f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_3dvar

   SUBROUTINE c__PDAF3_assim_offline_en3dvar(init_dim_obs_pdaf, obs_op_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_en3dvar(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_en3dvar

   SUBROUTINE c__PDAF3_assim_offline_en3dvar_estkf(init_dim_obs_pdaf,  &
      obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_en3dvar_estkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_en3dvar_estkf

   SUBROUTINE c__PDAF3_assim_offline_en3dvar_lestkf(init_dim_obs_pdaf,  &
      obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_en3dvar_lestkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_en3dvar_lestkf

   SUBROUTINE c__PDAF3_assim_offline_hyb3dvar(init_dim_obs_pdaf, obs_op_pdaf,  &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_hyb3dvar(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_hyb3dvar

   SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_estkf(init_dim_obs_pdaf,  &
      obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_hyb3dvar_estkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_estkf

   SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_lestkf(init_dim_obs_pdaf,  &
      obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_hyb3dvar_lestkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_lestkf

   SUBROUTINE c__PDAF3_assim_offline_3dvar_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prodrinva_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_3dvar_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prodrinva_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_3dvar_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_en3dvar_estkf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prodrinva_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_en3dvar_estkf_nondiagR(f__init_dim_obs_pdaf,  &
         f__obs_op_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_en3dvar_estkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_en3dvar_lestkf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prodrinva_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prodrinva_l_pdaf, init_n_domains_pdaf,  &
      init_dim_l_pdaf, init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide product R^-1 A and apply localizations
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_en3dvar_lestkf_nondiagR(f__init_dim_obs_pdaf,  &
         f__obs_op_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prodrinva_l_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_en3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_estkf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prodrinva_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf,  &
      cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_hyb3dvar_estkf_nondiagR(f__init_dim_obs_pdaf,  &
         f__obs_op_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf,  &
         f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_estkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_lestkf_nondiagR( &
      init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prodrinva_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide product R^-1 A and apply localizations
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_hyb3dvar_lestkf_nondiagR(f__init_dim_obs_pdaf,  &
         f__obs_op_pdaf, f__prodrinva_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf, f__cvt_pdaf,  &
         f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prodrinva_l_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_hyb3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate(collect_state_pdaf, distribute_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate

   SUBROUTINE c__PDAF3_assimilate_local(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Get local state from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      g2l_state_pdaf_c_ptr => g2l_state_pdaf
      l2g_state_pdaf_c_ptr => l2g_state_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_local(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__g2l_state_pdaf, f__l2g_state_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_local

   SUBROUTINE c__PDAF3_assimilate_global(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_global(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_global

   SUBROUTINE c__PDAF3_assimilate_lenkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, localize_pdaf,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: localize_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      localize_covar_pdaf_c_ptr => localize_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_lenkf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__localize_covar_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_lenkf

   SUBROUTINE c__PDAF3_assimilate_ensrf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      localize_serial_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply localization to HP and BXY for single observation
      procedure(c__localize_covar_serial_pdaf) :: localize_serial_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      localize_covar_serial_pdaf_c_ptr => localize_serial_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_ensrf(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__localize_covar_serial_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_ensrf

   SUBROUTINE c__PDAF3_assim_offline(init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline

   SUBROUTINE c__PDAF3_assim_offline_local(init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Get local state from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      g2l_state_pdaf_c_ptr => g2l_state_pdaf
      l2g_state_pdaf_c_ptr => l2g_state_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_local(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__g2l_state_pdaf, f__l2g_state_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_local

   SUBROUTINE c__PDAF3_assim_offline_global(init_dim_obs_pdaf, obs_op_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_global(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_global

   SUBROUTINE c__PDAF3_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf,  &
      localize_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: localize_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      localize_covar_pdaf_c_ptr => localize_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_lenkf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__localize_covar_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_lenkf

   SUBROUTINE c__PDAF3_assim_offline_ensrf(init_dim_obs_pdaf, obs_op_pdaf,  &
      localize_serial_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply localization to HP and BXY for single observation
      procedure(c__localize_covar_serial_pdaf) :: localize_serial_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      localize_covar_serial_pdaf_c_ptr => localize_serial_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_ensrf(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__localize_covar_serial_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_ensrf

   SUBROUTINE c__PDAF3_assim_offline_local_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      prodrinva_l_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_local_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prodrinva_l_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_local_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_global_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prodrinva_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_global_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prodrinva_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_global_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_lnetf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, likelihood_l_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      likelihood_l_pdaf_c_ptr => likelihood_l_pdaf

      call PDAF3_assim_offline_lnetf_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prepoststep_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__likelihood_l_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_lnetf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_lknetf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prodrinva_l_pdaf, prodrinva_hyb_l_pdaf,  &
      likelihood_l_pdaf, likelihood_hyb_l_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: prodrinva_hyb_l_pdaf
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdaf
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdaf) :: likelihood_hyb_l_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      prodrinva_hyb_l_pdaf_c_ptr => prodrinva_hyb_l_pdaf
      likelihood_l_pdaf_c_ptr => likelihood_l_pdaf
      likelihood_hyb_l_pdaf_c_ptr => likelihood_hyb_l_pdaf

      call PDAF3_assim_offline_lknetf_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prepoststep_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__prodrinva_l_pdaf, f__prodrinva_hyb_l_pdaf,  &
         f__likelihood_l_pdaf, f__likelihood_hyb_l_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_lknetf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_enkf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: add_obs_error_pdaf
      ! Initialize mean observation error variance
      procedure(c__init_obs_covar_pdaf) :: init_obscovar_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      add_obs_err_pdaf_c_ptr => add_obs_error_pdaf
      init_obs_covar_pdaf_c_ptr => init_obscovar_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_enkf_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__add_obs_err_pdaf, f__init_obs_covar_pdaf, f__prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAF3_assim_offline_enkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_lenkf_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, prepoststep_pdaf, localize_pdaf, add_obs_error_pdaf,  &
      init_obscovar_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply covariance localization
      procedure(c__localize_covar_pdaf) :: localize_pdaf
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: add_obs_error_pdaf
      ! Initialize mean observation error variance
      procedure(c__init_obs_covar_pdaf) :: init_obscovar_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      localize_covar_pdaf_c_ptr => localize_pdaf
      add_obs_err_pdaf_c_ptr => add_obs_error_pdaf
      init_obs_covar_pdaf_c_ptr => init_obscovar_pdaf

      call PDAF3_assim_offline_lenkf_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__prepoststep_pdaf, f__localize_covar_pdaf, f__add_obs_err_pdaf,  &
         f__init_obs_covar_pdaf, outflag)
   END SUBROUTINE c__PDAF3_assim_offline_lenkf_nondiagR

   SUBROUTINE c__PDAF3_assim_offline_nonlin_nondiagR(init_dim_obs_pdaf,  &
      obs_op_pdaf, likelihood_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Compute likelihood
      procedure(c__likelihood_pdaf) :: likelihood_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      likelihood_pdaf_c_ptr => likelihood_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_assim_offline_nonlin_nondiagR(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__likelihood_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assim_offline_nonlin_nondiagR

   SUBROUTINE c__PDAF3_assimilate_3dvar_all(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_3dvar_all(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_3dvar_all

   SUBROUTINE c__PDAF3_assimilate_3dvar(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_pdaf,  &
      cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_3dvar(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_3dvar

   SUBROUTINE c__PDAF3_assimilate_en3dvar(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf,  &
      init_dim_l_pdaf, init_dim_obs_l_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_en3dvar(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf,  &
         f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf,  &
         f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_en3dvar

   SUBROUTINE c__PDAF3_assimilate_en3dvar_estkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_en3dvar_estkf(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_en3dvar_estkf

   SUBROUTINE c__PDAF3_assimilate_en3dvar_lestkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf,  &
      init_dim_l_pdaf, init_dim_obs_l_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_en3dvar_lestkf(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAF3_assimilate_hyb3dvar(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_hyb3dvar(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf, f__cvt_adj_ens_pdaf,  &
         f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf, f__obs_op_adj_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_hyb3dvar

   SUBROUTINE c__PDAF3_assimilate_hyb3dvar_estkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_hyb3dvar_estkf(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_hyb3dvar_estkf

   SUBROUTINE c__PDAF3_assimilate_hyb3dvar_lestkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf,  &
      obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      cvt_ens_pdaf_c_ptr => cvt_ens_pdaf
      cvt_adj_ens_pdaf_c_ptr => cvt_adj_ens_pdaf
      cvt_pdaf_c_ptr => cvt_pdaf
      cvt_adj_pdaf_c_ptr => cvt_adj_pdaf
      obs_op_lin_pdaf_c_ptr => obs_op_lin_pdaf
      obs_op_adj_pdaf_c_ptr => obs_op_adj_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_hyb3dvar_lestkf(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__cvt_ens_pdaf,  &
         f__cvt_adj_ens_pdaf, f__cvt_pdaf, f__cvt_adj_pdaf, f__obs_op_lin_pdaf,  &
         f__obs_op_adj_pdaf, f__init_n_domains_p_pdaf, f__init_dim_l_pdaf,  &
         f__init_dim_obs_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAF3_assimilate_local_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      prodrinva_l_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_local_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prodrinva_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_local_nondiagR

   SUBROUTINE c__PDAF3_assimilate_global_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prodrinva_pdaf_c_ptr => prodrinva_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_global_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prodrinva_pdaf,  &
         f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_global_nondiagR

   SUBROUTINE c__PDAF3_assimilate_lnetf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      likelihood_l_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      likelihood_l_pdaf_c_ptr => likelihood_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_lnetf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__likelihood_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_lnetf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_lknetf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      prodrinva_l_pdaf, prodrinva_hyb_l_pdaf, likelihood_l_pdaf,  &
      likelihood_hyb_l_pdaf, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdaf
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: prodrinva_hyb_l_pdaf
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdaf
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdaf) :: likelihood_hyb_l_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      init_n_domains_p_pdaf_c_ptr => init_n_domains_pdaf
      init_dim_l_pdaf_c_ptr => init_dim_l_pdaf
      init_dim_obs_l_pdaf_c_ptr => init_dim_obs_l_pdaf
      prodrinva_l_pdaf_c_ptr => prodrinva_l_pdaf
      prodrinva_hyb_l_pdaf_c_ptr => prodrinva_hyb_l_pdaf
      likelihood_l_pdaf_c_ptr => likelihood_l_pdaf
      likelihood_hyb_l_pdaf_c_ptr => likelihood_hyb_l_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_lknetf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__init_n_domains_p_pdaf, f__init_dim_l_pdaf, f__init_dim_obs_l_pdaf,  &
         f__prodrinva_l_pdaf, f__prodrinva_hyb_l_pdaf, f__likelihood_l_pdaf,  &
         f__likelihood_hyb_l_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_lknetf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_enkf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,  &
      add_obs_error_pdaf, init_obscovar_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: add_obs_error_pdaf
      ! Initialize mean observation error variance
      procedure(c__init_obs_covar_pdaf) :: init_obscovar_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      add_obs_err_pdaf_c_ptr => add_obs_error_pdaf
      init_obs_covar_pdaf_c_ptr => init_obscovar_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_enkf_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__add_obs_err_pdaf, f__init_obs_covar_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAF3_assimilate_enkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_lenkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, localize_pdaf,  &
      add_obs_error_pdaf, init_obscovar_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply covariance localization
      procedure(c__localize_covar_pdaf) :: localize_pdaf
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: add_obs_error_pdaf
      ! Initialize mean observation error variance
      procedure(c__init_obs_covar_pdaf) :: init_obscovar_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      localize_covar_pdaf_c_ptr => localize_pdaf
      add_obs_err_pdaf_c_ptr => add_obs_error_pdaf
      init_obs_covar_pdaf_c_ptr => init_obscovar_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_lenkf_nondiagR(f__collect_state_pdaf,  f__distribute_state_pdaf, &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__prepoststep_pdaf, f__localize_covar_pdaf,  &
         f__add_obs_err_pdaf, f__init_obs_covar_pdaf, f__next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAF3_assimilate_lenkf_nondiagR

   SUBROUTINE c__PDAF3_assimilate_nonlin_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, likelihood_pdaf,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Compute likelihood
      procedure(c__likelihood_pdaf) :: likelihood_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      likelihood_pdaf_c_ptr => likelihood_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_assimilate_nonlin_nondiagR(f__collect_state_pdaf,  &
         f__distribute_state_pdaf, f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__likelihood_pdaf, f__prepoststep_pdaf, f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_assimilate_nonlin_nondiagR

   SUBROUTINE c__PDAF3_generate_obs(collect_state_pdaf, distribute_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, get_obs_pdaf, prepoststep_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Initialize observation vector
      procedure(c__get_obs_f_pdaf) :: get_obs_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      collect_state_pdaf_c_ptr => collect_state_pdaf
      distribute_state_pdaf_c_ptr => distribute_state_pdaf
      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      get_obs_f_pdaf_c_ptr  => get_obs_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf
      next_observation_pdaf_c_ptr => next_observation_pdaf

      call PDAF3_generate_obs(f__collect_state_pdaf, f__distribute_state_pdaf,  &
         f__init_dim_obs_pdaf, f__obs_op_pdaf, f__get_obs_f_pdaf, f__prepoststep_pdaf,  &
         f__next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAF3_generate_obs

   SUBROUTINE c__PDAF3_generate_obs_offline(init_dim_obs_pdaf, obs_op_pdaf,  &
      get_obs_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Initialize observation vector
      procedure(c__get_obs_f_pdaf) :: get_obs_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf

      init_dim_obs_pdaf_c_ptr => init_dim_obs_pdaf
      obs_op_pdaf_c_ptr => obs_op_pdaf
      get_obs_f_pdaf_c_ptr  => get_obs_pdaf
      prepoststep_pdaf_c_ptr => prepoststep_pdaf

      call PDAF3_generate_obs_offline(f__init_dim_obs_pdaf, f__obs_op_pdaf,  &
         f__get_obs_f_pdaf, f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF3_generate_obs_offline
end MODULE pdaf3_c_assim
