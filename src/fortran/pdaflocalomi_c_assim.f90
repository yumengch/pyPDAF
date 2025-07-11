
MODULE pdaflocalomi_c_assim
use iso_c_binding, only: c_double, c_int, c_bool, c_loc, c_char, c_null_char
use PDAF3
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAFlocalomi_assimilate(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf,  &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate(collect_state_pdaf, distribute_state_pdaf,  &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
         next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate

   SUBROUTINE c__PDAFlocalomi_assimilate_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdafomi, prodrinva_l_pdafomi, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdafomi
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_nondiagR(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
         prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdafomi, prodrinva_l_pdafomi, next_observation_pdaf,  &
         outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_nondiagR

   SUBROUTINE c__PDAFlocalomi_assimilate_lnetf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdafomi, likelihood_l_pdafomi, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_lnetf_nondiagR(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
         prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdafomi, likelihood_l_pdafomi, next_observation_pdaf,  &
         outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_lnetf_nondiagR

   SUBROUTINE c__PDAFlocalomi_assimilate_lknetf_nondiagR(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdafomi, prodrinva_l_pdafomi, prodrinva_hyb_l_pdafomi,  &
      likelihood_l_pdafomi, likelihood_hyb_l_pdafomi, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdafomi
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: prodrinva_hyb_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdaf) :: likelihood_hyb_l_pdafomi
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_lknetf_nondiagR(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
         prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdafomi, prodrinva_l_pdafomi, prodrinva_hyb_l_pdafomi,  &
         likelihood_l_pdafomi, likelihood_hyb_l_pdafomi, next_observation_pdaf,  &
         outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_lknetf_nondiagR

   SUBROUTINE c__PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf,  &
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
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
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
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf,  &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
         prepoststep_pdaf, next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf,  &
      distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf,  &
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
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
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
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf,  &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf,  &
         obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf,  &
         init_dim_l_pdaf, init_dim_obs_l_pdaf, prepoststep_pdaf,  &
         next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR( &
      collect_state_pdaf, distribute_state_pdaf, init_dim_obs_pdafomi,  &
      obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      obs_op_lin_pdafomi, obs_op_adj_pdafomi, prodrinva_l_pdafomi,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Provide product R^-1 A with localization
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdafomi
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
         prodrinva_pdafomi, cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi,  &
         obs_op_adj_pdafomi, prodrinva_l_pdafomi, init_n_domains_pdaf,  &
         init_dim_l_pdaf, init_dim_obs_l_pdafomi, prepoststep_pdaf,  &
         next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR( &
      collect_state_pdaf, distribute_state_pdaf, init_dim_obs_pdafomi,  &
      obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
      prodrinva_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdafomi, prepoststep_pdaf, next_observation_pdaf,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_l_pdaf) :: prodrinva_l_pdafomi
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf

      call PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf,  &
         distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,  &
         prodrinva_pdafomi, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf,  &
         cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
         prodrinva_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdafomi, prepoststep_pdaf, next_observation_pdaf, outflag)

   END SUBROUTINE c__PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR
END MODULE pdaflocalomi_c_assim
