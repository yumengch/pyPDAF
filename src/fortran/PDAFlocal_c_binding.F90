MODULE PDAFlocal_c_binding
use iso_c_binding, only: c_double, c_int
use U_PDAF_interface_c_binding
implicit none

contains

   SUBROUTINE c__PDAFlocal_set_indices(dim_l, map) bind(c)
      ! Dimension of local state vector
      integer(c_int), INTENT(in) :: dim_l
      ! Index array for mapping between local and global state vector
      integer(c_int), INTENT(in) :: map(dim_l)

      call PDAFlocal_set_indices(dim_l, map)
   end subroutine c__PDAFlocal_set_indices

   subroutine c__PDAFlocal_set_increment_weights(dim_l, weights) bind(c)
      ! Dimension of local state vector
      integer(c_int), INTENT(in) :: dim_l
      ! Weights array
      real(c_double), INTENT(in) :: weights(dim_l)

      call PDAFlocal_set_increment_weights(dim_l, weights)
   end subroutine c__PDAFlocal_set_increment_weights

   subroutine c__PDAFlocal_clear_increment_weights() bind(c)
      call PDAFlocal_clear_increment_weights()
   end subroutine c__PDAFlocal_clear_increment_weights


   subroutine c__PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l) bind(c)
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! PE-local full state dimension
      integer(c_int), INTENT(in) :: dim_p
      ! Local state dimension
      integer(c_int), INTENT(in) :: dim_l
      ! PE-local full state vector
      real(c_double), INTENT(in)    :: state_p(dim_p)
      ! State vector on local analysis domain
      real(c_double), INTENT(out)   :: state_l(dim_l)

      call PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l)
   end subroutine c__PDAFlocal_g2l_cb

   subroutine c__PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p) bind(c)
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Local state dimension
      integer(c_int), INTENT(in) :: dim_l
      ! PE-local full state dimension
      integer(c_int), INTENT(in) :: dim_p
      ! State vector on local analysis domain
      real(c_double), INTENT(in)    :: state_l(dim_l)
      ! PE-local full state vector
      real(c_double), INTENT(inout) :: state_p(dim_p)

      call PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p)
   end subroutine c__PDAFlocal_l2g_cb

   subroutine c__PDAFlocalomi_assimilate(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag

      call PDAFlocalomi_assimilate(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate

   subroutine c__PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
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
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag

      call PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_en3dvar_lestkf

   subroutine c__PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: prodRinvA_pdafomi
      ! Provide product R^-1 A with localization
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag

      call PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR

   subroutine c__PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
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
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag

      call PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_hyb3dvar_lestkf

   subroutine c__PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: prodRinvA_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR

   subroutine c__PDAFlocalomi_assimilate_lknetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
      likelihood_l_pdafomi, likelihood_hyb_l_pdafomi,  &
      next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodRinvA_hyb_l_pdaf) :: prodRinvA_hyb_l_pdafomi
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdaf) :: likelihood_hyb_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_assimilate_lknetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
      likelihood_l_pdafomi, likelihood_hyb_l_pdafomi,  &
      next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_lknetf_nondiagR

   subroutine c__PDAFlocalomi_assimilate_lnetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi,  &
      next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_assimilate_lnetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi,  &
      next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_lnetf_nondiagR

   subroutine c__PDAFlocalomi_assimilate_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
      next_observation_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product of inverse of R with matrix A
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_assimilate_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
      next_observation_pdaf, outflag)
   end subroutine c__PDAFlocalomi_assimilate_nondiagR

   subroutine c__PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      outflag)
   end subroutine c__PDAFlocalomi_put_state

   subroutine c__PDAFlocalomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
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
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, outflag)
   end subroutine c__PDAFlocalomi_put_state_en3dvar_lestkf

   subroutine c__PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: prodRinvA_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, outflag)
   end subroutine c__PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR

   subroutine c__PDAFlocalomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
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
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, outflag)
   end subroutine c__PDAFlocalomi_put_state_hyb3dvar_lestkf

   subroutine c__PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
      obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: prodRinvA_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
      obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
      prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, outflag)
   end subroutine c__PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR

   subroutine c__PDAFlocalomi_put_state_lknetf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
      likelihood_l_pdafomi, likelihood_hyb_l_pdafomi,  &
      outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodRinvA_hyb_l_pdaf) :: prodRinvA_hyb_l_pdafomi
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdaf) :: likelihood_hyb_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_lknetf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
      likelihood_l_pdafomi, likelihood_hyb_l_pdafomi,  &
      outflag)
   end subroutine c__PDAFlocalomi_put_state_lknetf_nondiagR

   subroutine c__PDAFlocalomi_put_state_lnetf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi,  &
      outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdaf) :: likelihood_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_lnetf_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi,  &
      outflag)
   end subroutine c__PDAFlocalomi_put_state_lnetf_nondiagR

   subroutine c__PDAFlocalomi_put_state_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
      outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdafomi
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product of inverse of R with matrix A
      procedure(c__prodRinvA_l_pdaf) :: prodRinvA_l_pdafomi
      ! Status flag
      integer(c_int), INTENT(inout) :: outflag
      call PDAFlocalomi_put_state_nondiagR(collect_state_pdaf, &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
      outflag)
   end subroutine c__PDAFlocalomi_put_state_nondiagR

   subroutine c__PDAFlocal_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l,  &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: U_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: U_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: U_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l,  &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_en3dvar_lestkf

   subroutine c__PDAFlocal_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l,  &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: U_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: U_cvt_adj_ens
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: U_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: U_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: U_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l,  &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_hyb3dvar_lestkf

   subroutine c__PDAFlocal_assimilate_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_lestkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_lestkf

   subroutine c__PDAFlocal_assimilate_letkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_letkf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_letkf

   subroutine c__PDAFlocal_assimilate_lknetf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_prodRinvA_hyb_l, U_init_n_domains_p, U_init_dim_l, &
      U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_likelihood_l, U_likelihood_hyb_l, &
      U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Provide product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodRinvA_hyb_l_pdaf) :: U_prodRinvA_hyb_l
      ! Compute likelihood
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! Compute likelihood with hybrid weight
      procedure(c__likelihood_hyb_l_pdaf) :: U_likelihood_hyb_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_lknetf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_prodRinvA_hyb_l, U_init_n_domains_p, U_init_dim_l, &
      U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_likelihood_l, U_likelihood_hyb_l, &
      U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_lknetf

   subroutine c__PDAFlocal_assimilate_lnetf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
      U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_lnetf(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
      U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_lnetf

   subroutine c__PDAFlocal_assimilate_lseik(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide number of forecast time steps until
      ! next assimilations, model physical time and
      ! end of assimilation cycles
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! distribute a state vector from pdaf to the model/any arrays
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_assimilate_lseik(U_collect_state, U_distribute_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
      U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_next_observation, outflag)
   end subroutine c__PDAFlocal_assimilate_lseik


   subroutine c__PDAFlocal_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: U_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: U_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: U_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag)
   end subroutine c__PDAFlocal_put_state_en3dvar_lestkf

   subroutine c__PDAFlocal_put_state_hyb3dvar_lestkf(U_collect_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: U_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: U_cvt_adj_ens
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: U_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: U_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_f_pdaf) :: U_init_obs_f
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_hyb3dvar_lestkf(U_collect_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag)
   end subroutine c__PDAFlocal_put_state_hyb3dvar_lestkf

   SUBROUTINE c__PDAFlocal_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag)
   end subroutine c__PDAFlocal_put_state_lestkf

   subroutine c__PDAFlocal_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag)
   end subroutine c__PDAFlocal_put_state_letkf

   subroutine c__PDAFlocal_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_prodRinvA_hyb_l, &
      U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, U_likelihood_l, U_likelihood_hyb_l, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Provide product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodRinvA_hyb_l_pdaf) :: U_prodRinvA_hyb_l
      ! Compute likelihood
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! Compute likelihood with hybrid weight
      procedure(c__likelihood_hyb_l_pdaf) :: U_likelihood_hyb_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_prodRinvA_hyb_l, &
      U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, U_likelihood_l, U_likelihood_hyb_l, outflag)
   end subroutine c__PDAFlocal_put_state_lknetf

   subroutine c__PDAFlocal_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      outflag)
   end subroutine c__PDAFlocal_put_state_lnetf

   SUBROUTINE c__PDAFlocal_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag) bind(c)
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: U_init_obs_l
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: U_init_obsvar_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! Preprocesse the ensemble before analysis
      ! and postprocess the ensemble before
      ! distributing to the model for next forecast
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Status flag
      integer(c_int), INTENT(out) :: outflag
      call PDAFlocal_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
      U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
      U_init_obsvar, U_init_obsvar_l, outflag)
   END SUBROUTINE c__PDAFlocal_put_state_lseik


END MODULE PDAFlocal_c_binding
