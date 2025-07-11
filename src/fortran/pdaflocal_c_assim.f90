MODULE pdaflocal_c_assim
use iso_c_binding, only: c_double, c_int, c_bool, c_loc, c_char, c_null_char
use PDAF
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAFlocal_assimilate_en3dvar_lestkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f,  &
      u_obs_op_f, u_init_obs_f, u_init_obs_l, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_obs,  &
      u_init_obsvar, u_init_obsvar_l, u_prepoststep, u_next_observation,  &
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

      call PDAFlocal_assimilate_en3dvar_lestkf(u_collect_state,  &
         u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj,  &
         u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_prepoststep,  &
         u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAFlocal_assimilate_hyb3dvar_lestkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
      u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_obs,  &
      u_init_obsvar, u_init_obsvar_l, u_prepoststep, u_next_observation,  &
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

      call PDAFlocal_assimilate_hyb3dvar_lestkf(u_collect_state,  &
         u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_cvt_ens, u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin,  &
         u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f, u_init_obs_f,  &
         u_init_obs_l, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
         u_init_dim_obs_l, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_prepoststep, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAFlocal_assimilate_lseik(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prepoststep, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
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
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      call PDAFlocal_assimilate_lseik(u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_lseik

   SUBROUTINE c__PDAFlocal_assimilate_letkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prepoststep, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
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
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      call PDAFlocal_assimilate_letkf(u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_letkf

   SUBROUTINE c__PDAFlocal_assimilate_lestkf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prepoststep, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
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
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      call PDAFlocal_assimilate_lestkf(u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_lestkf

   SUBROUTINE c__PDAFlocal_assimilate_lnetf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prepoststep, u_likelihood_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_obs, u_next_observation, outflag) bind(c)
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
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Routine to provide time step, time and dimensionof next observation
      procedure(c__next_observation_pdaf) :: u_next_observation

      call PDAFlocal_assimilate_lnetf(u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
         u_likelihood_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_obs, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_lnetf

   SUBROUTINE c__PDAFlocal_assimilate_lknetf(u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prepoststep, u_prodrinva_l, u_prodrinva_hyb_l, u_init_n_domains_p,  &
      u_init_dim_l, u_init_dim_obs_l, u_g2l_obs, u_init_obsvar,  &
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

      call PDAFlocal_assimilate_lknetf(u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep,  &
         u_prodrinva_l, u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l,  &
         u_init_dim_obs_l, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_likelihood_l, u_likelihood_hyb_l, u_next_observation, outflag)

   END SUBROUTINE c__PDAFlocal_assimilate_lknetf
END MODULE pdaflocal_c_assim
