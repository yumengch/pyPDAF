MODULE PDAF_c_binding
use iso_c_binding, only: c_double, c_int, c_bool, c_ptr, c_loc, c_char, c_null_char
use PDAF_analysis_utils
use U_PDAF_interface_c_binding

implicit none

contains

   SUBROUTINE c__PDAF_assimilate_3dvar(U_collect_state, U_distribute_state, &
                                       U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                       U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                       U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: U_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: U_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj

      call PDAF_assimilate_3dvar(U_collect_state, U_distribute_state, &
                                       U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                       U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                       U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_assimilate_3dvar

   SUBROUTINE c__PDAF_assimilate_en3dvar_estkf(U_collect_state, U_distribute_state, &
                                               U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                               U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
                                               U_init_obsvar, U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: U_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: U_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj

      call PDAF_assimilate_en3dvar_estkf(U_collect_state, U_distribute_state, &
                                               U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                               U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
                                               U_init_obsvar, U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_assimilate_en3dvar_estkf

   SUBROUTINE c__PDAF_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
                                                U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
                                                U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
                                                U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
                                                U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
                                                U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l

      call PDAF_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
                                                U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
                                                U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
                                                U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
                                                U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
                                                U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAF_assimilate_enkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_add_obs_error, &
         U_init_obs_covar, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: U_init_obs_covar
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: U_add_obs_error
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_enkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_add_obs_error, &
         U_init_obs_covar, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_enkf

   SUBROUTINE c__PDAF_assimilate_estkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_init_obsvar, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 HV
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_estkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_init_obsvar, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_estkf

   SUBROUTINE c__PDAF_assimilate_etkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_init_obsvar, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 HV
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_etkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_init_obsvar, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_etkf

   SUBROUTINE c__PDAF_assimilate_hyb3dvar_estkf(U_collect_state, U_distribute_state, &
                                                U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                                U_init_obsvar, U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
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
      call PDAF_assimilate_hyb3dvar_estkf(U_collect_state, U_distribute_state, &
                                                U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                                U_init_obsvar, U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_assimilate_hyb3dvar_estkf

   SUBROUTINE c__PDAF_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
                                                 U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                 U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                                 U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
                                                 U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
                                                 U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
                                                 U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      call PDAF_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
                                                 U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
                                                 U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
                                                 U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
                                                 U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
                                                 U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
                                                 U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAF_assimilate_lenkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_localize, &
         U_add_obs_error, U_init_obs_covar, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: U_init_obs_covar
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: U_localize
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: U_add_obs_error
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_lenkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_localize, &
         U_add_obs_error, U_init_obs_covar, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_lenkf

   SUBROUTINE c__PDAF_assimilate_lestkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_lestkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_lestkf

   SUBROUTINE c__PDAF_assimilate_letkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_letkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_letkf

   SUBROUTINE c__PDAF_assimilate_lnetf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
         U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_lnetf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
         U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_lnetf

   SUBROUTINE c__PDAF_assimilate_lknetf(U_collect_state, U_distribute_state, &
                                        U_init_dim_obs, U_obs_op, U_init_obs, &
                                        U_init_obs_l, U_prepoststep, &
                                        U_prodRinvA_l, U_prodRinvA_hyb_l, &
                                        U_init_n_domains_p, U_init_dim_l, &
                                        U_init_dim_obs_l, &
                                        U_g2l_state, U_l2g_state, U_g2l_obs, &
                                        U_init_obsvar, U_init_obsvar_l, &
                                        U_likelihood_l, U_likelihood_hyb_l, &
                                        U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      call PDAF_assimilate_lknetf(U_collect_state, U_distribute_state, &
                                  U_init_dim_obs, U_obs_op, U_init_obs, &
                                  U_init_obs_l, U_prepoststep, &
                                  U_prodRinvA_l, U_prodRinvA_hyb_l, &
                                  U_init_n_domains_p, U_init_dim_l, &
                                  U_init_dim_obs_l, &
                                  U_g2l_state, U_l2g_state, U_g2l_obs, &
                                  U_init_obsvar, U_init_obsvar_l, &
                                  U_likelihood_l, U_likelihood_hyb_l, &
                                  U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_lknetf

   SUBROUTINE c__PDAF_assimilate_lseik(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_lseik(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
         U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_lseik

   SUBROUTINE c__PDAF_assimilate_netf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
         U_likelihood, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: U_likelihood
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_netf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
         U_likelihood, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_netf

   SUBROUTINE c__PDAF_assimilate_pf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
         U_likelihood, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: U_likelihood
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_pf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
         U_likelihood, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_pf

   SUBROUTINE c__PDAF_assimilate_seek(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 HV
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_seek(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_seek

   SUBROUTINE c__PDAF_assimilate_seik(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 HV
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_assimilate_seik(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_seik

   SUBROUTINE c__PDAF_assimilate_prepost(U_collect_state, U_distribute_state, &
                                         U_prepoststep, U_next_observation, flag) &
                                         bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      call PDAF_assimilate_prepost(U_collect_state, U_distribute_state, &
                                 U_prepoststep, U_next_observation, flag)
   END SUBROUTINE c__PDAF_assimilate_prepost

   SUBROUTINE c__PDAF_deallocate() bind(c)
      CALL PDAF_deallocate()
   END SUBROUTINE c__PDAF_deallocate

   SUBROUTINE c__PDAF_diag_effsample(dim_sample, weights, effSample) bind(c)
      ! Sample size
      INTEGER(c_int), INTENT(in) :: dim_sample
      ! weights of the samples
      REAL(c_double), INTENT(in) :: weights(dim_sample)
      ! effecfive sample size
      REAL(c_double), INTENT(out) :: effSample

      CALL PDAF_diag_effsample(dim_sample, weights, effSample)
   END SUBROUTINE c__PDAF_diag_effsample

   SUBROUTINE c__PDAF_diag_ensstats(dim, dim_ens, element, &
         state, ens, skewness, kurtosis, status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Index of state vector/ensemble element to be used.
      ! If element=0, mean values over all elements are computed
      INTEGER(c_int), INTENT(in) :: element
      ! State vector (typically ensemble mean)
      REAL(c_double), INTENT(in) :: state(dim)
      ! State ensemble
      REAL(c_double), INTENT(in) :: ens(dim, dim_ens)
      ! Skewness of ensemble
      REAL(c_double), INTENT(out) :: skewness
      ! Kurtosis of ensemble
      REAL(c_double), INTENT(out) :: kurtosis
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status

      CALL PDAF_diag_ensstats(dim, dim_ens, element, &
         state, ens, skewness, kurtosis, status)
   END SUBROUTINE c__PDAF_diag_ensstats

   SUBROUTINE c__PDAF_diag_histogram(ncall, dim, dim_ens, element, &
         state, ens, hist, delta, status) bind(c)
      ! The number of calls used to increment the histogram and
      ! is needed to compute the delta-measure that
      ! describes the deviation from the ideal histogram.
      INTEGER(c_int), INTENT(in) :: ncall
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Element of vector used for histogram. If element=0, all elements are used
      INTEGER(c_int), INTENT(in) :: element
      ! Assumed truth
      REAL(c_double), INTENT(in) :: state(dim)
      ! Ensemble
      REAL(c_double), INTENT(in) :: ens(dim, dim_ens)
      ! Histogram about the state
      INTEGER(c_int), INTENT(inout) :: hist(dim_ens+1)
      ! deviation measure from flat histogram.
      ! It must be initialised to be 0
      REAL(c_double), INTENT(out) :: delta
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status

      CALL PDAF_diag_histogram(ncall, dim, dim_ens, element, &
         state, ens, hist, delta, status)
   END SUBROUTINE c__PDAF_diag_histogram

   SUBROUTINE c__PDAF_eofcovar(dim_state, nstates, nfields, dim_fields, offsets, &
         remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose, status) bind(c)
      ! Dimension of state vector
      INTEGER(c_int), INTENT(in) :: dim_state
      ! Number of state vectors, typically from different time steps
      INTEGER(c_int), INTENT(in) :: nstates
      ! Number of fields in state vector
      ! For example, if the state vector contains temperature and humidity,
      ! nfields=2, only used when `do_mv = 1`
      INTEGER(c_int), INTENT(in) :: nfields
      ! Size of each field, only used when `do_mv = 1`.
      ! Each field could be 2D or 3D so can have different sizes
      INTEGER(c_int), INTENT(in) :: dim_fields(nfields)
      ! Start position of each field, only used when `do_mv = 1`
      ! It starts from 1.
      INTEGER(c_int), INTENT(in) :: offsets(nfields)
      ! 1: Do multivariate scaling; 0: no scaling
      ! nfields, dim_fields and offsets are only used if do_mv=1
      INTEGER(c_int), INTENT(in) :: do_mv
      ! 1: subtract mean state from states (average over nstates dimension)
      ! before computing EOFs; 0: don't remove
      INTEGER(c_int), INTENT(in) :: remove_mstate
      ! State perturbations or an ensemble of state vectors
      REAL(c_double), INTENT(inout) :: states(dim_state, nstates)
      ! Standard deviation of field variability.
      ! Without multivariate scaling (do_mv=0), it is stddev = 1.0
      REAL(c_double), INTENT(out) :: stddev(nfields)
      ! Singular values divided by sqrt(nstates-1)
      REAL(c_double), INTENT(out) :: svals(nstates)
      ! Singular vectors
      REAL(c_double), INTENT(out) :: svec(dim_state, nstates)
      ! Mean state (only changed if remove_mstate=1)
      REAL(c_double), INTENT(inout) :: meanstate(dim_state)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status

      CALL PDAF_eofcovar(dim_state, nstates, nfields, dim_fields, offsets, &
         remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose, status)
   END SUBROUTINE c__PDAF_eofcovar

   SUBROUTINE c__PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_f

      CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)
   END SUBROUTINE c__PDAF_gather_dim_obs_f

   SUBROUTINE c__PDAF_gather_obs_f(obs_p, dimobs_p, obs_f, dimobs_f, status) bind(c)
      ! dimensions of PE local obs
      INTEGER(c_int), intent(in) :: dimobs_p
      ! dimension of full gathered obs
      INTEGER(c_int), intent(in) :: dimobs_f
      ! PE-local vector
      REAL(c_double), INTENT(in) :: obs_p(dimobs_p)
      ! Full gathered vector
      REAL(c_double), INTENT(out) :: obs_f(dimobs_f)
      ! Status flag:
      ! (0) no error;
      ! (1) when PDAF_gather_dim_obs_f not executed before
      INTEGER(c_int), INTENT(out) :: status

      CALL PDAF_gather_obs_f(obs_p, obs_f, status)
   END SUBROUTINE c__PDAF_gather_obs_f

   SUBROUTINE c__PDAF_gather_obs_f2(coords_p, dimobs_p, coords_f, dimobs_f, nrows, status) bind(c)
      ! dimensions of PE local obs
      INTEGER(c_int), intent(in) :: dimobs_p
      ! dimension of full gathered obs
      INTEGER(c_int), intent(in) :: dimobs_f
      ! Number of rows in array
      INTEGER(c_int), INTENT(in) :: nrows
      ! PE-local array
      REAL(c_double), INTENT(in) :: coords_p(nrows, dimobs_p)
      ! Full gathered array
      REAL(c_double), INTENT(out) :: coords_f(nrows, dimobs_f)
      ! Status flag:
      ! (0) no error;
      ! (1) when PDAF_gather dim_obs_f not executed before
      INTEGER(c_int), INTENT(out) :: status

      CALL PDAF_gather_obs_f2(coords_p, coords_f, nrows, status)
   END SUBROUTINE c__PDAF_gather_obs_f2

   SUBROUTINE c__PDAF_generate_obs(U_collect_state, U_distribute_state, &
         U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_init_obserr_f, U_prepoststep, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: U_get_obs_f
      ! Initialize vector of observation error standard deviations
      procedure(c__init_obserr_f_pdaf) :: U_init_obserr_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAF_generate_obs(U_collect_state, U_distribute_state, &
         U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_init_obserr_f, U_prepoststep, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAF_generate_obs

   SUBROUTINE c__PDAF_get_assim_flag(did_assim) bind(c)
      ! Flag: (1) for assimilation; (0) else
      INTEGER(c_int),INTENT(out) :: did_assim

      CALL PDAF_get_assim_flag(did_assim)
   END SUBROUTINE c__PDAF_get_assim_flag

   SUBROUTINE c__PDAF_get_ensstats(dims, c_skew_ptr, c_kurt_ptr, status) bind(c)
      ! dimension of pointer
      INTEGER(c_int), intent(out) :: dims(1)
      ! Skewness array
      type(c_ptr), INTENT(out) :: c_skew_ptr
      ! kurtosis array
      type(c_ptr), INTENT(out) :: c_kurt_ptr
      ! Status flag
      INTEGER(c_int), INTENT(out)       :: status

      REAL, POINTER :: skew_ptr(:)
      REAL, POINTER :: kurt_ptr(:)

      call PDAF_get_ensstats(skew_ptr, kurt_ptr, status)
      dims = shape(skew_ptr)
      c_skew_ptr = c_loc(skew_ptr(1))
      c_kurt_ptr = c_loc(kurt_ptr(1))
   END SUBROUTINE c__PDAF_get_ensstats

   SUBROUTINE c__PDAF_get_localfilter(lfilter) bind(c)
      ! Whether the filter is domain-localized (1) or not (0)
      INTEGER(c_int), INTENT(out) :: lfilter

      CALL PDAF_get_localfilter(lfilter)
   END SUBROUTINE c__PDAF_get_localfilter

   SUBROUTINE c__PDAF_get_memberid(memberid) bind(c)
      ! Index in the local ensemble
      INTEGER(c_int),INTENT(inout) :: memberid

      CALL PDAF_get_memberid(memberid)
   END SUBROUTINE c__PDAF_get_memberid

   SUBROUTINE c__PDAF_get_obsmemberid(memberid) bind(c)
      ! Index in the local observed ensemble
      INTEGER(c_int),INTENT(inout) :: memberid

      CALL PDAF_get_obsmemberid(memberid)
   END SUBROUTINE c__PDAF_get_obsmemberid

   SUBROUTINE c__PDAF_get_smootherens(c_sens_point, maxlag, dims, status) bind(c)
      ! A smoother array
      type(c_ptr), intent(out) :: c_sens_point
      ! Number of past timesteps processed in sens
      INTEGER(c_int), INTENT(out) :: maxlag
      ! dimension of pointer
      INTEGER(c_int), intent(out) :: dims(3)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status

      REAL, POINTER :: sens_point(:, :, :)

      CALL PDAF_get_smootherens(sens_point, maxlag, status)

      dims = shape(sens_point)
      c_sens_point = c_loc(sens_point(1, 1, 1))
   END SUBROUTINE c__PDAF_get_smootherens

   SUBROUTINE c__PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
         U_prepoststep, flag) bind(c)
      ! Flag and number of time steps
      INTEGER(c_int), INTENT(inout) :: steps
      ! current model time
      REAL(c_double), INTENT(out) :: time
      ! Whether to exit from forecasts
      INTEGER(c_int), INTENT(inout) :: doexit
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
         U_prepoststep, flag)
   END SUBROUTINE c__PDAF_get_state

   SUBROUTINE c__PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
         param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, in_filterpe, U_init_ens, in_screen, &
         flag) bind(c)
      ! Type of filter
      INTEGER(c_int), INTENT(in) :: filtertype
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Initial time step of assimilation
      INTEGER(c_int), INTENT(in) :: stepnull
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Integer parameter array
      INTEGER(c_int), INTENT(inout) :: param_int(dim_pint)
      ! Number of real parameter
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Real parameter array
      REAL(c_double), INTENT(inout) :: param_real(dim_preal)
      ! Model communicator
      INTEGER(c_int), INTENT(in) :: COMM_model
      ! Coupling communicator
      INTEGER(c_int), INTENT(in) :: COMM_couple
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: COMM_filter
      ! Id of my ensemble task
      INTEGER(c_int), INTENT(in) :: task_id
      ! Number of parallel model tasks
      INTEGER(c_int), INTENT(in) :: n_modeltasks
      ! Is my PE a filter-PE?
      LOGICAL(c_bool), INTENT(in) :: in_filterpe
      ! Control screen output:
      INTEGER(c_int), INTENT(in) :: in_screen
      ! Status flag, 0: no error, error codes:
      INTEGER(c_int), INTENT(out):: flag
      ! User-supplied routine for ensemble initialization
      procedure(c__init_ens_pdaf) :: U_init_ens
      logical :: filterpe
      filterpe = in_filterpe

      CALL PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
         param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, filterpe, U_init_ens, in_screen, &
         flag)
   END SUBROUTINE c__PDAF_init

   SUBROUTINE c__PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
         nrows, ncols, A, var_obs, weight, verbose) bind(c)
      ! Type of weight function
      INTEGER(c_int), INTENT(in) :: wtype
      ! Type of regulated weighting
      INTEGER(c_int), INTENT(in) :: rtype
      ! Cut-off radius (check https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE)
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius (check https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE)
      REAL(c_double), INTENT(in) :: sradius
      ! Distance to observation
      REAL(c_double), INTENT(in) :: distance
      ! Number of rows in matrix A
      INTEGER(c_int), INTENT(in) :: nrows
      ! Number of columns in matrix A
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix
      REAL(c_double), INTENT(in) :: A(nrows, ncols)
      ! Observation variance
      REAL(c_double), INTENT(in) :: var_obs
      ! Weights
      REAL(c_double), INTENT(out) :: weight
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose

      CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
         nrows, ncols, A, var_obs, weight, verbose)
   END SUBROUTINE c__PDAF_local_weight

   SUBROUTINE c__PDAF_print_info(printtype) bind(c)
      ! - X=1: Basic timers
      ! - X=3: Timers showing the time spent int he different call-back routines (this variant was added with PDAF 1.15)
      ! - X=4: More detailed timers about parts of the filter algorithm (before PDAF 1.15, this was timer level 3)
      ! - X=5: Very detailed timers about various operations in the filter algorithm (before PDAF 1.15, this was timer level 4)
      INTEGER(c_int), INTENT(in) :: printtype

      CALL PDAF_print_info(printtype)
   END SUBROUTINE c__PDAF_print_info

   SUBROUTINE c__PDAF_put_state_3dvar(U_collect_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: U_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: U_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: U_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: U_obs_op_adj
      call PDAF_put_state_3dvar(U_collect_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
         U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
         U_prepoststep, outflag)
   END SUBROUTINE c__PDAF_put_state_3dvar

   SUBROUTINE c__PDAF_put_state_en3dvar_estkf(U_collect_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_obsvar, U_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
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
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      call PDAF_put_state_en3dvar_estkf(U_collect_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
         U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
         U_init_obsvar, U_prepoststep, outflag)
   END SUBROUTINE c__PDAF_put_state_en3dvar_estkf

   SUBROUTINE c__PDAF_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      call PDAF_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prodRinvA, &
         U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
         U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
         U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
         U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_prepoststep, outflag)
   END SUBROUTINE c__PDAF_put_state_en3dvar_lestkf

   SUBROUTINE c__PDAF_put_state_enkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_init_obs, U_prepoststep, U_add_obs_err, U_init_obs_covar, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: U_add_obs_err
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: U_init_obs_covar

      CALL PDAF_put_state_enkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_init_obs, U_prepoststep, U_add_obs_err, U_init_obs_covar, flag)
   END SUBROUTINE c__PDAF_put_state_enkf

   SUBROUTINE c__PDAF_put_state_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA

      CALL PDAF_put_state_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
   END SUBROUTINE c__PDAF_put_state_estkf

   SUBROUTINE c__PDAF_put_state_etkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA

      CALL PDAF_put_state_etkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
   END SUBROUTINE c__PDAF_put_state_etkf

   SUBROUTINE c__PDAF_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
         U_get_obs_f, U_init_obserr_f, U_prepoststep, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: U_get_obs_f
      ! Initialize vector of observation errors
      procedure(c__init_obserr_f_pdaf) :: U_init_obserr_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
         U_get_obs_f, U_init_obserr_f, U_prepoststep, flag)
   END SUBROUTINE c__PDAF_put_state_generate_obs

   SUBROUTINE c__PDAF_put_state_hyb3dvar_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
      U_init_obs, U_prodRinvA, &
      U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
      U_init_obsvar, U_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
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
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: U_init_obsvar
      call PDAF_put_state_hyb3dvar_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prodRinvA, &
         U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
         U_init_obsvar, U_prepoststep, outflag)
   END SUBROUTINE c__PDAF_put_state_hyb3dvar_estkf

   SUBROUTINE c__PDAF_put_state_hyb3dvar_lestkf(U_collect_state, &
      U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
      U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
      U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
      U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
      U_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      call PDAF_put_state_hyb3dvar_lestkf(U_collect_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
         U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
         U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
         U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
         U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
         U_prepoststep, outflag)
   END SUBROUTINE c__PDAF_put_state_hyb3dvar_lestkf

   SUBROUTINE c__PDAF_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_init_obs, U_prepoststep, U_localize, U_add_obs_err, U_init_obs_covar, &
         flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: U_localize
      ! Add obs error covariance R to HPH in EnKF
      procedure(c__add_obs_err_pdaf) :: U_add_obs_err
      ! Initialize obs. error cov. matrix R in EnKF
      procedure(c__init_obs_covar_pdaf) :: U_init_obs_covar

      CALL PDAF_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_init_obs, U_prepoststep, U_localize, U_add_obs_err, U_init_obs_covar, &
         flag)
   END SUBROUTINE c__PDAF_put_state_lenkf

   SUBROUTINE c__PDAF_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag)
   END SUBROUTINE c__PDAF_put_state_lestkf

   SUBROUTINE c__PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag)
   END SUBROUTINE c__PDAF_put_state_letkf

   SUBROUTINE c__PDAF_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: U_likelihood_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         outflag)
   END SUBROUTINE c__PDAF_put_state_lnetf

   SUBROUTINE c__PDAF_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
                                       U_init_obs, U_init_obs_l, U_prepoststep, &
                                       U_prodRinvA_l, U_prodRinvA_hyb_l, &
                                       U_init_n_domains_p, &
                                       U_init_dim_l, U_init_dim_obs_l, &
                                       U_g2l_state, U_l2g_state, U_g2l_obs, &
                                       U_init_obsvar, U_init_obsvar_l, &
                                       U_likelihood_l, &
                                       U_likelihood_hyb_l, outflag) bind(c)
            ! Status flag
            INTEGER(c_int), INTENT(out) :: outflag
            ! Routine to collect a state vector
            procedure(c__collect_state_pdaf) :: U_collect_state
            ! Provide number of local analysis domains
            procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
            ! Init state dimension for local ana. domain
            procedure(c__init_dim_l_pdaf) :: U_init_dim_l
            ! Observation operator
            procedure(c__obs_op_pdaf) :: U_obs_op
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
            ! Get state on local ana. domain from full state
            procedure(c__g2l_state_pdaf) :: U_g2l_state
            ! Init full state from state on local analysis domain
            procedure(c__l2g_state_pdaf) :: U_l2g_state
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
            ! User supplied pre/poststep routine
            procedure(c__prepoststep_pdaf) :: U_prepoststep
            call PDAF_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
                                       U_init_obs, U_init_obs_l, U_prepoststep, &
                                       U_prodRinvA_l, U_prodRinvA_hyb_l, &
                                       U_init_n_domains_p, &
                                       U_init_dim_l, U_init_dim_obs_l, &
                                       U_g2l_state, U_l2g_state, U_g2l_obs, &
                                       U_init_obsvar, U_init_obsvar_l, &
                                       U_likelihood_l, &
                                       U_likelihood_hyb_l, outflag)
   END SUBROUTINE c__PDAF_put_state_lknetf

   SUBROUTINE c__PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: U_g2l_obs
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodRinvA_l_pdaf) :: U_prodRinvA_l
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
         U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
         U_init_obsvar, U_init_obsvar_l, flag)
   END SUBROUTINE c__PDAF_put_state_lseik

   SUBROUTINE c__PDAF_put_state_netf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_likelihood, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: U_likelihood

      CALL PDAF_put_state_netf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_likelihood, flag)
   END SUBROUTINE c__PDAF_put_state_netf

   SUBROUTINE c__PDAF_put_state_pf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_likelihood, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: U_likelihood

      CALL PDAF_put_state_pf(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_likelihood, flag)
   END SUBROUTINE c__PDAF_put_state_pf

   SUBROUTINE c__PDAF_put_state_prepost(U_collect_state, U_prepoststep, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAF_put_state_prepost(U_collect_state, U_prepoststep, flag)
   END SUBROUTINE c__PDAF_put_state_prepost

   SUBROUTINE c__PDAF_put_state_seek(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 HV
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA

      CALL PDAF_put_state_seek(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, flag)
   END SUBROUTINE c__PDAF_put_state_seek

   SUBROUTINE c__PDAF_put_state_seik(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide product R^-1 A
      procedure(c__prodRinvA_pdaf) :: U_prodRinvA

      CALL PDAF_put_state_seik(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
   END SUBROUTINE c__PDAF_put_state_seik

   subroutine c__PDAF_reset_forget(forget_in) bind(c)
      ! New value of forgetting factor
      REAL(c_double), INTENT(in) :: forget_in
      call PDAF_reset_forget(forget_in)
   END subroutine c__PDAF_reset_forget

   SUBROUTINE c__PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens, verbose, flag) bind(c)
      ! Size of state vector
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Array of EOF modes
      REAL(c_double), INTENT(inout) :: modes(dim, dim_ens-1)
      ! Vector of singular values
      REAL(c_double), INTENT(in) :: svals(dim_ens-1)
      ! PE-local model mean state
      REAL(c_double), INTENT(inout) :: state(dim)
      ! State ensemble
      REAL(c_double), INTENT(out) :: ens(dim, dim_ens)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      CALL PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens, verbose, flag)
   END SUBROUTINE c__PDAF_SampleEns

   SUBROUTINE c__PDAF_set_debug_flag(debugval) bind(c)
      ! Value of debugging flag; print debug information for >0
      INTEGER(c_int), INTENT(in)        :: debugval
      call PDAF_set_debug_flag(debugval)
   END SUBROUTINE c__PDAF_set_debug_flag

   SUBROUTINE c__PDAF_set_ens_pointer(c_ens_point, dims, status) bind(c)
      ! Pointer to smoother array
      type(c_ptr), intent(out) :: c_ens_point
      ! dimension of the pointer
      integer(c_int), intent(out) :: dims(2)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status
      ! Pointer to smoother array
      REAL, POINTER :: ens_point(:,:)
      CALL PDAF_set_ens_pointer(ens_point, status)
      c_ens_point = c_loc(ens_point(1,1))
      dims = shape(ens_point)
   END SUBROUTINE c__PDAF_set_ens_pointer

   SUBROUTINE c__PDAF_set_smootherens(c_sens_point, maxlag, dims, status) bind(c)
      ! Pointer to smoother array
      type(c_ptr), INTENT(out) :: c_sens_point
      ! Number of past timesteps processed in sens
      INTEGER(c_int), INTENT(in) :: maxlag
      ! dimension of the pointer
      integer(c_int), intent(out) :: dims(3)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status
      ! Pointer to smoother array
      REAL, POINTER :: sens_point(:,:,:)

      CALL PDAF_set_smootherens(sens_point, maxlag, status)
      c_sens_point = c_loc(sens_point(1,1,1))
      dims = shape(sens_point)
   END SUBROUTINE c__PDAF_set_smootherens

   SUBROUTINE c__PDAF_seik_TtimesA(rank, dim_col, A, B) bind(c)
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim_col
      ! Input matrix
      REAL(c_double), INTENT(in) :: A(rank, dim_col)
      ! Output matrix (TA)
      REAL(c_double), INTENT(out) :: B(rank+1, dim_col)

      CALL PDAF_seik_TtimesA(rank, dim_col, A, B)
   END SUBROUTINE c__PDAF_seik_TtimesA

   SUBROUTINE c__PDAF_etkf_Tleft(dim_ens, dim, A) bind(c)
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim
      ! Input/output matrix
      REAL(c_double), INTENT(inout) :: A(dim_ens, dim)

      CALL PDAF_etkf_Tleft(dim_ens, dim, A)
   END SUBROUTINE c__PDAF_etkf_Tleft

   SUBROUTINE c__PDAF_estkf_OmegaA(rank, dim_col, A, B) bind(c)
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim_col
      ! Input matrix
      REAL(c_double), INTENT(in) :: A(rank, dim_col)
      ! Output matrix (TA)
      REAL(c_double), INTENT(out) :: B(rank+1, dim_col)

      CALL PDAF_estkf_OmegaA(rank, dim_col, A, B)
   END SUBROUTINE c__PDAF_estkf_OmegaA

   SUBROUTINE c__PDAF_enkf_omega(seed, r, dim_ens, omega, norm, &
         otype, screen) bind(c)
      ! Seed for random number generation
      INTEGER(c_int), INTENT(in) :: seed(4)
      ! Approximated rank of covar matrix
      INTEGER(c_int), INTENT(in) :: r
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Random matrix
      REAL(c_double), INTENT(inout) :: omega(dim_ens,r)
      ! Norm for ensemble transformation
      REAL(c_double), INTENT(inout) :: norm
      ! Type of Omega:
      ! (1) Simple Gaussian random matrix
      ! (2) Columns of unit norm
      ! (3) Columns of norm dim_ens^(-1/2)
      ! (4) Projection orthogonal (1,..,1)^T
      ! (6) Combination of 2 and 4
      ! (7) Combination of 3 and 4
      ! (8) Rows of sum 0 and variance 1
      INTEGER(c_int), INTENT(in) :: otype
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      CALL PDAF_enkf_omega(seed, r, dim_ens, omega, norm, &
         otype, screen)
   END SUBROUTINE c__PDAF_enkf_omega

   SUBROUTINE c__PDAF_seik_omega(rank, omega, omegatype, screen) bind(c)
      ! Approximated rank of covar matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Matrix Omega
      REAL(c_double), INTENT(inout) :: omega(rank+1, rank)
      ! Select type of Omega:
      !   (1) generated from random vectors
      !   (0) generated from deterministic vectors (Householder)
      INTEGER(c_int), INTENT(in) :: omegatype
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      CALL PDAF_seik_omega(rank, omega, omegatype, screen)
   END SUBROUTINE c__PDAF_seik_omega

   SUBROUTINE c__PDAF_incremental(steps, U_dist_stateinc) bind(c)
      ! Time steps over which increment is distributed
      INTEGER(c_int), INTENT(in) :: steps
      ! Add state increment during integration
      procedure(c__dist_stateinc_pdaf) :: U_dist_stateinc

      CALL PDAF_incremental(steps, U_dist_stateinc)
   END SUBROUTINE c__PDAF_incremental

   SUBROUTINE c__PDAF_add_increment(dim_p, state_p) bind(c)
      ! State dimension
      INTEGER(c_int),INTENT(in) :: dim_p
      ! State vector
      REAL(c_double),INTENT(inout) :: state_p(dim_p)

      CALL PDAF_add_increment(dim_p, state_p)
   END SUBROUTINE c__PDAF_add_increment

   ! In documentation but not in the PDAF_interface module
   SUBROUTINE c__PDAF_local_weights(wtype, cradius, sradius, dim, distance, &
      weight, verbose) bind(c)
      ! Type of weight function
      ! (0): unit weight (=1 up to distance=cradius)
      ! (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
      ! (2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
      INTEGER(c_int), INTENT(in) :: wtype
      ! Parameter for cut-off
      REAL(c_double), INTENT(in)    :: cradius
      ! Support radius
      REAL(c_double), INTENT(in)    :: sradius
      ! Size of distance and weight arrays
      INTEGER(c_int), INTENT(in) :: dim
      ! Array holding distances
      REAL(c_double), INTENT(in)    :: distance(dim)
      ! Array for weights
      REAL(c_double), INTENT(out)   :: weight(dim)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      call  PDAF_local_weights(wtype, cradius, sradius, dim, distance, &
                               weight, verbose)
   END SUBROUTINE c__PDAF_local_weights

   SUBROUTINE c__PDAF_diag_CRPS(dim, dim_ens, element, oens, obs, &
      CRPS, reli, resol, uncert, status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! ID of element to be used. If element=0, mean values over all elements are computed
      INTEGER(c_int), INTENT(in) :: element
      ! State ensemble
      REAL(c_double), INTENT(in)    :: oens(dim, dim_ens)
      ! State ensemble
      REAL(c_double), INTENT(in)    :: obs(dim)
      ! CRPS
      REAL(c_double), INTENT(out)   :: CRPS
      ! Reliability
      REAL(c_double), INTENT(out)   :: reli
      ! resolution
      REAL(c_double), INTENT(out)   :: resol
      ! uncertainty
      REAL(c_double), INTENT(out)   :: uncert
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status
      call PDAF_diag_CRPS(dim, dim_ens, element, oens, obs, &
         CRPS, reli, resol, uncert, status)
   END SUBROUTINE c__PDAF_diag_CRPS

   SUBROUTINE c__PDAF_force_analysis() bind(c)

      call PDAF_force_analysis()
   END SUBROUTINE c__PDAF_force_analysis

   SUBROUTINE c__PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p, coords_f, &
        nrows, status) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Number of rows in array
      INTEGER(c_int), INTENT(in) :: nrows
      ! PE-local array
      REAL(c_double), INTENT(in)  :: coords_p(nrows, dim_obs_p)
      ! Full gathered array
      REAL(c_double), INTENT(out) :: coords_f(nrows, dim_obs_f)
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status
      call PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p, coords_f, &
        nrows, status)
   END SUBROUTINE c__PDAF_gather_obs_f2_flex

   SUBROUTINE c__PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, obs_f, status) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! PE-local vector
      REAL(c_double), INTENT(in)  :: obs_p(dim_obs_p)
      ! Full gathered vector
      REAL(c_double), INTENT(out) :: obs_f(dim_obs_f)
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status
      call PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, obs_f, status)
   END SUBROUTINE c__PDAF_gather_obs_f_flex

   SUBROUTINE c__PDAF_prepost(U_collect_state, U_distribute_state, &
      U_prepoststep, U_next_observation, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Routine to provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      call PDAF_prepost(U_collect_state, U_distribute_state, &
         U_prepoststep, U_next_observation, outflag)
   END SUBROUTINE c__PDAF_prepost

   SUBROUTINE c__PDAF_set_memberid(memberid) bind(c)
      ! Index in the local ensemble
      INTEGER(c_int),INTENT(inout) :: memberid
      call PDAF_set_memberid(memberid)
   END SUBROUTINE c__PDAF_set_memberid

   SUBROUTINE c__PDAF_set_comm_pdaf(in_COMM_pdaf) bind(c)
      ! MPI communicator for PDAF
     INTEGER(c_int),INTENT(in) :: in_COMM_pdaf

     call PDAF_set_comm_pdaf(in_COMM_pdaf)
   END SUBROUTINE c__PDAF_set_comm_pdaf

   ! added interface in PDAF V2.2.1
   SUBROUTINE c__PDAF_set_offline_mode(screen) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in)        :: screen
      call PDAF_set_offline_mode(screen)
   END SUBROUTINE c__PDAF_set_offline_mode

   ! in PDAF_analysis_utils.F90 in V2.2.1
   SUBROUTINE c__PDAF_print_domain_stats(n_domains_p) bind(c)
      ! Number of PE-local analysis domains
      INTEGER(c_int), INTENT(in) :: n_domains_p
      call PDAF_print_domain_stats(n_domains_p)
   END SUBROUTINE c__PDAF_print_domain_stats

   SUBROUTINE c__PDAF_init_local_obsstats() bind(c)
      call PDAF_init_local_obsstats()
   END SUBROUTINE c__PDAF_init_local_obsstats

   SUBROUTINE c__PDAF_incr_local_obsstats(dim_obs_l) bind(c)
      ! Number of locally assimilated observations
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      call PDAF_incr_local_obsstats(dim_obs_l)
   END SUBROUTINE c__PDAF_incr_local_obsstats

   SUBROUTINE c__PDAF_print_local_obsstats(screen, n_domains_with_obs) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! number of domains with observations
      INTEGER(c_int), OPTIONAL, INTENT(out) :: n_domains_with_obs
      call PDAF_print_local_obsstats(screen, n_domains_with_obs)
   END SUBROUTINE c__PDAF_print_local_obsstats

   SUBROUTINE c__PDAF_omit_obs_omi(dim_p, dim_obs_p, dim_ens, state_p, ens_p, &
                                   obs_p, U_init_obs, U_obs_op, compute_mean, screen) bind(c)
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: U_init_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! on exit: PE-local forecast mean state
      REAL(c_double), INTENT(inout) :: state_p(dim_p)
      ! PE-local state ensemble
      REAL(c_double), INTENT(in) :: ens_p(dim_p, dim_ens)
      ! PE-local observation vector
      REAL(c_double), INTENT(inout) :: obs_p(dim_obs_p)
      ! (1) compute mean; (0) state_p holds mean
      INTEGER(c_int), INTENT(in) :: compute_mean
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      call PDAF_omit_obs_omi(dim_p, dim_obs_p, dim_ens, state_p, ens_p, &
                              obs_p, U_init_obs, U_obs_op, compute_mean, screen)
   END SUBROUTINE c__PDAF_omit_obs_omi

   SUBROUTINE c__PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, &
                                      CRPS, reli, resol, uncert, status)!
      ! PE-local state dimension
      INTEGER, INTENT(in) :: dim
      ! Ensemble size
      INTEGER, INTENT(in) :: dim_ens
      ! ID of element to be used
      ! If element=0, mean values over all elements are computed
      INTEGER, INTENT(in) :: element
      ! State ensemble
      REAL, INTENT(in)    :: oens(dim, dim_ens)
      ! State ensemble
      REAL, INTENT(in)    :: obs(dim)
      ! CRPS
      REAL, INTENT(out)   :: CRPS
      ! Reliability
      REAL, INTENT(out)   :: reli
      ! resolution
      REAL, INTENT(out)   :: resol
      ! uncertainty
      REAL, INTENT(out)   :: uncert
      ! Status flag (0=success)
      INTEGER, INTENT(out) :: status

      call PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, &
                                CRPS, reli, resol, uncert, status)
   END SUBROUTINE c__PDAF_diag_CRPS_nompi
END MODULE PDAF_c_binding
