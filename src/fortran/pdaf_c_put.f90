MODULE pdaf_c_put
use PDAF
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAF_put_state_lestkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_lestkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
         u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
         u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag)

   END SUBROUTINE c__PDAF_put_state_lestkf

   SUBROUTINE c__PDAF_put_state_lenkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prepoststep, u_localize, u_add_obs_err,  &
      u_init_obs_covar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_lenkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_localize, u_add_obs_err,  &
         u_init_obs_covar, outflag)

   END SUBROUTINE c__PDAF_put_state_lenkf

   SUBROUTINE c__PDAF_put_state_lseik(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_lseik(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
         u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
         u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag)

   END SUBROUTINE c__PDAF_put_state_lseik

   SUBROUTINE c__PDAF_put_state_en3dvar_lestkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens,  &
      u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f, u_init_obs_f,  &
      u_init_obs_l, u_prodrinva_l, u_init_n_domains_p, u_init_dim_l,  &
      u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar,  &
      u_init_obsvar_l, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_en3dvar_lestkf(u_collect_state, u_init_dim_obs,  &
         u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens,  &
         u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f,  &
         u_init_obs_f, u_init_obs_l, u_prodrinva_l, u_init_n_domains_p,  &
         u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
         u_init_obsvar, u_init_obsvar_l, u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_en3dvar_lestkf

   SUBROUTINE c__PDAF_put_state_lknetf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_likelihood_l, u_likelihood_hyb_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_lknetf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
         u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_likelihood_l, u_likelihood_hyb_l, outflag)

   END SUBROUTINE c__PDAF_put_state_lknetf

   SUBROUTINE c__PDAF_put_state_hyb3dvar_estkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, u_init_obsvar, u_prepoststep,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_hyb3dvar_estkf(u_collect_state, u_init_dim_obs,  &
         u_obs_op, u_init_obs, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens,  &
         u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, u_init_obsvar,  &
         u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_hyb3dvar_estkf

   SUBROUTINE c__PDAF_put_state_lnetf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep, u_likelihood_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_lnetf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obs_l, u_prepoststep, u_likelihood_l,  &
         u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
         u_l2g_state, u_g2l_obs, outflag)

   END SUBROUTINE c__PDAF_put_state_lnetf

   SUBROUTINE c__PDAF_put_state_letkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_letkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obs_l, u_prepoststep, u_prodrinva_l,  &
         u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
         u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, outflag)

   END SUBROUTINE c__PDAF_put_state_letkf

   SUBROUTINE c__PDAF_put_state_en3dvar_estkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens,  &
      u_obs_op_lin, u_obs_op_adj, u_init_obsvar, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_en3dvar_estkf(u_collect_state, u_init_dim_obs,  &
         u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens,  &
         u_obs_op_lin, u_obs_op_adj, u_init_obsvar, u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_en3dvar_estkf

   SUBROUTINE c__PDAF_put_state_hyb3dvar_lestkf(u_collect_state,  &
      u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens,  &
      u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
      u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l, u_prodrinva_l,  &
      u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l, u_g2l_state,  &
      u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l, u_prepoststep,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_hyb3dvar_lestkf(u_collect_state, u_init_dim_obs,  &
         u_obs_op, u_init_obs, u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_cvt,  &
         u_cvt_adj, u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f,  &
         u_init_obs_f, u_init_obs_l, u_prodrinva_l, u_init_n_domains_p,  &
         u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
         u_init_obsvar, u_init_obsvar_l, u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_hyb3dvar_lestkf

   SUBROUTINE c__PDAF_put_state_netf(u_collect_state, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prepoststep, u_likelihood, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_netf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_likelihood, outflag)

   END SUBROUTINE c__PDAF_put_state_netf

   SUBROUTINE c__PDAF_put_state_etkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_etkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar, outflag)

   END SUBROUTINE c__PDAF_put_state_etkf

   SUBROUTINE c__PDAF_put_state_ensrf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_init_obsvars, u_localize_covar_serial,  &
      u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_ensrf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_init_obsvars, u_localize_covar_serial, u_prepoststep,  &
         outflag)

   END SUBROUTINE c__PDAF_put_state_ensrf

   SUBROUTINE c__PDAF_put_state_enkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prepoststep, u_add_obs_err, u_init_obs_covar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_enkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_add_obs_err, u_init_obs_covar, outflag)

   END SUBROUTINE c__PDAF_put_state_enkf

   SUBROUTINE c__PDAF_put_state_seik(u_collect_state, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_seik(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar, outflag)

   END SUBROUTINE c__PDAF_put_state_seik

   SUBROUTINE c__PDAF_put_state_generate_obs(u_collect_state, u_init_dim_obs_f,  &
      u_obs_op_f, u_init_obserr_f, u_get_obs_f, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_generate_obs(u_collect_state, u_init_dim_obs_f,  &
         u_obs_op_f, u_init_obserr_f, u_get_obs_f, u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_generate_obs

   SUBROUTINE c__PDAF_put_state_pf(u_collect_state, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prepoststep, u_likelihood, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_pf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_likelihood, outflag)

   END SUBROUTINE c__PDAF_put_state_pf

   SUBROUTINE c__PDAF_put_state_3dvar(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin,  &
      u_obs_op_adj, u_prepoststep, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_3dvar(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
         u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_3dvar

   SUBROUTINE c__PDAF_put_state_estkf(u_collect_state, u_init_dim_obs,  &
      u_obs_op, u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
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

      call PDAF_put_state_estkf(u_collect_state, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prepoststep, u_prodrinva, u_init_obsvar, outflag)

   END SUBROUTINE c__PDAF_put_state_estkf

   SUBROUTINE c__PDAF_put_state_prepost(u_collect_state, u_prepoststep,  &
      outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAF_put_state_prepost(u_collect_state, u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_put_state_prepost
END MODULE pdaf_c_put
