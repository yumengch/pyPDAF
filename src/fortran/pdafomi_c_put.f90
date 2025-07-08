MODULE pdafomi_c_put
use PDAF
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAFomi_put_state_local_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
      prodrinva_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdafomi_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product of inverse of R with matrix A
      procedure(c__prodrinva_l_pdafomi_pdaf) :: prodrinva_l_pdafomi
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf

      call PDAFomi_put_state_local_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
         prodrinva_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_local_nondiagR

   SUBROUTINE c__PDAFomi_put_state_global_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf_pdaf) :: prodrinva_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_global_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdaf, obs_op_pdaf, prodrinva_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_global_nondiagR

   SUBROUTINE c__PDAFomi_put_state_enkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, add_obs_error_pdafomi,  &
      init_obscovar_pdafomi, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Add observation error covariance matrix
      procedure(c__add_obs_error_pdafomi_pdaf) :: add_obs_error_pdafomi
      ! Initialize mean observation error variance
      procedure(c__init_obscovar_pdafomi_pdaf) :: init_obscovar_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_enkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, add_obs_error_pdafomi,  &
         init_obscovar_pdafomi, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_enkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_lenkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
      localize_covar_pdafomi, add_obs_error_pdafomi, init_obscovar_pdafomi,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdafomi_pdaf) :: localize_covar_pdafomi
      ! Provide product R^-1 A
      procedure(c__add_obs_error_pdafomi_pdaf) :: add_obs_error_pdafomi
      ! Initialize mean observation error variance
      procedure(c__init_obscovar_pdafomi_pdaf) :: init_obscovar_pdafomi

      call PDAFomi_put_state_lenkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
         localize_covar_pdafomi, add_obs_error_pdafomi, init_obscovar_pdafomi,  &
         outflag)

   END SUBROUTINE c__PDAFomi_put_state_lenkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_nonlin_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, likelihood_pdafomi,  &
      prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Compute likelihood
      procedure(c__likelihood_pdafomi_pdaf) :: likelihood_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_nonlin_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, likelihood_pdafomi,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_nonlin_nondiagR

   SUBROUTINE c__PDAFomi_put_state_lnetf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
      likelihood_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdafomi_pdaf) :: init_dim_obs_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdafomi_pdaf) :: likelihood_l_pdafomi
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf

      call PDAFomi_put_state_lnetf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
         likelihood_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_lnetf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_lknetf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
      prodrinva_l_pdafomi, prodrinva_hyb_l_pdafomi, likelihood_l_pdafomi,  &
      likelihood_hyb_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdafomi_pdaf) :: init_dim_obs_l_pdafomi
      ! Provide product R^-1 A on local analysis domain
      procedure(c__prodrinva_l_pdafomi_pdaf) :: prodrinva_l_pdafomi
      ! Product R^-1 A on local analysis domain with hybrid weight
      procedure(c__prodrinva_hyb_l_pdafomi_pdaf) :: prodrinva_hyb_l_pdafomi
      ! Compute likelihood and apply localization
      procedure(c__likelihood_l_pdafomi_pdaf) :: likelihood_l_pdafomi
      ! Compute likelihood and apply localization with tempering
      procedure(c__likelihood_hyb_l_pdafomi_pdaf) :: likelihood_hyb_l_pdafomi
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf

      call PDAFomi_put_state_lknetf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi,  &
         prodrinva_l_pdafomi, prodrinva_hyb_l_pdafomi, likelihood_l_pdafomi,  &
         likelihood_hyb_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_lknetf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf,  &
      obs_op_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf,  &
         obs_op_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_3dvar

   SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf(collect_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_en3dvar_estkf(collect_state_pdaf,  &
         init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
         obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf

   SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf(collect_state_pdaf,  &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf_pdaf) :: obs_op_f_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf_pdaf) :: init_dim_obs_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_en3dvar_lestkf(collect_state_pdaf,  &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
         obs_op_lin_pdaf, obs_op_adj_pdaf, init_n_domains_pdaf,  &
         init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf(collect_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf,  &
      cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! Apply ensemble control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint ensemble control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf_pdaf) :: obs_op_adj_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_hyb3dvar_estkf(collect_state_pdaf,  &
         init_dim_obs_pdaf, obs_op_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
         cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf(collect_state_pdaf,  &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
      cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf_pdaf) :: obs_op_f_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf_pdaf) :: init_dim_obs_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_hyb3dvar_lestkf(collect_state_pdaf,  &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf,  &
         cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf,  &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
         g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf

   SUBROUTINE c__PDAFomi_put_state_3dvar_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_pdaf,  &
      cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdafomi_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdafomi_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdafomi_pdaf) :: obs_op_adj_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_3dvar_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_pdaf,  &
         cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_3dvar_nondiagR

   SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
      prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdafomi_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdafomi_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdafomi_pdaf) :: obs_op_adj_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_en3dvar_estkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
         cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
      prodrinva_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf,  &
      init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdafomi_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdafomi_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdafomi_pdaf) :: obs_op_adj_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_l_pdafomi_pdaf) :: prodrinva_l_pdafomi
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdafomi_pdaf) :: init_dim_obs_l_pdafomi
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
         cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi,  &
         prodrinva_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf,  &
         prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi,  &
      obs_op_adj_pdafomi, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdafomi_pdaf) :: prodrinva_pdafomi
      ! Apply ensemble control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint ensemble control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdafomi_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdafomi_pdaf) :: obs_op_adj_pdafomi
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_hyb3dvar_estkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
         cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi,  &
         obs_op_adj_pdafomi, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf,  &
      init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
      cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi,  &
      obs_op_adj_pdafomi, prodrinva_l_pdafomi, init_n_domains_pdaf,  &
      init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf,  &
      prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_pdafomi_pdaf) :: init_dim_obs_pdafomi
      ! Full observation operator
      procedure(c__obs_op_pdafomi_pdaf) :: obs_op_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_pdafomi_pdaf) :: prodrinva_pdafomi
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdafomi_pdaf) :: obs_op_lin_pdafomi
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdafomi_pdaf) :: obs_op_adj_pdafomi
      ! Provide product R^-1 A
      procedure(c__prodrinva_l_pdafomi_pdaf) :: prodrinva_l_pdafomi
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdafomi_pdaf) :: init_dim_obs_l_pdafomi
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf,  &
         init_dim_obs_pdafomi, obs_op_pdafomi, prodrinva_pdafomi, cvt_ens_pdaf,  &
         cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi,  &
         obs_op_adj_pdafomi, prodrinva_l_pdafomi, init_n_domains_pdaf,  &
         init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf,  &
         l2g_state_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf_nondiagR

   SUBROUTINE c__PDAFomi_put_state_local(collect_state_pdaf,  &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf,  &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
      g2l_state_pdaf, l2g_state_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf_pdaf) :: obs_op_f_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_pdaf_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf_pdaf) :: init_dim_l_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf_pdaf) :: init_dim_obs_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf_pdaf) :: l2g_state_pdaf

      call PDAFomi_put_state_local(collect_state_pdaf, init_dim_obs_f_pdaf,  &
         obs_op_f_pdaf, prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,  &
         init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_local

   SUBROUTINE c__PDAFomi_put_state_global(collect_state_pdaf,  &
      init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_global(collect_state_pdaf, init_dim_obs_pdaf,  &
         obs_op_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_global

   SUBROUTINE c__PDAFomi_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdaf,  &
      obs_op_pdaf, prepoststep_pdaf, localize_covar_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf_pdaf) :: localize_covar_pdaf

      call PDAFomi_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdaf,  &
         obs_op_pdaf, prepoststep_pdaf, localize_covar_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_lenkf

   SUBROUTINE c__PDAFomi_put_state_ensrf(collect_state_pdaf, init_dim_obs_pdaf,  &
      obs_op_pdaf, localize_covar_serial_pdaf, prepoststep_pdaf, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf_pdaf) :: obs_op_pdaf
      ! Apply localization to HP and HXY
      procedure(c__localize_covar_serial_pdaf_pdaf) :: localize_covar_serial_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_ensrf(collect_state_pdaf, init_dim_obs_pdaf,  &
         obs_op_pdaf, localize_covar_serial_pdaf, prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAFomi_put_state_ensrf

   SUBROUTINE c__PDAFomi_put_state_generate_obs(collect_state_pdaf,  &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, get_obs_f_pdaf, prepoststep_pdaf,  &
      outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(out) :: outflag

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf_pdaf) :: collect_state_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf_pdaf) :: init_dim_obs_f_pdaf
      ! Observation operator
      procedure(c__obs_op_f_pdaf_pdaf) :: obs_op_f_pdaf
      ! Initialize observation vector
      procedure(c__get_obs_f_pdaf_pdaf) :: get_obs_f_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf_pdaf) :: prepoststep_pdaf

      call PDAFomi_put_state_generate_obs(collect_state_pdaf,  &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, get_obs_f_pdaf, prepoststep_pdaf,  &
         outflag)

   END SUBROUTINE c__PDAFomi_put_state_generate_obs
END MODULE pdafomi_c_put
