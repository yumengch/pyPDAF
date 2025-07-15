MODULE pdaf_c_internal
use iso_c_binding, only: c_int, c_double, c_bool, c_char
use pdaf_c_cb_interface
implicit none

contains
   SUBROUTINE c__PDAF_set_forget_local(domain, step, dim_obs_l, dim_ens, hx_l,  &
         hxbar_l, obs_l, u_init_obsvar_l, forget, aforget) bind(c)
         use PDAF_analysis_utils
         implicit none
         ! Current local analysis domain
         INTEGER(c_int), INTENT(in) :: domain
         ! Current time step
         INTEGER(c_int), INTENT(in) :: step
         ! Dimension of local observation vector
         INTEGER(c_int), INTENT(in) :: dim_obs_l
         ! Ensemble size
         INTEGER(c_int), INTENT(in) :: dim_ens
         ! Local observed ensemble
         REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
         ! Local observed state estimate
         REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
         ! Local observation vector
         REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
         ! Prescribed forgetting factor
         REAL(c_double), INTENT(in) :: forget
         ! Adaptive forgetting factor
         REAL(c_double), INTENT(out) :: aforget

         ! Initialize local mean obs. error variance
         procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l

         call PDAF_set_forget_local(domain, step, dim_obs_l, dim_ens, hx_l,  &
            hxbar_l, obs_l, u_init_obsvar_l, forget, aforget)

      END SUBROUTINE c__PDAF_set_forget_local

   SUBROUTINE c__PDAF_fcst_operations(step, u_collect_state,  &
      u_distribute_state, u_init_dim_obs, u_obs_op, u_init_obs, outflag) bind(c)
      use PDAF_forecast
      implicit none
      ! Time step in current forecast phase
      INTEGER(c_int), INTENT(in) :: step
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
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs

      call PDAF_fcst_operations(step, u_collect_state, u_distribute_state,  &
         u_init_dim_obs, u_obs_op, u_init_obs, outflag)

   END SUBROUTINE c__PDAF_fcst_operations

   SUBROUTINE c__PDAF_letkf_ana_T(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
      state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget,  &
      u_prodrinva_l, type_trans, screen, debug, flag) bind(c)
      use PDAF_letkf_analysis_T
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Local forecast state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! on exit: local weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hz_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: rndmat
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_letkf_ana_T(domain_p, step, dim_l, dim_obs_l, dim_ens, state_l,  &
         ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget, u_prodrinva_l,  &
         type_trans, screen, debug, flag)

   END SUBROUTINE c__PDAF_letkf_ana_T

   SUBROUTINE c__PDAFseik_update(step, dim_p, dim_obs_p, dim_ens, rank,  &
      state_p, uinv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_init_obsvar, u_prepoststep, screen, subtype, flag) bind(c)
      use PDAF_seik_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for SEIK analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFseik_update(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
         uinv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_init_obsvar, u_prepoststep, screen, subtype, flag)

   END SUBROUTINE c__PDAFseik_update

   SUBROUTINE c__PDAF3dvar_update(step, dim_p, dim_obs_p, dim_ens, dim_cvec,  &
      state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_prepoststep, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj, screen,  &
      subtype, flag) bind(c)
      use PDAF_3dvar_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Not used in 3D-Var
      REAL(c_double), DIMENSION(1, 1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for 3DVAR analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply control vector transform matrix
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF3dvar_update(step, dim_p, dim_obs_p, dim_ens, dim_cvec, state_p,  &
         ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_prepoststep, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj, screen,  &
         subtype, flag)

   END SUBROUTINE c__PDAF3dvar_update

   SUBROUTINE c__PDAFen3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, u_init_obsvar, screen, subtype, flag) bind(c)
      use PDAF_en3dvar_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Transform matrix
      REAL(c_double), DIMENSION(dim_ens-1, dim_ens-1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for 3DVAR analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
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

      call PDAFen3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens,  &
         dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens,  &
         u_obs_op_lin, u_obs_op_adj, u_init_obsvar, screen, subtype, flag)

   END SUBROUTINE c__PDAFen3dvar_update_estkf

   SUBROUTINE c__PDAFen3dvar_update_lestkf(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs,  &
      u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      screen, subtype, flag) bind(c)
      use PDAF_en3dvar_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Transform matrix
      REAL(c_double), DIMENSION(dim_ens-1, dim_ens-1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for 3DVAR analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
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

      call PDAFen3dvar_update_lestkf(step, dim_p, dim_obs_p, dim_ens,  &
         dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens,  &
         u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f,  &
         u_init_obs_f, u_init_obs_l, u_prodrinva_l, u_init_n_domains_p,  &
         u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
         u_init_obsvar, u_init_obsvar_l, screen, subtype, flag)

   END SUBROUTINE c__PDAFen3dvar_update_lestkf

   SUBROUTINE c__PDAFetkf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_init_obsvar, u_prepoststep, screen, subtype, dim_lag, sens_p,  &
      cnt_maxlag, flag) bind(c)
      use PDAF_etkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for ETKF analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFetkf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
         ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_init_obsvar, u_prepoststep, screen, subtype, dim_lag, sens_p,  &
         cnt_maxlag, flag)

   END SUBROUTINE c__PDAFetkf_update

   SUBROUTINE c__PDAF_netf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
      rndmat, t, type_forget, forget, type_winf, limit_winf, type_noise,  &
      noise_amp, hz_p, obs_p, u_likelihood, screen, debug, flag) bind(c)
      use PDAF_netf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Orthogonal random matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: rndmat
      ! Ensemble transform matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: t
      ! Type of forgetting factor
      INTEGER(c_int), INTENT(in) :: type_forget
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Type of weights inflation
      INTEGER(c_int), INTENT(in) :: type_winf
      ! Limit for weights inflation
      REAL(c_double), INTENT(in) :: limit_winf
      ! Type of pertubing noise
      INTEGER(c_int), INTENT(in) :: type_noise
      ! Amplitude of noise
      REAL(c_double), INTENT(in) :: noise_amp
      ! Temporary matrices for analysis
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: hz_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood

      call PDAF_netf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         rndmat, t, type_forget, forget, type_winf, limit_winf, type_noise,  &
         noise_amp, hz_p, obs_p, u_likelihood, screen, debug, flag)

   END SUBROUTINE c__PDAF_netf_ana

   SUBROUTINE c__PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, ens_p,  &
      rndmat, t, u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood, screen,  &
      flag) bind(c)
      use PDAF_netf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Orthogonal random matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: rndmat
      ! Ensemble transform matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: t
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood

      call PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, ens_p, rndmat,  &
         t, u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood, screen, flag)

   END SUBROUTINE c__PDAF_netf_smootherT

   SUBROUTINE c__PDAF_smoother_netf(dim_p, dim_ens, dim_lag, ainv, sens_p,  &
      cnt_maxlag, screen) bind(c)
      use PDAF_netf_analysis
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: ainv
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_smoother_netf(dim_p, dim_ens, dim_lag, ainv, sens_p,  &
         cnt_maxlag, screen)

   END SUBROUTINE c__PDAF_smoother_netf

   SUBROUTINE c__PDAF_lnetf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
      ens_l, hx_l, obs_l, rndmat, u_likelihood_l, type_forget, forget,  &
      type_winf, limit_winf, cnt_small_svals, eff_dimens, t, screen, debug,  &
      flag) bind(c)
      use PDAF_lnetf_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: rndmat
      ! Typ eof forgetting factor
      INTEGER(c_int), INTENT(in) :: type_forget
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Type of weights inflation
      INTEGER(c_int), INTENT(in) :: type_winf
      ! Limit for weights inflation
      REAL(c_double), INTENT(in) :: limit_winf
      ! Number of small eigen values
      INTEGER(c_int), INTENT(inout) :: cnt_small_svals
      ! Effective ensemble size
      REAL(c_double), DIMENSION(1), INTENT(inout) :: eff_dimens
      ! local ensemble transformation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: t
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l

      call PDAF_lnetf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, ens_l,  &
         hx_l, obs_l, rndmat, u_likelihood_l, type_forget, forget, type_winf,  &
         limit_winf, cnt_small_svals, eff_dimens, t, screen, debug, flag)

   END SUBROUTINE c__PDAF_lnetf_ana

   SUBROUTINE c__PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l,  &
      dim_ens, hx_f, rndmat, u_g2l_obs, u_init_obs_l, u_likelihood_l, screen,  &
      t, flag) bind(c)
      use PDAF_lnetf_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of full observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local full observed state ens.
      REAL(c_double), DIMENSION(dim_obs_f, dim_ens), INTENT(in) :: hx_f
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: rndmat
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! local ensemble transformation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: t
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l

      call PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l, dim_ens,  &
         hx_f, rndmat, u_g2l_obs, u_init_obs_l, u_likelihood_l, screen, t, flag)

   END SUBROUTINE c__PDAF_lnetf_smootherT

   SUBROUTINE c__PDAF_smoother_lnetf(domain_p, step, dim_p, dim_l, dim_ens,  &
      dim_lag, ainv, ens_l, sens_p, cnt_maxlag, u_g2l_state, u_l2g_state,  &
      screen) bind(c)
      use PDAF_lnetf_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: ainv
      ! local past ensemble (temporary)
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Get state on local ana. domain from global state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state

      call PDAF_smoother_lnetf(domain_p, step, dim_p, dim_l, dim_ens, dim_lag,  &
         ainv, ens_l, sens_p, cnt_maxlag, u_g2l_state, u_l2g_state, screen)

   END SUBROUTINE c__PDAF_smoother_lnetf

   SUBROUTINE c__PDAF_memcount_ini(ncounters) bind(c)
      use PDAF_memcounting
      implicit none
      ! Number of memory counters
      INTEGER(c_int), INTENT(in) :: ncounters


      call PDAF_memcount_ini(ncounters)

   END SUBROUTINE c__PDAF_memcount_ini

   SUBROUTINE c__PDAF_memcount_define(stortype, wordlength) bind(c)
      use PDAF_memcounting
      implicit none
      ! Type of variable
      CHARACTER(kind=c_char), INTENT(in) :: stortype
      ! Word length for chosen type
      INTEGER(c_int), INTENT(in) :: wordlength


      call PDAF_memcount_define(stortype, wordlength)

   END SUBROUTINE c__PDAF_memcount_define

   SUBROUTINE c__PDAF_memcount(id, stortype, dim) bind(c)
      use PDAF_memcounting
      implicit none
      ! Id of the counter
      INTEGER(c_int), INTENT(in) :: id
      ! Type of variable
      CHARACTER(kind=c_char), INTENT(in) :: stortype
      ! Dimension of allocated variable
      INTEGER(c_int), INTENT(in) :: dim


      call PDAF_memcount(id, stortype, dim)

   END SUBROUTINE c__PDAF_memcount

   SUBROUTINE c__PDAF_init_filters(type_filter, subtype, param_int, dim_pint,  &
      param_real, dim_preal, filterstr, ensemblefilter, fixedbasis, screen,  &
      flag) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Type of filter
      INTEGER(c_int), INTENT(in) :: type_filter
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Name of filter algorithm
      CHARACTER(kind=c_char), INTENT(out) :: filterstr
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_init_filters(type_filter, subtype, param_int, dim_pint,  &
         param_real, dim_preal, filterstr, ensemblefilter_out, fixedbasis_out, screen,  &
         flag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_init_filters

   SUBROUTINE c__PDAF_alloc_filters(filterstr, subtype, flag) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Name of filter algorithm
      CHARACTER(kind=c_char), INTENT(in) :: filterstr
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_alloc_filters(filterstr, subtype, flag)

   END SUBROUTINE c__PDAF_alloc_filters

   SUBROUTINE c__PDAF_configinfo_filters(subtype, verbose) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_configinfo_filters(subtype, verbose)

   END SUBROUTINE c__PDAF_configinfo_filters

   SUBROUTINE c__PDAF_options_filters(type_filter) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Type of filter
      INTEGER(c_int), INTENT(in) :: type_filter


      call PDAF_options_filters(type_filter)

   END SUBROUTINE c__PDAF_options_filters

   SUBROUTINE c__PDAF_print_info_filters(printtype) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_print_info_filters(printtype)

   END SUBROUTINE c__PDAF_print_info_filters

   SUBROUTINE c__PDAF_allreduce(val_p, val_g, mpitype, mpiop, status) bind(c)
      use PDAF_comm_obs
      implicit none
      ! PE-local value
      INTEGER(c_int), INTENT(in) :: val_p
      ! reduced global value
      INTEGER(c_int), INTENT(out) :: val_g
      ! MPI data type
      INTEGER(c_int), INTENT(in) :: mpitype
      ! MPI operator
      INTEGER(c_int), INTENT(in) :: mpiop
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_allreduce(val_p, val_g, mpitype, mpiop, status)

   END SUBROUTINE c__PDAF_allreduce

   SUBROUTINE c__PDAFlseik_update(step, dim_p, dim_obs_f, dim_ens, rank,  &
      state_p, uinv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_prepoststep, screen, subtype, flag) bind(c)
      use PDAF_lseik_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute product of R^(-1) with HV
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from global state
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

      call PDAFlseik_update(step, dim_p, dim_obs_f, dim_ens, rank, state_p,  &
         uinv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_prepoststep, screen, subtype, flag)

   END SUBROUTINE c__PDAFlseik_update

   SUBROUTINE c__PDAF_ensrf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_ensrf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_ensrf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_ensrf_init

   SUBROUTINE c__PDAF_ensrf_alloc(outflag) bind(c)
      use PDAF_ensrf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_ensrf_alloc(outflag)

   END SUBROUTINE c__PDAF_ensrf_alloc

   SUBROUTINE c__PDAF_ensrf_config(subtype, verbose) bind(c)
      use PDAF_ensrf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_ensrf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_ensrf_config

   SUBROUTINE c__PDAF_ensrf_set_iparam(id, value, flag) bind(c)
      use PDAF_ensrf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_ensrf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_ensrf_set_iparam

   SUBROUTINE c__PDAF_ensrf_set_rparam(id, value, flag) bind(c)
      use PDAF_ensrf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_ensrf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_ensrf_set_rparam

   SUBROUTINE c__PDAF_ensrf_options() bind(c)
      use PDAF_ensrf
      implicit none
      call PDAF_ensrf_options()

   END SUBROUTINE c__PDAF_ensrf_options

   SUBROUTINE c__PDAF_ensrf_memtime(printtype) bind(c)
      use PDAF_ensrf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_ensrf_memtime(printtype)

   END SUBROUTINE c__PDAF_ensrf_memtime

   SUBROUTINE c__PDAF_estkf_ana_fixed(step, dim_p, dim_obs_p, dim_ens, rank,  &
      state_p, ainv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
      type_sqrt, debug, flag) bind(c)
      use PDAF_estkf_analysis_fixed
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! on exit: PE-local forecast mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix A - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: ainv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hl_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 with some matrix
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_estkf_ana_fixed(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
         ainv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
         type_sqrt, debug, flag)

   END SUBROUTINE c__PDAF_estkf_ana_fixed

   SUBROUTINE c__PDAF_etkf_ana_fixed(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ainv, ens_p, hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen, debug,  &
      flag) bind(c)
      use PDAF_etkf_analysis_fixed
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! on exit: weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hz_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_etkf_ana_fixed(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
         ens_p, hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen, debug, flag)

   END SUBROUTINE c__PDAF_etkf_ana_fixed

   SUBROUTINE c__PDAFestkf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
      u_init_obsvar, u_prepoststep, screen, subtype, envar_mode, dim_lag,  &
      sens_p, cnt_maxlag, flag) bind(c)
      use PDAF_estkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of transform matrix A
      REAL(c_double), DIMENSION(dim_ens-1, dim_ens-1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Flag whether routine is called from 3DVar for special functionality
      INTEGER(c_int), INTENT(in) :: envar_mode
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for ESTKF analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFestkf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
         ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_prodrinva,  &
         u_init_obsvar, u_prepoststep, screen, subtype, envar_mode, dim_lag,  &
         sens_p, cnt_maxlag, flag)

   END SUBROUTINE c__PDAFestkf_update

   SUBROUTINE c__PDAFlknetf_update_step(step, dim_p, dim_obs_f, dim_ens,  &
      state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_likelihood_l, u_likelihood_hyb_l, u_prepoststep, screen, subtype,  &
      flag) bind(c)
      use PDAF_lknetf_update_step

      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute product of R^(-1) with HV with hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: u_prodrinva_hyb_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from global state
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFlknetf_update_step(step, dim_p, dim_obs_f, dim_ens, state_p,  &
         ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
         u_prodrinva_hyb_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_likelihood_l, u_likelihood_hyb_l, u_prepoststep, screen, subtype, flag)

   END SUBROUTINE c__PDAFlknetf_update_step

   SUBROUTINE c__PDAFletkf_update(step, dim_p, dim_obs_f, dim_ens, state_p,  &
      ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_prepoststep, screen, subtype, dim_lag, sens_p, cnt_maxlag, flag) bind(c)
      use PDAF_letkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
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

      call PDAFletkf_update(step, dim_p, dim_obs_f, dim_ens, state_p, ainv,  &
         ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_prepoststep, screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

   END SUBROUTINE c__PDAFletkf_update

   SUBROUTINE c__PDAF_lseik_ana_trans(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, rank, state_l, uinv_l, ens_l, hl_l, hxbar_l, obs_l, omegat_in,  &
      forget, u_prodrinva_l, nm1vsn, type_sqrt, screen, debug, flag) bind(c)
      use PDAF_lseik_analysis_trans
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! on exit: state on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! Inverse of matrix U - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hl_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Matrix Omega
      REAL(c_double), DIMENSION(rank, dim_ens), INTENT(inout) :: omegat_in
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Whether covariance is normalized with 1/N or 1/(N-1)
      INTEGER(c_int), INTENT(in) :: nm1vsn
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_lseik_ana_trans(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         rank, state_l, uinv_l, ens_l, hl_l, hxbar_l, obs_l, omegat_in, forget,  &
         u_prodrinva_l, nm1vsn, type_sqrt, screen, debug, flag)

   END SUBROUTINE c__PDAF_lseik_ana_trans

   SUBROUTINE c__PDAF_en3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cvec_p,  &
      dim_obs_p, ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, screen) bind(c)
      use PDAF_en3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_en3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cvec_p,  &
         dim_obs_p, ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens,  &
         u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, screen)

   END SUBROUTINE c__PDAF_en3dvar_optim_lbfgs

   SUBROUTINE c__PDAF_en3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cvec_p,  &
      dim_obs_p, ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, screen) bind(c)
      use PDAF_en3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_en3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cvec_p,  &
         dim_obs_p, ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens,  &
         u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, screen)

   END SUBROUTINE c__PDAF_en3dvar_optim_cgplus

   SUBROUTINE c__PDAF_en3dvar_optim_cg(step, dim_p, dim_ens, dim_cvec_p,  &
      dim_obs_p, ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, screen) bind(c)
      use PDAF_en3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_en3dvar_optim_cg(step, dim_p, dim_ens, dim_cvec_p, dim_obs_p,  &
         ens_p, obs_p, dy_p, v_p, u_prodrinva, u_cvt_ens, u_cvt_adj_ens,  &
         u_obs_op_lin, u_obs_op_adj, opt_parallel, screen)

   END SUBROUTINE c__PDAF_en3dvar_optim_cg

   SUBROUTINE c__PDAF_en3dvar_costf_cvt(step, iter, dim_p, dim_ens, dim_cvec_p,  &
      dim_obs_p, ens_p, obs_p, dy_p, v_p, j_tot, gradj, u_prodrinva, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel) bind(c)
      use PDAF_en3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Optimization iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(in) :: v_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: PE-local gradient of J
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: gradj
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_en3dvar_costf_cvt(step, iter, dim_p, dim_ens, dim_cvec_p,  &
         dim_obs_p, ens_p, obs_p, dy_p, v_p, j_tot, gradj, u_prodrinva,  &
         u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel)

   END SUBROUTINE c__PDAF_en3dvar_costf_cvt

   SUBROUTINE c__PDAF_en3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens,  &
      dim_cvec_p, dim_obs_p, ens_p, obs_p, dy_p, v_p, d_p, j_tot, gradj,  &
      hessjd, u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel) bind(c)
      use PDAF_en3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Optimization iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(in) :: v_p
      ! CG descent direction
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: d_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: gradient of J
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: gradj
      ! on exit: Hessian of J times d_p
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: hessjd
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_en3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens, dim_cvec_p,  &
         dim_obs_p, ens_p, obs_p, dy_p, v_p, d_p, j_tot, gradj, hessjd,  &
         u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj,  &
         opt_parallel)

   END SUBROUTINE c__PDAF_en3dvar_costf_cg_cvt

   SUBROUTINE c__PDAF_gather_ens(dim_p, dim_ens_p, ens, screen) bind(c)
      use PDAF_communicate_ens
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: ens
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_gather_ens(dim_p, dim_ens_p, ens, screen)

   END SUBROUTINE c__PDAF_gather_ens

   SUBROUTINE c__PDAF_scatter_ens(dim_p, dim_ens_p, ens, state, screen) bind(c)
      use PDAF_communicate_ens
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: ens
      ! PE-local state vector (for SEEK)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_scatter_ens(dim_p, dim_ens_p, ens, state, screen)

   END SUBROUTINE c__PDAF_scatter_ens

   SUBROUTINE c__PDAF_hyb3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cv_par_p,  &
      dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
      u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel, beta_3dvar, screen) bind(c)
      use PDAF_hyb3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized)
      INTEGER(c_int), INTENT(in) :: dim_cv_par_p
      ! Size of control vector (ensemble)
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(inout) :: v_par_p
      ! Control vector (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: v_ens_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: beta_3dvar
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_hyb3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cv_par_p,  &
         dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
         u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
         u_obs_op_adj, opt_parallel, beta_3dvar, screen)

   END SUBROUTINE c__PDAF_hyb3dvar_optim_lbfgs

   SUBROUTINE c__PDAF_hyb3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cv_par_p,  &
      dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
      u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel, beta_3dvar, screen) bind(c)
      use PDAF_hyb3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized)
      INTEGER(c_int), INTENT(in) :: dim_cv_par_p
      ! Size of control vector (ensemble)
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(inout) :: v_par_p
      ! Control vector (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: v_ens_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: beta_3dvar
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_hyb3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cv_par_p,  &
         dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
         u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
         u_obs_op_adj, opt_parallel, beta_3dvar, screen)

   END SUBROUTINE c__PDAF_hyb3dvar_optim_cgplus

   SUBROUTINE c__PDAF_hyb3dvar_optim_cg(step, dim_p, dim_ens, dim_cv_par_p,  &
      dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
      u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel, beta_3dvar, screen) bind(c)
      use PDAF_hyb3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized)
      INTEGER(c_int), INTENT(in) :: dim_cv_par_p
      ! Size of control vector (ensemble)
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(inout) :: v_par_p
      ! Control vector (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: v_ens_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: beta_3dvar
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_hyb3dvar_optim_cg(step, dim_p, dim_ens, dim_cv_par_p,  &
         dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
         u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
         u_obs_op_adj, opt_parallel, beta_3dvar, screen)

   END SUBROUTINE c__PDAF_hyb3dvar_optim_cg

   SUBROUTINE c__PDAF_hyb3dvar_costf_cvt(step, iter, dim_p, dim_ens, dim_cv_p,  &
      dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p,  &
      v_ens_p, v_p, j_tot, gradj, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens,  &
      u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, beta) bind(c)
      use PDAF_hyb3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Optimization iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (full)
      INTEGER(c_int), INTENT(in) :: dim_cv_p
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cv_par_p
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(inout) :: v_par_p
      ! Control vector (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: v_ens_p
      ! Control vector (full)
      REAL(c_double), DIMENSION(dim_cv_p), INTENT(in) :: v_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: PE-local gradient of J (full)
      REAL(c_double), DIMENSION(dim_cv_p), INTENT(out) :: gradj
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: beta

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_hyb3dvar_costf_cvt(step, iter, dim_p, dim_ens, dim_cv_p,  &
         dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p,  &
         v_ens_p, v_p, j_tot, gradj, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens,  &
         u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, opt_parallel, beta)

   END SUBROUTINE c__PDAF_hyb3dvar_costf_cvt

   SUBROUTINE c__PDAF_hyb3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens,  &
      dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p,  &
      v_ens_p, d_par_p, d_ens_p, j_tot, gradj_par, gradj_ens, hessjd_par,  &
      hessjd_ens, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens,  &
      u_obs_op_lin, u_obs_op_adj, opt_parallel, beta) bind(c)
      use PDAF_hyb3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Optimization iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cv_par_p
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(in) :: v_par_p
      ! Control vector (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(in) :: v_ens_p
      ! CG descent direction (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(inout) :: d_par_p
      ! CG descent direction (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: d_ens_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: gradient of J (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(out) :: gradj_par
      ! on exit: gradient of J (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(out) :: gradj_ens
      ! on exit: Hessian of J times d_p (parameterized part)
      REAL(c_double), DIMENSION(dim_cv_par_p), INTENT(out) :: hessjd_par
      ! on exit: Hessian of J times d_p (ensemble part)
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(out) :: hessjd_ens
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: beta

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAF_hyb3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens, dim_cv_par_p,  &
         dim_cv_ens_p, dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p,  &
         d_par_p, d_ens_p, j_tot, gradj_par, gradj_ens, hessjd_par, hessjd_ens,  &
         u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
         u_obs_op_adj, opt_parallel, beta)

   END SUBROUTINE c__PDAF_hyb3dvar_costf_cg_cvt

   SUBROUTINE c__PDAF_print_version() bind(c)
      use PDAF_info
      implicit none
      call PDAF_print_version()

   END SUBROUTINE c__PDAF_print_version

   SUBROUTINE c__PDAFen3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec_ens, state_p, ens_p, state_inc_p, hxbar_p, obs_p, u_prodrinva,  &
      u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj, screen, type_opt,  &
      debug, flag) bind(c)
      use PDAF_en3dvar_analysis_cvt
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local state analysis increment
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_inc_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of minimizer for 3DVar
      INTEGER(c_int), INTENT(in) :: type_opt
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAFen3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens,  &
         dim_cvec_ens, state_p, ens_p, state_inc_p, hxbar_p, obs_p,  &
         u_prodrinva, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin, u_obs_op_adj,  &
         screen, type_opt, debug, flag)

   END SUBROUTINE c__PDAFen3dvar_analysis_cvt

   SUBROUTINE c__PDAF_sisort(n, veca) bind(c)
      use PDAF_diag
      implicit none
      !
      INTEGER(c_int), INTENT(in) :: n
      !
      REAL(c_double), DIMENSION(n), INTENT(inout) :: veca


      call PDAF_sisort(n, veca)

   END SUBROUTINE c__PDAF_sisort

   ! SUBROUTINE c__PDAF_unbiased_moments_from_summed_residuals(dim_ens, dim_p,  &
   !    kmax, sum_expo_resid, moments) bind(c)
   !    ! number of ensemble members/samples
   !    INTEGER(c_int), INTENT(in) :: dim_ens
   !    ! local size of the state
   !    INTEGER(c_int), INTENT(in) :: dim_p
   !    ! maximum order of central moment that is computed, maximum is 4
   !    INTEGER(c_int), INTENT(in) :: kmax
   !    ! sum of exponentiated residulals
   !    REAL(c_double), DIMENSION(dim_p, kmax), INTENT(in) :: sum_expo_resid
   !    ! The columns contain the moments of the ensemble
   !    REAL(c_double), DIMENSION(dim_p, kmax), INTENT(inout) :: moments


   !    call PDAF_unbiased_moments_from_summed_residuals(dim_ens, dim_p, kmax,  &
   !       sum_expo_resid, moments)

   ! END SUBROUTINE c__PDAF_unbiased_moments_from_summed_residuals

   ! SUBROUTINE c__PDAF_biased_moments_from_summed_residuals(dim_ens, dim_p,  &
   !    kmax, sum_expo_resid, moments) bind(c)
   !    ! number of ensemble members/samples
   !    INTEGER(c_int), INTENT(in) :: dim_ens
   !    ! local size of the state
   !    INTEGER(c_int), INTENT(in) :: dim_p
   !    ! maximum order of central moment that is computed, maximum is 4
   !    INTEGER(c_int), INTENT(in) :: kmax
   !    ! sum of exponentiated residulals
   !    REAL(c_double), DIMENSION(dim_p, kmax), INTENT(in) :: sum_expo_resid
   !    ! The columns contain the moments of the ensemble
   !    REAL(c_double), DIMENSION(dim_p, kmax), INTENT(inout) :: moments


   !    call PDAF_biased_moments_from_summed_residuals(dim_ens, dim_p, kmax,  &
   !       sum_expo_resid, moments)

   ! END SUBROUTINE c__PDAF_biased_moments_from_summed_residuals

   SUBROUTINE c__PDAF_enkf_ana_rlm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
      state_p, ens_p, hzb, hx_p, hxbar_p, obs_p, u_add_obs_err,  &
      u_init_obs_covar, screen, debug, flag) bind(c)
      use PDAF_enkf_analysis_rlm
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank to be considered for inversion of HPH
      INTEGER(c_int), INTENT(in) :: rank_ana
      ! PE-local ensemble mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Ensemble tranformation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: hzb
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: hx_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar

      call PDAF_enkf_ana_rlm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
         state_p, ens_p, hzb, hx_p, hxbar_p, obs_p, u_add_obs_err,  &
         u_init_obs_covar, screen, debug, flag)

   END SUBROUTINE c__PDAF_enkf_ana_rlm

   SUBROUTINE c__PDAF_smoother_enkf(dim_p, dim_ens, dim_lag, ainv, sens_p,  &
      cnt_maxlag, forget, screen) bind(c)
      use PDAF_enkf_analysis_rlm
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: ainv
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_smoother_enkf(dim_p, dim_ens, dim_lag, ainv, sens_p,  &
         cnt_maxlag, forget, screen)

   END SUBROUTINE c__PDAF_smoother_enkf

   SUBROUTINE c__PDAFensrf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obsvars,  &
      u_localize_covar_serial, u_prepoststep, screen, subtype, flag) bind(c)
      use PDAF_ensrf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Specification of filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Initialize vector of observation error variances
      procedure(c__init_obsvar_pdaf) :: u_init_obsvars
      ! Apply localization for single-observation vectors
      procedure(c__localize_covar_serial_pdaf) :: u_localize_covar_serial
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFensrf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_init_obsvars,  &
         u_localize_covar_serial, u_prepoststep, screen, subtype, flag)

   END SUBROUTINE c__PDAFensrf_update

   SUBROUTINE c__PDAF_pf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
      type_resample, type_winf, limit_winf, type_noise, noise_amp, hz_p, obs_p,  &
      u_likelihood, screen, debug, flag) bind(c)
      use PDAF_pf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local forecast mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Type of resampling scheme
      INTEGER(c_int), INTENT(in) :: type_resample
      ! Type of weights inflation
      INTEGER(c_int), INTENT(in) :: type_winf
      ! Limit for weights inflation
      REAL(c_double), INTENT(in) :: limit_winf
      ! Type of pertubing noise
      INTEGER(c_int), INTENT(in) :: type_noise
      ! Amplitude of noise
      REAL(c_double), INTENT(in) :: noise_amp
      ! Temporary matrices for analysis
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: hz_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood

      call PDAF_pf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         type_resample, type_winf, limit_winf, type_noise, noise_amp, hz_p,  &
         obs_p, u_likelihood, screen, debug, flag)

   END SUBROUTINE c__PDAF_pf_ana

   SUBROUTINE c__PDAF_pf_resampling(method, nin, nout, weights, ids,  &
      screen) bind(c)
      use PDAF_pf_analysis
      implicit none
      ! Choose resampling method
      INTEGER(c_int), INTENT(in) :: method
      ! number of particles
      INTEGER(c_int), INTENT(in) :: nin
      ! number of particles to be resampled
      INTEGER(c_int), INTENT(in) :: nout
      ! Weights
      REAL(c_double), DIMENSION(Nin), INTENT(in) :: weights
      ! Indices of resampled ensmeble states
      INTEGER(c_int), DIMENSION(Nout), INTENT(out) :: ids
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_pf_resampling(method, nin, nout, weights, ids, screen)

   END SUBROUTINE c__PDAF_pf_resampling

   SUBROUTINE c__PDAF_mvnormalize(mode, dim_state, dim_field, offset, ncol,  &
      states, stddev, status) bind(c)
      use PDAF_sample
      implicit none
      ! Mode: (1) normalize, (2) re-scale
      INTEGER(c_int), INTENT(in) :: mode
      ! Dimension of state vector
      INTEGER(c_int), INTENT(in) :: dim_state
      ! Dimension of a field in state vector
      INTEGER(c_int), INTENT(in) :: dim_field
      ! Offset of field in state vector
      INTEGER(c_int), INTENT(in) :: offset
      ! Number of columns in array states
      INTEGER(c_int), INTENT(in) :: ncol
      ! State vector array
      REAL(c_double), DIMENSION(dim_state, ncol), INTENT(inout) :: states
      ! Standard deviation of field
      REAL(c_double), INTENT(inout) :: stddev
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_mvnormalize(mode, dim_state, dim_field, offset, ncol, states,  &
         stddev, status)

   END SUBROUTINE c__PDAF_mvnormalize

   SUBROUTINE c__PDAF_3dvar_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_3dvar
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_3dvar_init

   SUBROUTINE c__PDAF_3dvar_alloc(subtype, outflag) bind(c)
      use PDAF_3dvar
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_3dvar_alloc(subtype, outflag)

   END SUBROUTINE c__PDAF_3dvar_alloc

   SUBROUTINE c__PDAF_3dvar_config(subtype, verbose) bind(c)
      use PDAF_3dvar
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_3dvar_config(subtype, verbose)

   END SUBROUTINE c__PDAF_3dvar_config

   SUBROUTINE c__PDAF_3dvar_set_iparam(id, value, flag) bind(c)
      use PDAF_3dvar
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_3dvar_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_3dvar_set_iparam

   SUBROUTINE c__PDAF_3dvar_set_rparam(id, value, flag) bind(c)
      use PDAF_3dvar
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_3dvar_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_3dvar_set_rparam

   SUBROUTINE c__PDAF_3dvar_options() bind(c)
      use PDAF_3dvar
      implicit none
      call PDAF_3dvar_options()

   END SUBROUTINE c__PDAF_3dvar_options

   SUBROUTINE c__PDAF_3dvar_memtime(printtype) bind(c)
      use PDAF_3dvar
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_3dvar_memtime(printtype)

   END SUBROUTINE c__PDAF_3dvar_memtime


   SUBROUTINE c__PDAF_reset_dim_ens(dim_ens_in, outflag) bind(c)
      use PDAF_set
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: dim_ens_in
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_reset_dim_ens(dim_ens_in, outflag)

   END SUBROUTINE c__PDAF_reset_dim_ens

   SUBROUTINE c__PDAF_reset_dim_p(dim_p_in, outflag) bind(c)
      use PDAF_set
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: dim_p_in
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_reset_dim_p(dim_p_in, outflag)

   END SUBROUTINE c__PDAF_reset_dim_p

   SUBROUTINE c__PDAF_3dvar_optim_lbfgs(step, dim_p, dim_cvec_p, dim_obs_p,  &
      obs_p, dy_p, v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel, screen) bind(c)
      use PDAF_3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

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

      call PDAF_3dvar_optim_lbfgs(step, dim_p, dim_cvec_p, dim_obs_p, obs_p,  &
         dy_p, v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
         opt_parallel, screen)

   END SUBROUTINE c__PDAF_3dvar_optim_lbfgs

   SUBROUTINE c__PDAF_3dvar_optim_cgplus(step, dim_p, dim_cvec_p, dim_obs_p,  &
      obs_p, dy_p, v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin,  &
      u_obs_op_adj, opt_parallel, screen) bind(c)
      use PDAF_3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

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

      call PDAF_3dvar_optim_cgplus(step, dim_p, dim_cvec_p, dim_obs_p, obs_p,  &
         dy_p, v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
         opt_parallel, screen)

   END SUBROUTINE c__PDAF_3dvar_optim_cgplus

   SUBROUTINE c__PDAF_3dvar_optim_cg(step, dim_p, dim_cvec_p, dim_obs_p, obs_p,  &
      dy_p, v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
      opt_parallel, screen) bind(c)
      use PDAF_3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: v_p
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

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

      call PDAF_3dvar_optim_cg(step, dim_p, dim_cvec_p, dim_obs_p, obs_p, dy_p,  &
         v_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
         opt_parallel, screen)

   END SUBROUTINE c__PDAF_3dvar_optim_cg

   SUBROUTINE c__PDAF_3dvar_costf_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p,  &
      obs_p, dy_p, v_p, j_tot, gradj, u_prodrinva, u_cvt, u_cvt_adj,  &
      u_obs_op_lin, u_obs_op_adj, opt_parallel) bind(c)
      use PDAF_3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Optimization iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(in) :: v_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: PE-local gradient of J
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: gradj
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel

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

      call PDAF_3dvar_costf_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p,  &
         obs_p, dy_p, v_p, j_tot, gradj, u_prodrinva, u_cvt, u_cvt_adj,  &
         u_obs_op_lin, u_obs_op_adj, opt_parallel)

   END SUBROUTINE c__PDAF_3dvar_costf_cvt

   SUBROUTINE c__PDAF_3dvar_costf_cg_cvt(step, iter, dim_p, dim_cvec_p,  &
      dim_obs_p, obs_p, dy_p, v_p, d_p, j_tot, gradj, hessjd, u_prodrinva,  &
      u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj, opt_parallel) bind(c)
      use PDAF_3dvar_optim
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! CG iteration
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Background innovation
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: dy_p
      ! Control vector
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(in) :: v_p
      ! CG descent direction
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(inout) :: d_p
      ! on exit: Value of cost function
      REAL(c_double), INTENT(out) :: j_tot
      ! on exit: gradient of J
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: gradj
      ! on exit: Hessian of J times d_p
      REAL(c_double), DIMENSION(dim_cvec_p), INTENT(out) :: hessjd
      ! Whether to use a decomposed control vector
      INTEGER(c_int), INTENT(in) :: opt_parallel

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

      call PDAF_3dvar_costf_cg_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p,  &
         obs_p, dy_p, v_p, d_p, j_tot, gradj, hessjd, u_prodrinva, u_cvt,  &
         u_cvt_adj, u_obs_op_lin, u_obs_op_adj, opt_parallel)

   END SUBROUTINE c__PDAF_3dvar_costf_cg_cvt

   SUBROUTINE c__PDAF_lknetf_analysis_T(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, state_l, ainv_l, ens_l, hx_l, hxbar_l, obs_l, rndmat, forget,  &
      u_prodrinva_l, u_init_obsvar_l, u_likelihood_l, screen, type_forget,  &
      eff_dimens, type_hyb, hyb_g, hyb_k, gamma, skew_mabs, kurt_mabs,  &
      flag) bind(c)
      use PDAF_lknetf_analysis_sync
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! local forecast state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! on exit: local weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! local observed state ens.
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
      ! local observed ens. mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: rndmat
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of forgetting factor
      INTEGER(c_int), INTENT(in) :: type_forget
      ! Effective ensemble size
      REAL(c_double), DIMENSION(1), INTENT(inout) :: eff_dimens
      ! Type of hybrid weight
      INTEGER(c_int), INTENT(in) :: type_hyb
      ! Prescribed hybrid weight for state transformation
      REAL(c_double), INTENT(in) :: hyb_g
      ! Scale factor kappa (for type_hyb 3 and 4)
      REAL(c_double), INTENT(in) :: hyb_k
      ! Hybrid weight for state transformation
      REAL(c_double), DIMENSION(1), INTENT(inout) :: gamma
      ! Mean absolute skewness
      REAL(c_double), DIMENSION(1), INTENT(inout) :: skew_mabs
      ! Mean absolute kurtosis
      REAL(c_double), DIMENSION(1), INTENT(inout) :: kurt_mabs
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l
      ! Provide likelihood of an ensemble state
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l

      call PDAF_lknetf_analysis_T(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         state_l, ainv_l, ens_l, hx_l, hxbar_l, obs_l, rndmat, forget,  &
         u_prodrinva_l, u_init_obsvar_l, u_likelihood_l, screen, type_forget,  &
         eff_dimens, type_hyb, hyb_g, hyb_k, gamma, skew_mabs, kurt_mabs, flag)

   END SUBROUTINE c__PDAF_lknetf_analysis_T

   SUBROUTINE c__PDAF_get_ensstats(skew_ptr, kurt_ptr, status) bind(c)
      use PDAF_get
      implicit none
      ! Pointer to skewness array
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: skew_ptr
      ! Pointer to kurtosis array
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: kurt_ptr
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_get_ensstats(skew_ptr, kurt_ptr, status)
   END SUBROUTINE c__PDAF_get_ensstats

   SUBROUTINE c__PDAF_estkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_estkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_estkf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_estkf_init

   SUBROUTINE c__PDAF_estkf_alloc(outflag) bind(c)
      use PDAF_estkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_estkf_alloc(outflag)

   END SUBROUTINE c__PDAF_estkf_alloc

   SUBROUTINE c__PDAF_estkf_config(subtype, verbose) bind(c)
      use PDAF_estkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_estkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_estkf_config

   SUBROUTINE c__PDAF_estkf_set_iparam(id, value, flag) bind(c)
      use PDAF_estkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_estkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_estkf_set_iparam

   SUBROUTINE c__PDAF_estkf_set_rparam(id, value, flag) bind(c)
      use PDAF_estkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_estkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_estkf_set_rparam

   SUBROUTINE c__PDAF_estkf_options() bind(c)
      use PDAF_estkf
      implicit none
      call PDAF_estkf_options()

   END SUBROUTINE c__PDAF_estkf_options

   SUBROUTINE c__PDAF_estkf_memtime(printtype) bind(c)
      use PDAF_estkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_estkf_memtime(printtype)

   END SUBROUTINE c__PDAF_estkf_memtime

   SUBROUTINE c__PDAF_gen_obs(step, dim_p, dim_obs_f, dim_ens, state_p, ainv,  &
      ens_p, u_init_dim_obs_f, u_obs_op_f, u_get_obs_f, u_init_obserr_f,  &
      u_prepoststep, screen, flag) bind(c)
      use PDAF_generate_obs_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: u_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: u_obs_op_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: u_get_obs_f
      ! Initialize vector of observation error standard deviations
      procedure(c__init_obserr_f_pdaf) :: u_init_obserr_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAF_gen_obs(step, dim_p, dim_obs_f, dim_ens, state_p, ainv, ens_p,  &
         u_init_dim_obs_f, u_obs_op_f, u_get_obs_f, u_init_obserr_f,  &
         u_prepoststep, screen, flag)

   END SUBROUTINE c__PDAF_gen_obs

   SUBROUTINE c__PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, state_p, ens_p,  &
      u_init_dim_obs, u_obs_op, u_init_obs, screen, debug, do_ens_mean,  &
      do_init_dim, do_hx, do_hxbar, do_init_obs) bind(c)
      use PDAFobs
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(inout) :: dim_obs_p
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Whether to compute ensemble mean
      LOGICAL(c_bool), INTENT(in) :: do_ens_mean
      ! Whether to call U_init_dim_obs
      LOGICAL(c_bool), INTENT(in) :: do_init_dim
      ! Whether to initialize HX_p
      LOGICAL(c_bool), INTENT(in) :: do_hx
      ! Whether to initialize HXbar
      LOGICAL(c_bool), INTENT(in) :: do_hxbar
      ! Whether to initialize obs_p
      LOGICAL(c_bool), INTENT(in) :: do_init_obs

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs

      logical :: do_ens_mean_in
      logical :: do_init_dim_in
      logical :: do_hx_in
      logical :: do_hxbar_in
      logical :: do_init_obs_in

      ! Convert logicals to C-compatible logicals
      do_ens_mean_in = do_ens_mean
      do_init_dim_in = do_init_dim
      do_hx_in = do_hx
      do_hxbar_in = do_hxbar
      do_init_obs_in = do_init_obs

      call PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, state_p, ens_p,  &
         u_init_dim_obs, u_obs_op, u_init_obs, screen, debug, do_ens_mean_in,  &
         do_init_dim_in, do_hx_in, do_hxbar_in, do_init_obs_in)

   END SUBROUTINE c__PDAFobs_init

   SUBROUTINE c__PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f,  &
      dim_ens, u_init_dim_obs_l, u_g2l_obs, u_init_obs_l, debug) bind(c)
      use PDAFobs
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Size of local observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_l
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug

      ! Init. dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l

      call PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f, dim_ens,  &
         u_init_dim_obs_l, u_g2l_obs, u_init_obs_l, debug)

   END SUBROUTINE c__PDAFobs_init_local

   SUBROUTINE c__PDAFobs_init_obsvars(step, dim_obs_p, u_init_obsvars) bind(c)
      use PDAFobs
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p

      ! Initialize vector of observation error variances
      procedure(c__init_obsvars_pdaf) :: u_init_obsvars

      call PDAFobs_init_obsvars(step, dim_obs_p, u_init_obsvars)

   END SUBROUTINE c__PDAFobs_init_obsvars

   SUBROUTINE c__PDAFobs_dealloc() bind(c)
      use PDAFobs
      implicit none
      call PDAFobs_dealloc()

   END SUBROUTINE c__PDAFobs_dealloc

   SUBROUTINE c__PDAFobs_dealloc_local() bind(c)
      use PDAFobs
      implicit none
      call PDAFobs_dealloc_local()

   END SUBROUTINE c__PDAFobs_dealloc_local

   SUBROUTINE c__PDAF_NETF_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_NETF
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_NETF_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_NETF_init

   SUBROUTINE c__PDAF_netf_alloc(outflag) bind(c)
      use PDAF_NETF
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_netf_alloc(outflag)

   END SUBROUTINE c__PDAF_netf_alloc

   SUBROUTINE c__PDAF_netf_config(subtype, verbose) bind(c)
      use PDAF_NETF
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_netf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_netf_config

   SUBROUTINE c__PDAF_netf_set_iparam(id, value, flag) bind(c)
      use PDAF_NETF
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_netf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_netf_set_iparam

   SUBROUTINE c__PDAF_netf_set_rparam(id, value, flag) bind(c)
      use PDAF_NETF
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_netf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_netf_set_rparam

   SUBROUTINE c__PDAF_netf_options() bind(c)
      use PDAF_NETF
      implicit none
      call PDAF_netf_options()

   END SUBROUTINE c__PDAF_netf_options

   SUBROUTINE c__PDAF_netf_memtime(printtype) bind(c)
      use PDAF_NETF
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_netf_memtime(printtype)

   END SUBROUTINE c__PDAF_netf_memtime

   SUBROUTINE c__PDAF_lenkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_lenkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
      ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out

   END SUBROUTINE c__PDAF_lenkf_init

   SUBROUTINE c__PDAF_lenkf_alloc(outflag) bind(c)
      use PDAF_lenkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_lenkf_alloc(outflag)

   END SUBROUTINE c__PDAF_lenkf_alloc

   SUBROUTINE c__PDAF_lenkf_config(subtype, verbose) bind(c)
      use PDAF_lenkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_lenkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_lenkf_config

   SUBROUTINE c__PDAF_lenkf_set_iparam(id, value, flag) bind(c)
      use PDAF_lenkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lenkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_lenkf_set_iparam

   SUBROUTINE c__PDAF_lenkf_set_rparam(id, value, flag) bind(c)
      use PDAF_lenkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lenkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_lenkf_set_rparam

   SUBROUTINE c__PDAF_lenkf_options() bind(c)
      use PDAF_lenkf
      implicit none
      call PDAF_lenkf_options()

   END SUBROUTINE c__PDAF_lenkf_options

   SUBROUTINE c__PDAF_lenkf_memtime(printtype) bind(c)
      use PDAF_lenkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_lenkf_memtime(printtype)

   END SUBROUTINE c__PDAF_lenkf_memtime

   SUBROUTINE c__PDAF_lseik_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_lseik
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_lseik_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_lseik_init

   SUBROUTINE c__PDAF_lseik_alloc(outflag) bind(c)
      use PDAF_lseik
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_lseik_alloc(outflag)

   END SUBROUTINE c__PDAF_lseik_alloc

   SUBROUTINE c__PDAF_lseik_config(subtype, verbose) bind(c)
      use PDAF_lseik
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_lseik_config(subtype, verbose)

   END SUBROUTINE c__PDAF_lseik_config

   SUBROUTINE c__PDAF_lseik_set_iparam(id, value, flag) bind(c)
      use PDAF_lseik
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lseik_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_lseik_set_iparam

   SUBROUTINE c__PDAF_lseik_set_rparam(id, value, flag) bind(c)
      use PDAF_lseik
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lseik_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_lseik_set_rparam

   SUBROUTINE c__PDAF_lseik_options() bind(c)
      use PDAF_lseik
      implicit none
      call PDAF_lseik_options()

   END SUBROUTINE c__PDAF_lseik_options

   SUBROUTINE c__PDAF_lseik_memtime(printtype) bind(c)
      use PDAF_lseik
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_lseik_memtime(printtype)

   END SUBROUTINE c__PDAF_lseik_memtime


   SUBROUTINE c__PDAF_etkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_etkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_etkf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
      ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_etkf_init

   SUBROUTINE c__PDAF_etkf_alloc(outflag) bind(c)
      use PDAF_etkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_etkf_alloc(outflag)

   END SUBROUTINE c__PDAF_etkf_alloc

   SUBROUTINE c__PDAF_etkf_config(subtype, verbose) bind(c)
      use PDAF_etkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_etkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_etkf_config

   SUBROUTINE c__PDAF_etkf_set_iparam(id, value, flag) bind(c)
      use PDAF_etkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_etkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_etkf_set_iparam

   SUBROUTINE c__PDAF_etkf_set_rparam(id, value, flag) bind(c)
      use PDAF_etkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_etkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_etkf_set_rparam

   SUBROUTINE c__PDAF_etkf_options() bind(c)
      use PDAF_etkf
      implicit none
      call PDAF_etkf_options()

   END SUBROUTINE c__PDAF_etkf_options

   SUBROUTINE c__PDAF_etkf_memtime(printtype) bind(c)
      use PDAF_etkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_etkf_memtime(printtype)

   END SUBROUTINE c__PDAF_etkf_memtime

   SUBROUTINE c__PDAFlenkf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ens_p, u_init_dim_obs, u_obs_op, u_add_obs_err, u_init_obs,  &
      u_init_obs_covar, u_prepoststep, u_localize, screen, subtype, flag) bind(c)
      use PDAF_lenkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Specification of filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: u_localize

      call PDAFlenkf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         u_init_dim_obs, u_obs_op, u_add_obs_err, u_init_obs, u_init_obs_covar,  &
         u_prepoststep, u_localize, screen, subtype, flag)

   END SUBROUTINE c__PDAFlenkf_update

   SUBROUTINE c__PDAF_PF_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_PF
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_PF_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
      ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_PF_init

   SUBROUTINE c__PDAF_pf_alloc(outflag) bind(c)
      use PDAF_pf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_pf_alloc(outflag)

   END SUBROUTINE c__PDAF_pf_alloc

   SUBROUTINE c__PDAF_pf_config(subtype, verbose) bind(c)
      use PDAF_pf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_pf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_pf_config

   SUBROUTINE c__PDAF_pf_set_iparam(id, value, flag) bind(c)
      use PDAF_pf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_pf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_pf_set_iparam

   SUBROUTINE c__PDAF_pf_set_rparam(id, value, flag) bind(c)
      use PDAF_pf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_pf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_pf_set_rparam

   SUBROUTINE c__PDAF_pf_options() bind(c)
      use PDAF_pf
      implicit none
      call PDAF_pf_options()

   END SUBROUTINE c__PDAF_pf_options

   SUBROUTINE c__PDAF_pf_memtime(printtype) bind(c)
      use PDAF_pf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_pf_memtime(printtype)

   END SUBROUTINE c__PDAF_pf_memtime

   SUBROUTINE c__PDAF_lknetf_ana_letkfT(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget,  &
      u_prodrinva_hyb_l, u_init_obsvar_l, gamma, screen, type_forget, flag) bind(c)
      use PDAF_lknetf_analysis_step
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! local forecast state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! local weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! PE-local full observed state ens.
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hz_l
      ! local observed ens. mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: rndmat
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Hybrid weight for state transformation
      REAL(c_double), DIMENSION(1), INTENT(inout) :: gamma
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of forgetting factor
      INTEGER(c_int), INTENT(in) :: type_forget
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain including hybrid weight
      procedure(c__prodrinva_hyb_l_pdaf) :: u_prodrinva_hyb_l
      ! Initialize local mean observation error variance
      procedure(c__init_obsvar_l_pdaf) :: u_init_obsvar_l

      call PDAF_lknetf_ana_letkfT(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget,  &
         u_prodrinva_hyb_l, u_init_obsvar_l, gamma, screen, type_forget, flag)

   END SUBROUTINE c__PDAF_lknetf_ana_letkfT

   SUBROUTINE c__PDAF_lknetf_ana_lnetf(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, ens_l, hx_l, rndmat, obs_l, u_likelihood_hyb_l, cnt_small_svals,  &
      n_eff_all, gamma, screen, flag) bind(c)
      use PDAF_lknetf_analysis_step
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! local observed state ens.
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: rndmat
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Number of small eigen values
      INTEGER(c_int), INTENT(inout) :: cnt_small_svals
      ! Effective ensemble size
      REAL(c_double), DIMENSION(1), INTENT(inout) :: n_eff_all
      ! Hybrid weight for state transformation
      REAL(c_double), DIMENSION(1), INTENT(inout) :: gamma
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Compute observation likelihood for an ensemble member with hybrid weight
      procedure(c__likelihood_hyb_l_pdaf) :: u_likelihood_hyb_l

      call PDAF_lknetf_ana_lnetf(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         ens_l, hx_l, rndmat, obs_l, u_likelihood_hyb_l, cnt_small_svals,  &
         n_eff_all, gamma, screen, flag)

   END SUBROUTINE c__PDAF_lknetf_ana_lnetf

   SUBROUTINE c__PDAF_enkf_ana_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
      state_p, ens_p, hx_p, hxbar_p, obs_p, u_add_obs_err, u_init_obs_covar,  &
      screen, debug, flag) bind(c)
      use PDAF_enkf_analysis_rsm
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank to be considered for inversion of HPH
      INTEGER(c_int), INTENT(in) :: rank_ana
      ! PE-local ensemble mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: hx_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar

      call PDAF_enkf_ana_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
         state_p, ens_p, hx_p, hxbar_p, obs_p, u_add_obs_err, u_init_obs_covar,  &
         screen, debug, flag)

   END SUBROUTINE c__PDAF_enkf_ana_rsm

   SUBROUTINE c__PDAF_lknetf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_lknetf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_lknetf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_lknetf_init

   SUBROUTINE c__PDAF_lknetf_alloc(outflag) bind(c)
      use PDAF_lknetf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_lknetf_alloc(outflag)

   END SUBROUTINE c__PDAF_lknetf_alloc

   SUBROUTINE c__PDAF_lknetf_config(subtype, verbose) bind(c)
      use PDAF_lknetf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_lknetf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_lknetf_config

   SUBROUTINE c__PDAF_lknetf_set_iparam(id, value, flag) bind(c)
      use PDAF_lknetf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lknetf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_lknetf_set_iparam

   SUBROUTINE c__PDAF_lknetf_set_rparam(id, value, flag) bind(c)
      use PDAF_lknetf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lknetf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_lknetf_set_rparam

   SUBROUTINE c__PDAF_lknetf_options() bind(c)
      use PDAF_lknetf
      implicit none
      call PDAF_lknetf_options()

   END SUBROUTINE c__PDAF_lknetf_options

   SUBROUTINE c__PDAF_lknetf_memtime(printtype) bind(c)
      use PDAF_lknetf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_lknetf_memtime(printtype)

   END SUBROUTINE c__PDAF_lknetf_memtime

   SUBROUTINE c__PDAF_lknetf_alpha_neff(dim_ens, weights, hlimit, alpha) bind(c)
      use PDAF_lknetf
      implicit none
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Weights
      REAL(c_double), DIMENSION(dim_ens), INTENT(in) :: weights
      ! Minimum of n_eff / N
      REAL(c_double), INTENT(in) :: hlimit
      ! hybrid weight
      REAL(c_double), INTENT(inout) :: alpha


      call PDAF_lknetf_alpha_neff(dim_ens, weights, hlimit, alpha)

   END SUBROUTINE c__PDAF_lknetf_alpha_neff

   SUBROUTINE c__PDAF_lknetf_compute_gamma(domain_p, step, dim_obs_l, dim_ens,  &
      hx_l, hxbar_l, obs_l, type_hyb, hyb_g, hyb_k, gamma, n_eff_out,  &
      skew_mabs, kurt_mabs, u_likelihood_l, screen, flag) bind(c)
      use PDAF_lknetf, only: PDAF_lknetf_compute_gamma
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! local observed state ens.
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
      ! local mean observed ensemble
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Type of hybrid weight
      INTEGER(c_int), INTENT(in) :: type_hyb
      ! Prescribed hybrid weight for state transformation
      REAL(c_double), INTENT(in) :: hyb_g
      ! Hybrid weight for covariance transformation
      REAL(c_double), INTENT(in) :: hyb_k
      ! Hybrid weight for state transformation
      REAL(c_double), DIMENSION(1), INTENT(inout) :: gamma
      ! Effective ensemble size
      REAL(c_double), DIMENSION(1), INTENT(inout) :: n_eff_out
      ! Mean absolute skewness
      REAL(c_double), DIMENSION(1), INTENT(inout) :: skew_mabs
      ! Mean absolute kurtosis
      REAL(c_double), DIMENSION(1), INTENT(inout) :: kurt_mabs
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l

      call PDAF_lknetf_compute_gamma(domain_p, step, dim_obs_l, dim_ens, hx_l,  &
         hxbar_l, obs_l, type_hyb, hyb_g, hyb_k, gamma, n_eff_out, skew_mabs,  &
         kurt_mabs, u_likelihood_l, screen, flag)

   END SUBROUTINE c__PDAF_lknetf_compute_gamma

   SUBROUTINE c__PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, hx_l,  &
      hxbar_l, weights, type_hyb, hyb_g, hyb_k, gamma, n_eff_out, maskew,  &
      makurt, screen, flag) bind(c)
      use PDAF_lknetf, only: PDAF_lknetf_set_gamma
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! local observed state ens.
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(in) :: hx_l
      ! local mean observed ensemble
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Weight vector
      REAL(c_double), DIMENSION(dim_ens), INTENT(in) :: weights
      ! Type of hybrid weight
      INTEGER(c_int), INTENT(in) :: type_hyb
      ! Prescribed hybrid weight for state transformation
      REAL(c_double), INTENT(in) :: hyb_g
      ! Scale factor kappa (for type_hyb 3 and 4)
      REAL(c_double), INTENT(in) :: hyb_k
      ! Hybrid weight for state transformation
      REAL(c_double), DIMENSION(1), INTENT(inout) :: gamma
      ! Effective ensemble size
      REAL(c_double), DIMENSION(1), INTENT(inout) :: n_eff_out
      ! Mean absolute skewness
      REAL(c_double), DIMENSION(1), INTENT(inout) :: maskew
      ! Mean absolute kurtosis
      REAL(c_double), DIMENSION(1), INTENT(inout) :: makurt
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, hx_l, hxbar_l,  &
         weights, type_hyb, hyb_g, hyb_k, gamma, n_eff_out, maskew, makurt,  &
         screen, flag)

   END SUBROUTINE c__PDAF_lknetf_set_gamma

   SUBROUTINE c__PDAF_lknetf_reset_gamma(gamma_in) bind(c)
      use PDAF_lknetf
      implicit none
      ! Prescribed hybrid weight
      REAL(c_double), INTENT(in) :: gamma_in


      call PDAF_lknetf_reset_gamma(gamma_in)

   END SUBROUTINE c__PDAF_lknetf_reset_gamma

   SUBROUTINE c__PDAFhyb3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec, dim_cvec_ens, beta_3dvar, state_p, ens_p, state_inc_p, hxbar_p,  &
      obs_p, u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens,  &
      u_obs_op_lin, u_obs_op_adj, screen, type_opt, debug, flag) bind(c)
      use PDAF_hyb3dvar_analysis_cvt
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! Hybrid weight for hybrid 3D-Var
      REAL(c_double), INTENT(in) :: beta_3dvar
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local state analysis increment
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_inc_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of minimizer for 3DVar
      INTEGER(c_int), INTENT(in) :: type_opt
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! Apply control vector transform matrix to control vector (parameterized)
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix (parameterized)
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Apply control vector transform matrix to control vector (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj

      call PDAFhyb3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens, dim_cvec,  &
         dim_cvec_ens, beta_3dvar, state_p, ens_p, state_inc_p, hxbar_p, obs_p,  &
         u_prodrinva, u_cvt, u_cvt_adj, u_cvt_ens, u_cvt_adj_ens, u_obs_op_lin,  &
         u_obs_op_adj, screen, type_opt, debug, flag)

   END SUBROUTINE c__PDAFhyb3dvar_analysis_cvt

   SUBROUTINE c__PDAF3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_cvec,  &
      state_p, hxbar_p, obs_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin,  &
      u_obs_op_adj, screen, type_opt, debug, flag) bind(c)
      use PDAF_3dvar_analysis_cvt
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of minimizer for 3DVar
      INTEGER(c_int), INTENT(in) :: type_opt
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

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

      call PDAF3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_cvec, state_p,  &
         hxbar_p, obs_p, u_prodrinva, u_cvt, u_cvt_adj, u_obs_op_lin,  &
         u_obs_op_adj, screen, type_opt, debug, flag)

   END SUBROUTINE c__PDAF3dvar_analysis_cvt

   SUBROUTINE c__PDAF_lestkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_lestkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out


      call PDAF_lestkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out

   END SUBROUTINE c__PDAF_lestkf_init

   SUBROUTINE c__PDAF_lestkf_alloc(outflag) bind(c)
      use PDAF_lestkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_lestkf_alloc(outflag)

   END SUBROUTINE c__PDAF_lestkf_alloc

   SUBROUTINE c__PDAF_lestkf_config(subtype, verbose) bind(c)
      use PDAF_lestkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_lestkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_lestkf_config

   SUBROUTINE c__PDAF_lestkf_set_iparam(id, value, flag) bind(c)
      use PDAF_lestkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lestkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_lestkf_set_iparam

   SUBROUTINE c__PDAF_lestkf_set_rparam(id, value, flag) bind(c)
      use PDAF_lestkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lestkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_lestkf_set_rparam

   SUBROUTINE c__PDAF_lestkf_options() bind(c)
      use PDAF_lestkf
      implicit none
      call PDAF_lestkf_options()

   END SUBROUTINE c__PDAF_lestkf_options

   SUBROUTINE c__PDAF_lestkf_memtime(printtype) bind(c)
      use PDAF_lestkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_lestkf_memtime(printtype)

   END SUBROUTINE c__PDAF_lestkf_memtime

   SUBROUTINE c__PDAF_seik_ana(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
      uinv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, debug, flag) bind(c)
      use PDAF_seik_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of eigenvalue matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble (perturbations)
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hl_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_seik_ana(step, dim_p, dim_obs_p, dim_ens, rank, state_p, uinv,  &
         ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, debug, flag)

   END SUBROUTINE c__PDAF_seik_ana

   SUBROUTINE c__PDAF_seik_resample(subtype, dim_p, dim_ens, rank, uinv,  &
      state_p, enst_p, type_sqrt, type_trans, nm1vsn, screen, flag) bind(c)
      use PDAF_seik_analysis
      implicit none
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local ensemble times T
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: enst_p
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Flag which normalization of P ist used in SEIK
      INTEGER(c_int), INTENT(in) :: nm1vsn
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_seik_resample(subtype, dim_p, dim_ens, rank, uinv, state_p,  &
         enst_p, type_sqrt, type_trans, nm1vsn, screen, flag)

   END SUBROUTINE c__PDAF_seik_resample

   SUBROUTINE c__PDAF_lseik_ana(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
      rank, state_l, uinv_l, ens_l, hl_l, hxbar_l, obs_l, forget,  &
      u_prodrinva_l, screen, debug, flag) bind(c)
      use PDAF_lseik_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! State on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(in) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hl_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_lseik_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, rank,  &
         state_l, uinv_l, ens_l, hl_l, hxbar_l, obs_l, forget, u_prodrinva_l,  &
         screen, debug, flag)

   END SUBROUTINE c__PDAF_lseik_ana

   SUBROUTINE c__PDAF_lseik_resample(domain_p, subtype, dim_l, dim_ens, rank,  &
      uinv_l, state_l, ens_l, omegat_in, type_sqrt, screen, flag) bind(c)
      use PDAF_lseik_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Specification of filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv_l
      ! Local model state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Matrix Omega
      REAL(c_double), DIMENSION(rank, dim_ens), INTENT(inout) :: omegat_in
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_lseik_resample(domain_p, subtype, dim_l, dim_ens, rank, uinv_l,  &
         state_l, ens_l, omegat_in, type_sqrt, screen, flag)

   END SUBROUTINE c__PDAF_lseik_resample

   SUBROUTINE c__PDAF_prepost(u_collect_state, u_distribute_state,  &
      u_prepoststep, u_next_observation, outflag) bind(c)
      use PDAFprepost
      implicit none
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

      call PDAF_prepost(u_collect_state, u_distribute_state, u_prepoststep,  &
         u_next_observation, outflag)

   END SUBROUTINE c__PDAF_prepost

   SUBROUTINE c__PDAFenkf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ens_p, u_init_dim_obs, u_obs_op, u_add_obs_err, u_init_obs,  &
      u_init_obs_covar, u_prepoststep, screen, subtype, dim_lag, sens_p,  &
      cnt_maxlag, flag) bind(c)
      use PDAF_enkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Specification of filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFenkf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         u_init_dim_obs, u_obs_op, u_add_obs_err, u_init_obs, u_init_obs_covar,  &
         u_prepoststep, screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

   END SUBROUTINE c__PDAFenkf_update

   SUBROUTINE c__PDAF_init_parallel(dim_ens, ensemblefilter, fixedbasis,  &
      comm_model, in_comm_filter, in_comm_couple, in_n_modeltasks, in_task_id,  &
      screen, flag) bind(c)
      use PDAF_mod_parallel
      implicit none
      ! Rank of covar matrix/ensemble size
      INTEGER(c_int), INTENT(inout) :: dim_ens
      ! Is the filter ensemble-based?
      LOGICAL(c_bool), INTENT(in) :: ensemblefilter
      ! Run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(in) :: fixedbasis
      ! Model communicator (not shared)
      INTEGER(c_int), INTENT(in) :: comm_model
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: in_comm_filter
      ! Coupling communicator
      INTEGER(c_int), INTENT(in) :: in_comm_couple
      ! Number of model tasks
      INTEGER(c_int), INTENT(in) :: in_n_modeltasks
      ! Task ID of current PE
      INTEGER(c_int), INTENT(in) :: in_task_id
      ! Whether screen information is shown
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      logical :: ensemblefilter_in, fixedbasis_in

      ensemblefilter_in = ensemblefilter
      fixedbasis_in = fixedbasis
      call PDAF_init_parallel(dim_ens, ensemblefilter_in, fixedbasis_in, comm_model,  &
         in_comm_filter, in_comm_couple, in_n_modeltasks, in_task_id, screen, flag)

   END SUBROUTINE c__PDAF_init_parallel

   SUBROUTINE c__PDAF_seik_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_seik
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_seik_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_seik_init

   SUBROUTINE c__PDAF_seik_alloc(outflag) bind(c)
      use PDAF_seik
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_seik_alloc(outflag)

   END SUBROUTINE c__PDAF_seik_alloc

   SUBROUTINE c__PDAF_seik_config(subtype, verbose) bind(c)
      use PDAF_seik
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_seik_config(subtype, verbose)

   END SUBROUTINE c__PDAF_seik_config

   SUBROUTINE c__PDAF_seik_set_iparam(id, value, flag) bind(c)
      use PDAF_seik
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_seik_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_seik_set_iparam

   SUBROUTINE c__PDAF_seik_set_rparam(id, value, flag) bind(c)
      use PDAF_seik
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_seik_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_seik_set_rparam

   SUBROUTINE c__PDAF_seik_options() bind(c)
      use PDAF_seik
      implicit none
      call PDAF_seik_options()

   END SUBROUTINE c__PDAF_seik_options

   SUBROUTINE c__PDAF_seik_memtime(printtype) bind(c)
      use PDAF_seik
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_seik_memtime(printtype)

   END SUBROUTINE c__PDAF_seik_memtime

   SUBROUTINE c__PDAFnetf_update(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood,  &
      u_prepoststep, screen, subtype, dim_lag, sens_p, cnt_maxlag, flag) bind(c)
      use PDAF_netf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFnetf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
         ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood,  &
         u_prepoststep, screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

   END SUBROUTINE c__PDAFnetf_update

   SUBROUTINE c__PDAF_seik_ana_newT(step, dim_p, dim_obs_p, dim_ens, rank,  &
      state_p, uinv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
      debug, flag) bind(c)
      use PDAF_seik_analysis_newT
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of eigenvalue matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble (perturbations)
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hl_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_seik_ana_newT(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
         uinv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen, debug,  &
         flag)

   END SUBROUTINE c__PDAF_seik_ana_newT

   SUBROUTINE c__PDAF_seik_resample_newT(subtype, dim_p, dim_ens, rank, uinv,  &
      state_p, ens_p, type_sqrt, type_trans, nm1vsn, screen, flag) bind(c)
      use PDAF_seik_analysis_newT
      implicit none
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Flag which normalization of P ist used in SEIK
      INTEGER(c_int), INTENT(in) :: nm1vsn
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_seik_resample_newT(subtype, dim_p, dim_ens, rank, uinv,  &
         state_p, ens_p, type_sqrt, type_trans, nm1vsn, screen, flag)

   END SUBROUTINE c__PDAF_seik_resample_newT

   SUBROUTINE c__PDAF_lenkf_ana_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
      state_p, ens_p, hx_p, hxbar_p, obs_p, u_add_obs_err, u_init_obs_covar,  &
      u_localize, screen, debug, flag) bind(c)
      use PDAF_lenkf_analysis_rsm
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank to be considered for inversion of HPH
      INTEGER(c_int), INTENT(in) :: rank_ana
      ! PE-local ensemble mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: hx_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Add observation error covariance matrix
      procedure(c__add_obs_err_pdaf) :: u_add_obs_err
      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: u_localize

      call PDAF_lenkf_ana_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana,  &
         state_p, ens_p, hx_p, hxbar_p, obs_p, u_add_obs_err, u_init_obs_covar,  &
         u_localize, screen, debug, flag)

   END SUBROUTINE c__PDAF_lenkf_ana_rsm

   SUBROUTINE c__PDAF_lestkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
      rank, state_l, ainv_l, ens_l, hl_l, hxbar_l, obs_l, omegat_in, forget,  &
      u_prodrinva_l, envar_mode, type_sqrt, ta, screen, debug, flag) bind(c)
      use PDAF_lestkf_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! state on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! Inverse of matrix U - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hl_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Matrix Omega
      REAL(c_double), DIMENSION(rank, dim_ens), INTENT(in) :: omegat_in
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Flag whether routine is called from 3DVar for special functionality
      INTEGER(c_int), INTENT(in) :: envar_mode
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Ensemble transformation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ta
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_lestkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, rank,  &
         state_l, ainv_l, ens_l, hl_l, hxbar_l, obs_l, omegat_in, forget,  &
         u_prodrinva_l, envar_mode, type_sqrt, ta, screen, debug, flag)

   END SUBROUTINE c__PDAF_lestkf_ana

   SUBROUTINE c__PDAFlestkf_update(step, dim_p, dim_obs_f, dim_ens, rank,  &
      state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_prepoststep, screen, subtype, envar_mode, dim_lag, sens_p, cnt_maxlag,  &
      flag) bind(c)
      use PDAF_lestkf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Flag whether routine is called from 3DVar for special functionality
      INTEGER(c_int), INTENT(in) :: envar_mode
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute product of R^(-1) with HV
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from global state
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

      call PDAFlestkf_update(step, dim_p, dim_obs_f, dim_ens, rank, state_p,  &
         ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_prepoststep, screen, subtype, envar_mode, dim_lag, sens_p,  &
         cnt_maxlag, flag)

   END SUBROUTINE c__PDAFlestkf_update

   SUBROUTINE c__PDAF_LNETF_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_LNETF
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_LNETF_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_LNETF_init

   SUBROUTINE c__PDAF_lnetf_alloc(outflag) bind(c)
      use PDAF_LNETF
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_lnetf_alloc(outflag)

   END SUBROUTINE c__PDAF_lnetf_alloc

   SUBROUTINE c__PDAF_lnetf_config(subtype, verbose) bind(c)
      use PDAF_LNETF
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_lnetf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_lnetf_config

   SUBROUTINE c__PDAF_lnetf_set_iparam(id, value, flag) bind(c)
      use PDAF_LNETF
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lnetf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_lnetf_set_iparam

   SUBROUTINE c__PDAF_lnetf_set_rparam(id, value, flag) bind(c)
      use PDAF_LNETF
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_lnetf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_lnetf_set_rparam

   SUBROUTINE c__PDAF_lnetf_options() bind(c)
      use PDAF_LNETF
      implicit none
      call PDAF_lnetf_options()

   END SUBROUTINE c__PDAF_lnetf_options

   SUBROUTINE c__PDAF_lnetf_memtime(printtype) bind(c)
      use PDAF_LNETF
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_lnetf_memtime(printtype)

   END SUBROUTINE c__PDAF_lnetf_memtime

   SUBROUTINE c__PDAF_enkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_enkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out
      call PDAF_enkf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_enkf_init

   SUBROUTINE c__PDAF_enkf_alloc(outflag) bind(c)
      use PDAF_enkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_enkf_alloc(outflag)

   END SUBROUTINE c__PDAF_enkf_alloc

   SUBROUTINE c__PDAF_enkf_config(subtype, verbose) bind(c)
      use PDAF_enkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_enkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_enkf_config

   SUBROUTINE c__PDAF_enkf_set_iparam(id, value, flag) bind(c)
      use PDAF_enkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_enkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_enkf_set_iparam

   SUBROUTINE c__PDAF_enkf_set_rparam(id, value, flag) bind(c)
      use PDAF_enkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_enkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_enkf_set_rparam

   SUBROUTINE c__PDAF_enkf_options() bind(c)
      use PDAF_enkf
      implicit none
      call PDAF_enkf_options()

   END SUBROUTINE c__PDAF_enkf_options

   SUBROUTINE c__PDAF_enkf_memtime(printtype) bind(c)
      use PDAF_enkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_enkf_memtime(printtype)

   END SUBROUTINE c__PDAF_enkf_memtime

   SUBROUTINE c__PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p,  &
      resid) bind(c)
      use PDAF_enkf
      implicit none
      ! Global observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local residual matrix
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: resid_p
      ! Global residual matrix
      REAL(c_double), DIMENSION(dim_obs, dim_ens), INTENT(out) :: resid


      call PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)

   END SUBROUTINE c__PDAF_enkf_gather_resid

   SUBROUTINE c__PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens,  &
      obsens_p, obs_p, u_init_obs_covar, screen, flag) bind(c)
      use PDAF_enkf
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Local dimension of current observation
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local obs. ensemble
      REAL(c_double), DIMENSION(dim_obs_p,dim_ens), INTENT(out) :: obsens_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize observation error covariance matrix
      procedure(c__init_obs_covar_pdaf) :: u_init_obs_covar

      call PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens, obsens_p,  &
         obs_p, u_init_obs_covar, screen, flag)

   END SUBROUTINE c__PDAF_enkf_obs_ensemble

   SUBROUTINE c__PDAFpf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
      ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood, u_prepoststep,  &
      screen, subtype, flag) bind(c)
      use PDAF_pf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_pdaf) :: u_likelihood
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFpf_update(step, dim_p, dim_obs_p, dim_ens, state_p, ainv, ens_p,  &
         u_init_dim_obs, u_obs_op, u_init_obs, u_likelihood, u_prepoststep,  &
         screen, subtype, flag)

   END SUBROUTINE c__PDAFpf_update

   SUBROUTINE c__PDAF_generate_rndmat(dim, rndmat, mattype) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Size of matrix rndmat
      INTEGER(c_int), INTENT(in) :: dim
      ! Matrix
      REAL(c_double), DIMENSION(dim, dim), INTENT(out) :: rndmat
      ! Select type of random matrix:
      INTEGER(c_int), INTENT(in) :: mattype


      call PDAF_generate_rndmat(dim, rndmat, mattype)

   END SUBROUTINE c__PDAF_generate_rndmat

   SUBROUTINE c__PDAF_print_domain_stats(n_domains_p) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Number of PE-local analysis domains
      INTEGER(c_int), INTENT(in) :: n_domains_p


      call PDAF_print_domain_stats(n_domains_p)

   END SUBROUTINE c__PDAF_print_domain_stats

   SUBROUTINE c__PDAF_init_local_obsstats() bind(c)
      use PDAF_analysis_utils
      implicit none
      call PDAF_init_local_obsstats()

   END SUBROUTINE c__PDAF_init_local_obsstats

   SUBROUTINE c__PDAF_incr_local_obsstats(dim_obs_l) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Number of locally assimilated observations
      INTEGER(c_int), INTENT(in) :: dim_obs_l


      call PDAF_incr_local_obsstats(dim_obs_l)

   END SUBROUTINE c__PDAF_incr_local_obsstats

   SUBROUTINE c__PDAF_print_local_obsstats(screen, n_domains_with_obs) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      !
      INTEGER(c_int), INTENT(out) :: n_domains_with_obs


      call PDAF_print_local_obsstats(screen, n_domains_with_obs)

   END SUBROUTINE c__PDAF_print_local_obsstats

   SUBROUTINE c__PDAF_seik_matrixT(dim, dim_ens, a) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! dimension of states
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Input/output matrix
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(inout) :: a


      call PDAF_seik_matrixT(dim, dim_ens, a)

   END SUBROUTINE c__PDAF_seik_matrixT

   SUBROUTINE c__PDAF_seik_TtimesA(rank, dim_col, a, b) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim_col
      ! Input matrix
      REAL(c_double), DIMENSION(rank, dim_col), INTENT(in) :: a
      ! Output matrix (TA)
      REAL(c_double), DIMENSION(rank+1, dim_col), INTENT(out) :: b


      call PDAF_seik_TtimesA(rank, dim_col, a, b)

   END SUBROUTINE c__PDAF_seik_TtimesA

   SUBROUTINE c__PDAF_seik_Omega(rank, omega, omegatype, screen) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Approximated rank of covar matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Matrix Omega
      REAL(c_double), DIMENSION(rank+1, rank), INTENT(inout) :: omega
      ! Select type of Omega:
      INTEGER(c_int), INTENT(in) :: omegatype
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_seik_Omega(rank, omega, omegatype, screen)

   END SUBROUTINE c__PDAF_seik_Omega

   SUBROUTINE c__PDAF_seik_Uinv(rank, uinv) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv


      call PDAF_seik_Uinv(rank, uinv)

   END SUBROUTINE c__PDAF_seik_Uinv

   SUBROUTINE c__PDAF_ens_Omega(seed, r, dim_ens, omega, norm, otype,  &
      screen) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Seed for random number generation
      INTEGER(c_int), DIMENSION(4), INTENT(in) :: seed
      ! Approximated rank of covar matrix
      INTEGER(c_int), INTENT(in) :: r
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Random matrix
      REAL(c_double), DIMENSION(dim_ens,r), INTENT(inout) :: omega
      ! Norm for ensemble transformation
      REAL(c_double), INTENT(inout) :: norm
      ! Type of Omega:
      INTEGER(c_int), INTENT(in) :: otype
      ! Control verbosity
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_ens_Omega(seed, r, dim_ens, omega, norm, otype, screen)

   END SUBROUTINE c__PDAF_ens_Omega

   SUBROUTINE c__PDAF_estkf_OmegaA(rank, dim_col, a, b) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim_col
      ! Input matrix
      REAL(c_double), DIMENSION(rank, dim_col), INTENT(in) :: a
      ! Output matrix (TA)
      REAL(c_double), DIMENSION(rank+1, dim_col), INTENT(out) :: b


      call PDAF_estkf_OmegaA(rank, dim_col, a, b)

   END SUBROUTINE c__PDAF_estkf_OmegaA

   SUBROUTINE c__PDAF_estkf_AOmega(dim, dim_ens, a) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! dimension of states
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Input/output matrix
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(inout) :: a


      call PDAF_estkf_AOmega(dim, dim_ens, a)

   END SUBROUTINE c__PDAF_estkf_AOmega

   SUBROUTINE c__PDAF_subtract_rowmean(dim, dim_ens, a) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! dimension of states
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Input/output matrix
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(inout) :: a


      call PDAF_subtract_rowmean(dim, dim_ens, a)

   END SUBROUTINE c__PDAF_subtract_rowmean

   SUBROUTINE c__PDAF_subtract_colmean(dim_ens, dim, a) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of columns in A and B
      INTEGER(c_int), INTENT(in) :: dim
      ! Input/output matrix
      REAL(c_double), DIMENSION(dim_ens, dim), INTENT(inout) :: a


      call PDAF_subtract_colmean(dim_ens, dim, a)

   END SUBROUTINE c__PDAF_subtract_colmean

   SUBROUTINE c__PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p,  &
      type_noise, noise_amp, screen) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Number of particles
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector (not filled)
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Ensemble array
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Type of noise
      INTEGER(c_int), INTENT(in) :: type_noise
      ! Noise amplitude
      REAL(c_double), INTENT(in) :: noise_amp
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p, type_noise,  &
         noise_amp, screen)

   END SUBROUTINE c__PDAF_add_particle_noise

   SUBROUTINE c__PDAF_inflate_weights(screen, dim_ens, alpha, weights) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Minimum limit of n_eff / N
      REAL(c_double), INTENT(in) :: alpha
      ! weights (before and after inflation)
      REAL(c_double), DIMENSION(dim_ens), INTENT(inout) :: weights


      call PDAF_inflate_weights(screen, dim_ens, alpha, weights)

   END SUBROUTINE c__PDAF_inflate_weights

   SUBROUTINE c__PDAF_inflate_ens(dim, dim_ens, meanstate, ens, forget,  &
      do_ensmean) bind(c)
      use PDAF_analysis_utils
      implicit none
      ! dimension of states
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! state vector to hold ensemble mean
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: meanstate
      ! Input/output ensemble matrix
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(inout) :: ens
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Whether to compute the ensemble mean state
      LOGICAL(c_bool), INTENT(in) :: do_ensmean

      logical :: do_ensmean_in
      do_ensmean_in = do_ensmean
      call PDAF_inflate_ens(dim, dim_ens, meanstate, ens, forget, do_ensmean_in)

   END SUBROUTINE c__PDAF_inflate_ens

   SUBROUTINE c__PDAF_alloc(dim_p, dim_ens, dim_ens_task, dim_es, dim_bias_p,  &
      dim_lag, statetask, outflag) bind(c)
      use pdaf_utils
      implicit none
      ! Size of state vector
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Ensemble size handled by a model task
      INTEGER(c_int), INTENT(in) :: dim_ens_task
      ! Dimension of error space (size of Ainv)
      INTEGER(c_int), INTENT(in) :: dim_es
      ! Size of bias vector
      INTEGER(c_int), INTENT(in) :: dim_bias_p
      ! Smoother lag
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Task ID forecasting a single state
      INTEGER(c_int), INTENT(in) :: statetask
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_alloc(dim_p, dim_ens, dim_ens_task, dim_es, dim_bias_p,  &
         dim_lag, statetask, outflag)

   END SUBROUTINE c__PDAF_alloc

   SUBROUTINE c__PDAF_smoothing(dim_p, dim_ens, dim_lag, ainv, sens_p,  &
      cnt_maxlag, forget, screen) bind(c)
      use PDAF_smoother
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: ainv
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_smoothing(dim_p, dim_ens, dim_lag, ainv, sens_p, cnt_maxlag,  &
         forget, screen)

   END SUBROUTINE c__PDAF_smoothing

   SUBROUTINE c__PDAF_smoothing_local(domain_p, step, dim_p, dim_l, dim_ens,  &
      dim_lag, ainv, ens_l, sens_p, cnt_maxlag, u_g2l_state, u_l2g_state,  &
      forget, screen) bind(c)
      use PDAF_smoother
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! Weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(in) :: ainv
      ! local past ensemble (temporary)
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Get state on local ana. domain from global state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state

      call PDAF_smoothing_local(domain_p, step, dim_p, dim_l, dim_ens, dim_lag,  &
         ainv, ens_l, sens_p, cnt_maxlag, u_g2l_state, u_l2g_state, forget, screen)

   END SUBROUTINE c__PDAF_smoothing_local

   SUBROUTINE c__PDAF_smoother_shift(dim_p, dim_ens, dim_lag, ens_p, sens_p,  &
      cnt_maxlag, screen) bind(c)
      use PDAF_smoother
      implicit none
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Number of past time instances for smoother
      INTEGER(c_int), INTENT(in) :: dim_lag
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, 1), INTENT(inout) :: ens_p
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count available number of time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_smoother_shift(dim_p, dim_ens, dim_lag, ens_p, sens_p,  &
         cnt_maxlag, screen)

   END SUBROUTINE c__PDAF_smoother_shift

   SUBROUTINE c__PDAFlknetf_update_sync(step, dim_p, dim_obs_f, dim_ens,  &
      state_p, ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
      u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
      u_likelihood_l, u_prepoststep, screen, subtype, flag) bind(c)
      use PDAF_lknetf_update_sync
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute product of R^(-1) with HV
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from global state
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
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFlknetf_update_sync(step, dim_p, dim_obs_f, dim_ens, state_p,  &
         ainv, ens_p, u_init_dim_obs, u_obs_op, u_init_obs, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         u_likelihood_l, u_prepoststep, screen, subtype, flag)

   END SUBROUTINE c__PDAFlknetf_update_sync

   SUBROUTINE c__PDAF_etkf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
      ens_p, hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen, type_trans,  &
      debug, flag) bind(c)
      use PDAF_etkf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! on exit: weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hz_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_etkf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ainv, ens_p,  &
         hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen, type_trans, debug,  &
         flag)

   END SUBROUTINE c__PDAF_etkf_ana

   SUBROUTINE c__PDAF_letkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
      state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget,  &
      u_prodrinva_l, type_trans, screen, debug, flag) bind(c)
      use PDAF_letkf_analysis
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Local forecast state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! on exit: local weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hz_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Global random rotation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: rndmat
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_letkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, state_l,  &
         ainv_l, ens_l, hz_l, hxbar_l, obs_l, rndmat, forget, u_prodrinva_l,  &
         type_trans, screen, debug, flag)

   END SUBROUTINE c__PDAF_letkf_ana

   SUBROUTINE c__PDAF_letkf_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_letkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_letkf_init(subtype, param_int, dim_pint, param_real, dim_preal,  &
         ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_letkf_init

   SUBROUTINE c__PDAF_letkf_alloc(outflag) bind(c)
      use PDAF_letkf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_letkf_alloc(outflag)

   END SUBROUTINE c__PDAF_letkf_alloc

   SUBROUTINE c__PDAF_letkf_config(subtype, verbose) bind(c)
      use PDAF_letkf
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_letkf_config(subtype, verbose)

   END SUBROUTINE c__PDAF_letkf_config

   SUBROUTINE c__PDAF_letkf_set_iparam(id, value, flag) bind(c)
      use PDAF_letkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_letkf_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_letkf_set_iparam

   SUBROUTINE c__PDAF_letkf_set_rparam(id, value, flag) bind(c)
      use PDAF_letkf
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_letkf_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_letkf_set_rparam

   SUBROUTINE c__PDAF_letkf_options() bind(c)
      use PDAF_letkf
      implicit none
      call PDAF_letkf_options()

   END SUBROUTINE c__PDAF_letkf_options

   SUBROUTINE c__PDAF_letkf_memtime(printtype) bind(c)
      use PDAF_letkf
      implicit none
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_letkf_memtime(printtype)

   END SUBROUTINE c__PDAF_letkf_memtime

   SUBROUTINE c__PDAF_estkf_ana(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
      ainv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
      envar_mode, type_sqrt, type_trans, ta, debug, flag) bind(c)
      use PDAF_estkf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(inout) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! on exit: PE-local forecast mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix A - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: ainv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hl_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag whether routine is called from 3DVar for special functionality
      INTEGER(c_int), INTENT(in) :: envar_mode
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Ensemble transformation matrix
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ta
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 with some matrix
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_estkf_ana(step, dim_p, dim_obs_p, dim_ens, rank, state_p, ainv,  &
         ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen, envar_mode,  &
         type_sqrt, type_trans, ta, debug, flag)

   END SUBROUTINE c__PDAF_estkf_ana

   SUBROUTINE c__PDAF_ensrf_ana(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ens_p, hx_p, hxbar_p, obs_p, var_obs_p, u_localize_covar_serial, screen,  &
      debug) bind(c)
      use PDAF_ensrf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local ensemble mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hx_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(inout) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! PE-local vector of observation eror variances
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: var_obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug

      ! Apply localization for single-observation vectors
      procedure(c__localize_covar_serial_pdaf) :: u_localize_covar_serial

      call PDAF_ensrf_ana(step, dim_p, dim_obs_p, dim_ens, state_p, ens_p,  &
         hx_p, hxbar_p, obs_p, var_obs_p, u_localize_covar_serial, screen, debug)

   END SUBROUTINE c__PDAF_ensrf_ana

   SUBROUTINE c__PDAF_ensrf_ana_2step(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ens_p, hx_p, hxbar_p, obs_p, var_obs_p, u_localize_covar_serial, screen,  &
      debug) bind(c)
      use PDAF_ensrf_analysis
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of state ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local ensemble mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hx_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(inout) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! PE-local vector of observation eror variances
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: var_obs_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug

      ! Apply localization for single-observation vectors
      procedure(c__localize_covar_serial_pdaf) :: u_localize_covar_serial

      call PDAF_ensrf_ana_2step(step, dim_p, dim_obs_p, dim_ens, state_p,  &
         ens_p, hx_p, hxbar_p, obs_p, var_obs_p, u_localize_covar_serial,  &
         screen, debug)

   END SUBROUTINE c__PDAF_ensrf_ana_2step

   SUBROUTINE c__PDAFlnetf_update(step, dim_p, dim_obs_f, dim_ens, state_p,  &
      ainv, ens_p, u_obs_op, u_init_dim_obs, u_init_obs, u_init_obs_l,  &
      u_likelihood_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
      u_g2l_state, u_l2g_state, u_g2l_obs, u_prepoststep, screen, subtype,  &
      dim_lag, sens_p, cnt_maxlag, flag) bind(c)
      use PDAF_lnetf_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_f
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: dim_lag
      ! PE-local smoother ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens, dim_lag), INTENT(inout) :: sens_p
      ! Count number of past time steps for smoothing
      INTEGER(c_int), INTENT(inout) :: cnt_maxlag
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Initialize PE-local observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Init. observation vector on local analysis domain
      procedure(c__init_obs_l_pdaf) :: u_init_obs_l
      ! Compute observation likelihood for an ensemble member
      procedure(c__likelihood_l_pdaf) :: u_likelihood_l
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: u_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: u_init_dim_l
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: u_init_dim_obs_l
      ! Get state on local ana. domain from global state
      procedure(c__g2l_state_pdaf) :: u_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: u_l2g_state
      ! Restrict full obs. vector to local analysis domain
      procedure(c__g2l_obs_pdaf) :: u_g2l_obs
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      call PDAFlnetf_update(step, dim_p, dim_obs_f, dim_ens, state_p, ainv,  &
         ens_p, u_obs_op, u_init_dim_obs, u_init_obs, u_init_obs_l,  &
         u_likelihood_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_prepoststep, screen, subtype,  &
         dim_lag, sens_p, cnt_maxlag, flag)

   END SUBROUTINE c__PDAFlnetf_update

   SUBROUTINE c__PDAF_seik_ana_trans(step, dim_p, dim_obs_p, dim_ens, rank,  &
      state_p, uinv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
      type_sqrt, type_trans, nm1vsn, debug, flag) bind(c)
      use PDAF_seik_analysis_trans
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! PE-local forecast mean state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Inverse of matrix U - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: uinv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble (perturbations)
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hl_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Type of normalization in covariance matrix computation
      INTEGER(c_int), INTENT(in) :: nm1vsn
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_seik_ana_trans(step, dim_p, dim_obs_p, dim_ens, rank, state_p,  &
         uinv, ens_p, hl_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
         type_sqrt, type_trans, nm1vsn, debug, flag)

   END SUBROUTINE c__PDAF_seik_ana_trans

   SUBROUTINE c__PDAFhyb3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec, dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens, u_cvt,  &
      u_cvt_adj, u_obs_op_lin, u_obs_op_adj, u_init_obsvar, screen, subtype,  &
      flag) bind(c)
      use PDAF_hyb3dvar_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Transform matrix
      REAL(c_double), DIMENSION(dim_ens-1, dim_ens-1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for 3DVAR analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Apply control vector transform matrix
      procedure(c__cvt_pdaf) :: u_cvt
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: u_cvt_adj
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: u_obs_op_lin
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: u_obs_op_adj
      ! Initialize mean observation error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar

      call PDAFhyb3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens, dim_cvec,  &
         dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op,  &
         u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens,  &
         u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj, u_init_obsvar, screen,  &
         subtype, flag)

   END SUBROUTINE c__PDAFhyb3dvar_update_estkf

   SUBROUTINE c__PDAFhyb3dvar_update_lestkf(step, dim_p, dim_obs_p, dim_ens,  &
      dim_cvec, dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs, u_obs_op,  &
      u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens, u_cvt_adj_ens, u_cvt,  &
      u_cvt_adj, u_obs_op_lin, u_obs_op_adj, u_init_dim_obs_f, u_obs_op_f,  &
      u_init_obs_f, u_init_obs_l, u_prodrinva_l, u_init_n_domains_p,  &
      u_init_dim_l, u_init_dim_obs_l, u_g2l_state, u_l2g_state, u_g2l_obs,  &
      u_init_obsvar, u_init_obsvar_l, screen, subtype, flag) bind(c)
      use PDAF_hyb3dvar_update
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(out) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Size of control vector (parameterized part)
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! Size of control vector (ensemble part)
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! Transform matrix for LESKTF
      REAL(c_double), DIMENSION(dim_ens-1, dim_ens-1), INTENT(inout) :: ainv
      ! PE-local ensemble matrix
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Filter subtype
      INTEGER(c_int), INTENT(in) :: subtype
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: u_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: u_obs_op
      ! Initialize observation vector
      procedure(c__init_obs_pdaf) :: u_init_obs
      ! Provide product R^-1 A for 3DVAR analysis
      procedure(c__prodrinva_pdaf) :: u_prodrinva
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep
      ! Apply control vector transform matrix (ensemble)
      procedure(c__cvt_ens_pdaf) :: u_cvt_ens
      ! Apply adjoint control vector transform matrix (ensemble var)
      procedure(c__cvt_adj_ens_pdaf) :: u_cvt_adj_ens
      ! Apply control vector transform matrix
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

      call PDAFhyb3dvar_update_lestkf(step, dim_p, dim_obs_p, dim_ens,  &
         dim_cvec, dim_cvec_ens, state_p, ainv, ens_p, u_init_dim_obs,  &
         u_obs_op, u_init_obs, u_prodrinva, u_prepoststep, u_cvt_ens,  &
         u_cvt_adj_ens, u_cvt, u_cvt_adj, u_obs_op_lin, u_obs_op_adj,  &
         u_init_dim_obs_f, u_obs_op_f, u_init_obs_f, u_init_obs_l,  &
         u_prodrinva_l, u_init_n_domains_p, u_init_dim_l, u_init_dim_obs_l,  &
         u_g2l_state, u_l2g_state, u_g2l_obs, u_init_obsvar, u_init_obsvar_l,  &
         screen, subtype, flag)

   END SUBROUTINE c__PDAFhyb3dvar_update_lestkf

   SUBROUTINE c__PDAF_lestkf_ana_fixed(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, rank, state_l, ainv_l, ens_l, hl_l, hxbar_l, obs_l, forget,  &
      u_prodrinva_l, type_sqrt, screen, debug, flag) bind(c)
      use PDAF_lestkf_analysis_fixed
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! state on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! Inverse of matrix U - temporary use only
      REAL(c_double), DIMENSION(rank, rank), INTENT(inout) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hl_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Type of square-root of A
      INTEGER(c_int), INTENT(in) :: type_sqrt
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_lestkf_ana_fixed(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         rank, state_l, ainv_l, ens_l, hl_l, hxbar_l, obs_l, forget,  &
         u_prodrinva_l, type_sqrt, screen, debug, flag)

   END SUBROUTINE c__PDAF_lestkf_ana_fixed

   SUBROUTINE c__PDAF_genobs_init(subtype, param_int, dim_pint, param_real,  &
      dim_preal, ensemblefilter, fixedbasis, verbose, outflag) bind(c)
      use PDAF_genobs
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameters
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Is the chosen filter ensemble-based?
      LOGICAL(c_bool), INTENT(out) :: ensemblefilter
      ! Does the filter run with fixed error-space basis?
      LOGICAL(c_bool), INTENT(out) :: fixedbasis
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      logical :: ensemblefilter_out, fixedbasis_out

      call PDAF_genobs_init(subtype, param_int, dim_pint, param_real,  &
         dim_preal, ensemblefilter_out, fixedbasis_out, verbose, outflag)
      ensemblefilter = ensemblefilter_out
      fixedbasis = fixedbasis_out
   END SUBROUTINE c__PDAF_genobs_init

   SUBROUTINE c__PDAF_genobs_alloc(outflag) bind(c)
      use PDAF_genobs
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag


      call PDAF_genobs_alloc(outflag)

   END SUBROUTINE c__PDAF_genobs_alloc

   SUBROUTINE c__PDAF_genobs_config(subtype, verbose) bind(c)
      use PDAF_genobs
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_genobs_config(subtype, verbose)

   END SUBROUTINE c__PDAF_genobs_config

   SUBROUTINE c__PDAF_genobs_set_iparam(id, value, flag) bind(c)
      use PDAF_genobs
      implicit none
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_genobs_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_genobs_set_iparam

   SUBROUTINE c__PDAF_genobs_options() bind(c)
      use PDAF_genobs
      implicit none
      call PDAF_genobs_options()

   END SUBROUTINE c__PDAF_genobs_options

   SUBROUTINE c__PDAF_etkf_ana_T(step, dim_p, dim_obs_p, dim_ens, state_p,  &
      ainv, ens_p, hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen,  &
      type_trans, debug, flag) bind(c)
      use PDAF_etkf_analysis_T
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! on exit: PE-local forecast state
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: state_p
      ! on exit: weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(inout) :: ens_p
      ! PE-local observed ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(inout) :: hz_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: hxbar_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Forgetting factor
      REAL(c_double), INTENT(in) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Type of ensemble transformation
      INTEGER(c_int), INTENT(in) :: type_trans
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A
      procedure(c__prodrinva_pdaf) :: u_prodrinva

      call PDAF_etkf_ana_T(step, dim_p, dim_obs_p, dim_ens, state_p, ainv,  &
         ens_p, hz_p, hxbar_p, obs_p, forget, u_prodrinva, screen, type_trans,  &
         debug, flag)

   END SUBROUTINE c__PDAF_etkf_ana_T

   SUBROUTINE c__PDAF_letkf_ana_fixed(domain_p, step, dim_l, dim_obs_l,  &
      dim_ens, state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, forget,  &
      u_prodrinva_l, screen, debug, flag) bind(c)
      use PDAF_letkf_analysis_fixed
      implicit none
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! State dimension on local analysis domain
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Size of obs. vector on local ana. domain
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Local forecast state
      REAL(c_double), DIMENSION(dim_l), INTENT(inout) :: state_l
      ! on exit: local weight matrix for ensemble transformation
      REAL(c_double), DIMENSION(dim_ens, dim_ens), INTENT(out) :: ainv_l
      ! Local state ensemble
      REAL(c_double), DIMENSION(dim_l, dim_ens), INTENT(inout) :: ens_l
      ! Local observed state ensemble (perturbation)
      REAL(c_double), DIMENSION(dim_obs_l, dim_ens), INTENT(inout) :: hz_l
      ! Local observed ensemble mean
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: hxbar_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Forgetting factor
      REAL(c_double), INTENT(inout) :: forget
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen
      ! Flag for writing debug output
      INTEGER(c_int), INTENT(in) :: debug
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag

      ! Provide product R^-1 A for local analysis domain
      procedure(c__prodrinva_l_pdaf) :: u_prodrinva_l

      call PDAF_letkf_ana_fixed(domain_p, step, dim_l, dim_obs_l, dim_ens,  &
         state_l, ainv_l, ens_l, hz_l, hxbar_l, obs_l, forget, u_prodrinva_l,  &
         screen, debug, flag)

   END SUBROUTINE c__PDAF_letkf_ana_fixed
end module pdaf_c_internal

