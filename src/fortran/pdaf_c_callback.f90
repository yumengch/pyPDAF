MODULE pdaf_c_callback
use PDAF

implicit none

contains
   SUBROUTINE c__PDAFomi_init_obs_f_cb(step, dim_obs_f, observation_f) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of full observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out) :: observation_f


      call PDAFomi_init_obs_f_cb(step, dim_obs_f, observation_f)

   END SUBROUTINE c__PDAFomi_init_obs_f_cb

   SUBROUTINE c__PDAFomi_init_obsvar_cb(step, dim_obs_p, obs_p, meanvar) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Mean observation error variance
      REAL(c_double), INTENT(out) :: meanvar


      call PDAFomi_init_obsvar_cb(step, dim_obs_p, obs_p, meanvar)

   END SUBROUTINE c__PDAFomi_init_obsvar_cb

   SUBROUTINE c__PDAFomi_init_obsvars_f_cb(step, dim_obs_f, var_f) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of full observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! vector of observation error variances
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out) :: var_f


      call PDAFomi_init_obsvars_f_cb(step, dim_obs_f, var_f)

   END SUBROUTINE c__PDAFomi_init_obsvars_f_cb

   SUBROUTINE c__PDAFomi_g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l,  &
      ostate_f, ostate_l) bind(c)
      use iso_c_binding

      ! Index of current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of full PE-local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Dimension of local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Full PE-local obs.ervation vector
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(in) :: ostate_f
      ! Observation vector on local domain
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(out) :: ostate_l


      call PDAFomi_g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l, ostate_f,  &
         ostate_l)

   END SUBROUTINE c__PDAFomi_g2l_obs_cb

   SUBROUTINE c__PDAFomi_init_obs_l_cb(domain_p, step, dim_obs_l,  &
      observation_l) bind(c)
      use iso_c_binding

      ! Index of current local analysis domain index
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(out) :: observation_l


      call PDAFomi_init_obs_l_cb(domain_p, step, dim_obs_l, observation_l)

   END SUBROUTINE c__PDAFomi_init_obs_l_cb

   SUBROUTINE c__PDAFomi_init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l,  &
      meanvar_l) bind(c)
      use iso_c_binding

      ! Index of current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Local dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Local observation vector
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Mean local observation error variance
      REAL(c_double), INTENT(out) :: meanvar_l


      call PDAFomi_init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l, meanvar_l)

   END SUBROUTINE c__PDAFomi_init_obsvar_l_cb

   SUBROUTINE c__PDAFomi_prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l,  &
      a_l, c_l) bind(c)
      use iso_c_binding

      ! Index of current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Local vector of observations
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Input matrix
      REAL(c_double), DIMENSION(dim_obs_l, rank), INTENT(inout) :: a_l
      ! Output matrix
      REAL(c_double), DIMENSION(dim_obs_l, rank), INTENT(out) :: c_l


      call PDAFomi_prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l, a_l, c_l)

   END SUBROUTINE c__PDAFomi_prodRinvA_l_cb

   SUBROUTINE c__PDAFomi_likelihood_l_cb(domain_p, step, dim_obs_l, obs_l,  &
      resid_l, lhood_l) bind(c)
      use iso_c_binding

      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! PE-local vector of observations
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Input vector of residuum
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(inout) :: resid_l
      ! Output vector - log likelihood
      REAL(c_double), INTENT(out) :: lhood_l


      call PDAFomi_likelihood_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l,  &
         lhood_l)

   END SUBROUTINE c__PDAFomi_likelihood_l_cb

   SUBROUTINE c__PDAFomi_prodRinvA_hyb_l_cb(domain_p, step, dim_obs_l, rank,  &
      obs_l, alpha, a_l, c_l) bind(c)
      use iso_c_binding

      ! Index of current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: rank
      ! Local vector of observations
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: alpha
      ! Input matrix
      REAL(c_double), DIMENSION(dim_obs_l, rank), INTENT(inout) :: a_l
      ! Output matrix
      REAL(c_double), DIMENSION(dim_obs_l, rank), INTENT(out) :: c_l


      call PDAFomi_prodRinvA_hyb_l_cb(domain_p, step, dim_obs_l, rank, obs_l,  &
         alpha, a_l, c_l)

   END SUBROUTINE c__PDAFomi_prodRinvA_hyb_l_cb

   SUBROUTINE c__PDAFomi_likelihood_hyb_l_cb(domain_p, step, dim_obs_l, obs_l,  &
      resid_l, alpha, lhood_l) bind(c)
      use iso_c_binding

      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! PE-local vector of observations
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(in) :: obs_l
      ! Input vector of residuum
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(inout) :: resid_l
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: alpha
      ! Output vector - log likelihood
      REAL(c_double), INTENT(out) :: lhood_l


      call PDAFomi_likelihood_hyb_l_cb(domain_p, step, dim_obs_l, obs_l,  &
         resid_l, alpha, lhood_l)

   END SUBROUTINE c__PDAFomi_likelihood_hyb_l_cb

   SUBROUTINE c__PDAFomi_prodRinvA_cb(step, dim_obs_p, ncol, obs_p, a_p,  &
      c_p) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of PE-local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Number of columns in A_p and C_p
      INTEGER(c_int), INTENT(in) :: ncol
      ! PE-local vector of observations
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Input matrix
      REAL(c_double), DIMENSION(dim_obs_p, ncol), INTENT(in) :: a_p
      ! Output matrix
      REAL(c_double), DIMENSION(dim_obs_p, ncol), INTENT(out) :: c_p


      call PDAFomi_prodRinvA_cb(step, dim_obs_p, ncol, obs_p, a_p, c_p)

   END SUBROUTINE c__PDAFomi_prodRinvA_cb

   SUBROUTINE c__PDAFomi_likelihood_cb(step, dim_obs, obs, resid, lhood) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! PE-local vector of observations
      REAL(c_double), DIMENSION(dim_obs), INTENT(in) :: obs
      ! Input vector of residuum
      REAL(c_double), DIMENSION(dim_obs), INTENT(in) :: resid
      ! Output vector - log likelihood
      REAL(c_double), INTENT(out) :: lhood


      call PDAFomi_likelihood_cb(step, dim_obs, obs, resid, lhood)

   END SUBROUTINE c__PDAFomi_likelihood_cb

   SUBROUTINE c__PDAFomi_add_obs_error_cb(step, dim_obs_p, c_p) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of PE-local observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Matrix to which R is added
      REAL(c_double), DIMENSION(dim_obs_p,dim_obs_p), INTENT(inout) :: c_p


      call PDAFomi_add_obs_error_cb(step, dim_obs_p, c_p)

   END SUBROUTINE c__PDAFomi_add_obs_error_cb

   SUBROUTINE c__PDAFomi_init_obscovar_cb(step, dim_obs, dim_obs_p, covar,  &
      m_state_p, isdiag) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! PE-local dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Observation error covar. matrix
      REAL(c_double), DIMENSION(dim_obs,dim_obs), INTENT(out) :: covar
      ! Observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: m_state_p
      ! Whether matrix R is diagonal
      LOGICAL(c_bool), INTENT(out) :: isdiag


      call PDAFomi_init_obscovar_cb(step, dim_obs, dim_obs_p, covar, m_state_p,  &
         isdiag)

   END SUBROUTINE c__PDAFomi_init_obscovar_cb

   SUBROUTINE c__PDAFomi_init_obserr_f_cb(step, dim_obs_f, obs_f, obserr_f) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Full dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(in) :: obs_f
      ! Full observation error stddev
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out) :: obserr_f


      call PDAFomi_init_obserr_f_cb(step, dim_obs_f, obs_f, obserr_f)

   END SUBROUTINE c__PDAFomi_init_obserr_f_cb

   SUBROUTINE c__PDAFomi_localize_covar_cb(dim_p, dim_obs, hp_p, hph) bind(c)
      use iso_c_binding

      ! Process-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Number of observations
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Process-local part of matrix HP
      REAL(c_double), DIMENSION(dim_obs, dim_p), INTENT(inout) :: hp_p
      ! Matrix HPH
      REAL(c_double), DIMENSION(dim_obs, dim_obs), INTENT(inout) :: hph


      call PDAFomi_localize_covar_cb(dim_p, dim_obs, hp_p, hph)

   END SUBROUTINE c__PDAFomi_localize_covar_cb

   SUBROUTINE c__PDAFomi_localize_covar_serial_cb(iobs, dim_p, dim_obs, hp_p,  &
      hxy_p) bind(c)
      use iso_c_binding

      ! Index of current observation
      INTEGER(c_int), INTENT(in) :: iobs
      ! Process-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Number of observations
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Process-local part of matrix HP for observation iobs
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: hp_p
      ! Process-local part of matrix HX(HX_all) for full observations
      REAL(c_double), DIMENSION(dim_obs), INTENT(inout) :: hxy_p


      call PDAFomi_localize_covar_serial_cb(iobs, dim_p, dim_obs, hp_p, hxy_p)

   END SUBROUTINE c__PDAFomi_localize_covar_serial_cb

   SUBROUTINE c__PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, resid_l,  &
      obs_l) bind(c)
      use iso_c_binding

      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! PE-local dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs_l
      ! Input vector of residuum
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(inout) :: resid_l
      ! Input vector of local observations
      REAL(c_double), DIMENSION(dim_obs_l), INTENT(inout) :: obs_l


      call PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, resid_l, obs_l)

   END SUBROUTINE c__PDAFomi_omit_by_inno_l_cb

   SUBROUTINE c__PDAFomi_omit_by_inno_cb(dim_obs_f, resid_f, obs_f) bind(c)
      use iso_c_binding

      ! Full dimension of obs. vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Input vector of residuum
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(inout) :: resid_f
      ! Input vector of full observations
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(inout) :: obs_f


      call PDAFomi_omit_by_inno_cb(dim_obs_f, resid_f, obs_f)

   END SUBROUTINE c__PDAFomi_omit_by_inno_cb

   SUBROUTINE c__PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l,  &
      state_l) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! PE-local full state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local full state vector
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: state_p
      ! Local state dimension
      INTEGER(c_int), INTENT(in) :: dim_l
      ! State vector on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(out) :: state_l


      call PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l)

   END SUBROUTINE c__PDAFlocal_g2l_cb

   SUBROUTINE c__PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p,  &
      state_p) bind(c)
      use iso_c_binding

      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Current local analysis domain
      INTEGER(c_int), INTENT(in) :: domain_p
      ! Local state dimension
      INTEGER(c_int), INTENT(in) :: dim_l
      ! State vector on local analysis domain
      REAL(c_double), DIMENSION(dim_l), INTENT(in) :: state_l
      ! PE-local full state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! PE-local full state vector
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p


      call PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p)

   END SUBROUTINE c__PDAFlocal_l2g_cb
END MODULE pdaf_c_callback
