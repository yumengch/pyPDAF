module pdaf_c_cb_interface
implicit none

abstract interface
   SUBROUTINE c__add_obs_err_pdaf(step, dim_obs_p, C_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      IMPLICIT NONE
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Matrix to that observation covariance R is added
      REAL(c_double), DIMENSION(dim_obs_p,dim_obs_p), INTENT(inout) :: C_p
   END SUBROUTINE c__add_obs_err_pdaf

   SUBROUTINE c__init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, flag) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! filter type given in PDAF_init
      integer(c_int), intent(in) :: filtertype
      ! PE-local state dimension given by PDAF_init
      integer(c_int), intent(in) :: dim_p
      ! number of ensemble members
      integer(c_int), intent(in) :: dim_ens
      ! PE-local model state
      ! This array must be filled with the initial
      ! state of the model for SEEK, but it is not
      ! used for ensemble-based filters.
      !
      ! One can still make use of this array within
      ! this function.
      real(c_double), DIMENSION(dim_p), intent(inout) :: state_p
      ! This array is the inverse of matrix
      ! formed by right singular vectors of error
      ! covariance matrix of ensemble perturbations.
      !
      ! This array has to be filled in SEEK, but it is
      ! not used for ensemble-based filters.
      ! Nevertheless, one can still make use of this
      ! array within this function e.g.,
      ! for generating an initial ensemble perturbation
      ! from a given covariance matrix.
      !
      ! Dimension of this array is determined by the
      ! filter type.
      !
      ! - (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
      ! - (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
      ! - (1, 1) for (L)EnKF, particle filters and gen_obs
      real(c_double), DIMENSION(dim_ens - 1, dim_ens-1), intent(inout) :: uinv
      ! PE-local ensemble
      real(c_double), DIMENSION(dim_p, dim_ens), intent(inout) :: ens_p
      ! pdaf status flag
      integer(c_int), intent(inout) :: flag
   end subroutine c__init_ens_pdaf

   subroutine c__next_observation_pdaf(stepnow, nsteps, doexit, time) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! the current time step given by PDAF
      integer(c_int), intent(in)  :: stepnow
      ! number of forecast time steps until next assimilation;
      ! this can also be interpreted as
      ! number of assimilation function calls
      ! to perform a new assimilation
      integer(c_int), intent(out) :: nsteps
      ! whether to exit forecasting (1 for exit)
      integer(c_int), intent(out) :: doexit
      ! current model (physical) time
      real(c_double), intent(out) :: time
   end subroutine c__next_observation_pdaf

   subroutine c__collect_state_pdaf(dim_p, state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      ! local state vector
      real(c_double), DIMENSION(dim_p), intent(inout) :: state_p
   end subroutine c__collect_state_pdaf

   subroutine c__distribute_state_pdaf(dim_p, state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! PE-local state dimension
      integer(c_int), intent(in) :: dim_p
      ! PE-local state vector
      real(c_double), DIMENSION(dim_p), intent(inout) :: state_p
   end subroutine c__distribute_state_pdaf

   subroutine c__prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_l, &
      dim_obs_p, state_p, uinv, ens_p, flag) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! current time step
      ! (negative for call before analysis/preprocessing)
      integer(c_int), intent(in) :: step
      ! PE-local state vector dimension
      integer(c_int), intent(in) :: dim_p
      ! number of ensemble members
      integer(c_int), intent(in) :: dim_ens
      ! number of ensemble members run serially
      ! on each model task
      integer(c_int), intent(in) :: dim_ens_l
      ! PE-local dimension of observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! pe-local forecast/analysis state
      ! (the array 'state_p' is generally not
      ! initialised in the case of ESTKF/ETKF/EnKF/SEIK,
      ! so it can be used freely here.)
      real(c_double), DIMENSION(dim_p), intent(inout) :: state_p
      ! Inverse of the transformation matrix in ETKF and ESKTF;
      ! inverse of matrix formed by right singular vectors of error
      ! covariance matrix of ensemble perturbations in SEIK/SEEK.
      ! not used in EnKF.
      real(c_double), DIMENSION(dim_ens-1, dim_ens-1), intent(inout) :: uinv
      ! PE-local ensemble
      real(c_double), DIMENSION(dim_p, dim_ens), intent(inout) :: ens_p
      ! pdaf status flag
      integer(c_int), intent(in) :: flag
   end subroutine c__prepoststep_pdaf

   subroutine c__init_dim_obs_pdaf(step, dim_obs_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! current time step
      integer(c_int), intent(in)    :: step
      ! dimension of observation vector
      integer(c_int), intent(out) :: dim_obs_p
   end subroutine c__init_dim_obs_pdaf

   SUBROUTINE c__init_obs_pdaf(step, dim_obs_p, observation_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Size of the observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! Vector of observations
      real(c_double), DIMENSION(dim_obs_p), intent(out) :: observation_p
   END SUBROUTINE c__init_obs_pdaf

   SUBROUTINE c__init_obs_covar_pdaf(step, dim_obs, dim_obs_p, covar, obs_p, isdiag) bind(c)
      use iso_c_binding, only: c_double, c_int, c_bool
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Global size of observation vector
      integer(c_int), intent(in) :: dim_obs
      ! Size of process-local observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! Observation error covariance matrix
      real(c_double), intent(out), dimension(dim_obs,dim_obs) :: covar
      ! Process-local vector of observations
      real(c_double), intent(in), dimension(dim_obs_p) :: obs_p
      logical(c_bool), intent(out) :: isdiag
   END SUBROUTINE c__init_obs_covar_pdaf

   SUBROUTINE c__init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Size of observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! Vector of observations
      real(c_double), intent(in), dimension(dim_obs_p) :: obs_p
      ! Mean observation error variance
      real(c_double), intent(out) :: meanvar
   END SUBROUTINE c__init_obsvar_pdaf

   SUBROUTINE c__init_obsvars_pdaf(step, dim_obs_f, var_f) bind(c)
      use iso_c_binding, only: c_double, c_int
      IMPLICIT NONE
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Dimension of full observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! vector of observation error variances
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out) :: var_f
   END SUBROUTINE c__init_obsvars_pdaf

   SUBROUTINE c__prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of observations at current time step (i.e. the size of the observation vector)
      integer(c_int), intent(in) :: dim_obs_p
      ! Number of the columns in the matrix processes here.
      ! This is usually the ensemble size minus one
      ! (or the rank of the initial covariance matrix)
      integer(c_int), intent(in) :: rank
      ! Vector of observations
      real(c_double), intent(in), dimension(dim_obs_p) :: obs_p
      ! Input matrix provided by PDAF
      real(c_double), intent(in), dimension(dim_obs_p, rank) :: A_p
      ! Output matrix
      real(c_double), intent(out), dimension(dim_obs_p, rank) :: C_p
   END SUBROUTINE c__prodRinvA_pdaf

   SUBROUTINE c__obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Size of state vector
      ! (local part in case of parallel decomposed state)
      integer(c_int), intent(in) :: dim_p
      ! Size of PE-local observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! Model state vector
      real(c_double), intent(in), dimension(dim_p) :: state_p
      ! Observed state vector
      ! (i.e. the result after applying the observation operator to state_p)
      real(c_double), intent(inout), dimension(dim_obs_p) :: m_state_p
   END SUBROUTINE c__obs_op_pdaf

   SUBROUTINE c__g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, mstate_l) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Size of full observation vector for model sub-domain
      integer(c_int), intent(in) :: dim_obs_f
      ! Size of observation vector for local analysis domain
      integer(c_int), intent(in) :: dim_obs_l
      ! Full observation vector for model sub-domain
      integer(c_int), intent(in), dimension(dim_obs_f) :: mstate_f
      ! Observation vector for local analysis domain
      integer(c_int), intent(out), dimension(dim_obs_l) :: mstate_l
   END SUBROUTINE c__g2l_obs_pdaf

   subroutine c__g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none
      ! current time step
      integer(c_int), intent(in) :: step
      ! current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! pe-local full state dimension
      integer(c_int), intent(in) :: dim_p
      ! local state dimension
      integer(c_int), intent(in) :: dim_l
      ! pe-local full state vector
      real(c_double), dimension(dim_p), intent(in)    :: state_p
      ! state vector on local analysis domain
      real(c_double), dimension(dim_l), intent(out)   :: state_l
   end subroutine c__g2l_state_pdaf

   subroutine c__init_dim_l_pdaf(step, domain_p, dim_l) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      ! current time step
      integer(c_int), intent(in)  :: step
      ! current local analysis domain
      integer(c_int), intent(in)  :: domain_p
      ! local state dimension
      integer(c_int), intent(out) :: dim_l
   end subroutine c__init_dim_l_pdaf

   subroutine c__init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      ! index of current local analysis domain
      integer(c_int), intent(in)  :: domain_p
      ! current time step
      integer(c_int), intent(in)  :: step
      ! full dimension of observation vector
      integer(c_int), intent(in)  :: dim_obs_f
      ! local dimension of observation vector
      integer(c_int), intent(out) :: dim_obs_l
   end subroutine c__init_dim_obs_l_pdaf

   subroutine c__init_n_domains_p_pdaf(step, n_domains_p) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      ! current time step
      integer(c_int), intent(in)  :: step
      ! pe-local number of analysis domains
      integer(c_int), intent(out) :: n_domains_p
   end subroutine c__init_n_domains_p_pdaf


   SUBROUTINE c__init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Local size of the observation vector
      integer(c_int), intent(in) :: dim_obs_l
      ! Local vector of observations
      real(c_double), intent(out), dimension(dim_obs_l) :: observation_l
   END SUBROUTINE c__init_obs_l_pdaf

   SUBROUTINE c__init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, dim_obs_p, meanvar_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Local dimension of observation vector
      integer(c_int), intent(in) :: dim_obs_l
      ! Dimension of local observation vector
      integer(c_int), intent(in) :: dim_obs_p
      ! Local observation vector
      real(c_double), intent(in), dimension(dim_obs_p) :: obs_l
      ! Mean local observation error variance
      real(c_double), intent(out) :: meanvar_l
   END SUBROUTINE c__init_obsvar_l_pdaf

   SUBROUTINE c__init_obserr_f_pdaf(step, dim_obs_f, obs_f, obserr_f) bind(c)
      use iso_c_binding, only: c_double, c_int
      IMPLICIT NONE
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Full dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(in)    :: obs_f
      ! Full observation error stddev
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out)   :: obserr_f
   END SUBROUTINE c__init_obserr_f_pdaf

   subroutine c__l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none
      ! current time step
      integer(c_int), intent(in) :: step
      ! current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! local state dimension
      integer(c_int), intent(in) :: dim_l
      ! pe-local full state dimension
      integer(c_int), intent(in) :: dim_p
      ! state vector on local analysis domain
      real(c_double), DIMENSION(dim_l), intent(in)    :: state_l
      ! pe-local full state vector
      real(c_double), DIMENSION(dim_p), intent(inout) :: state_p
   end subroutine c__l2g_state_pdaf

   SUBROUTINE c__prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer(c_int), intent(in) :: dim_obs_l
      ! Number of the columns in the matrix processes here.
      ! This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
      integer(c_int), intent(in) :: rank
      ! Local vector of observations
      real(c_double), intent(in), dimension(dim_obs_l) :: obs_l
      ! Input matrix provided by PDAF
      real(c_double), intent(inout), dimension(dim_obs_l, rank) :: A_l
      ! Output matrix
      real(c_double), intent(out), dimension(dim_obs_l, rank) :: C_l
   END SUBROUTINE c__prodRinvA_l_pdaf

   subroutine c__localize_covar_pdaf(dim_p, dim_obs, hp_p, hph) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none
      ! pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      ! number of observations
      integer(c_int), intent(in) :: dim_obs
      ! pe local part of matrix hp
      real(c_double), DIMENSION(dim_obs, dim_p), intent(inout) :: hp_p
      ! matrix hph
      real(c_double), DIMENSION(dim_obs, dim_obs), intent(inout) :: hph
   end subroutine c__localize_covar_pdaf

   SUBROUTINE c__localize_covar_serial_pdaf(iobs, dim_p, dim_obs, HP_p, HXY_p) bind(c)
      use iso_c_binding, only: c_int, c_double
      IMPLICIT NONE
      ! Index of current observation
      INTEGER(c_int), INTENT(in) :: iobs
      ! Process-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Number of observations
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Process-local part of matrix HP for observation iobs
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: HP_p
      ! Process-local part of matrix HX(HX_all) for full observations
      REAL(c_double), DIMENSION(dim_obs), INTENT(inout) :: HXY_p
   END SUBROUTINE c__localize_covar_serial_pdaf

   SUBROUTINE c__likelihood_pdaf(step, dim_obs_p, obs_p, resid, likely) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of observations at current time step (i.e. the size of the observation vector)
      integer(c_int), intent(in) :: dim_obs_p
      ! Vector of observations
      real(c_double), intent(in), dimension(dim_obs_p) :: obs_p
      ! Input vector holding the residual
      real(c_double), intent(in), dimension(dim_obs_p) :: resid
      ! Output value of the likelihood
      real(c_double), intent(out) :: likely
   END SUBROUTINE c__likelihood_pdaf

   SUBROUTINE c__likelihood_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer(c_int), intent(in) :: dim_obs_l
      ! Local vector of observations
      real(c_double), intent(in), dimension(dim_obs_l) :: obs_l
      ! nput vector holding the local residual
      real(c_double), intent(inout), dimension(dim_obs_l) :: resid_l
      ! Output value of the local likelihood
      real(c_double), intent(out) :: likely_l
   END SUBROUTINE c__likelihood_l_pdaf

   SUBROUTINE c__get_obs_f_pdaf(step, dim_obs_f, observation_f) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Size of the full observation vector
      integer(c_int), intent(in) :: dim_obs_f
      ! Full vector of synthetic observations (process-local)
      real(c_double), intent(in), dimension(dim_obs_f) :: observation_f
   END SUBROUTINE c__get_obs_f_pdaf

   SUBROUTINE c__cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vcv_p, cv_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Iteration of optimization
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! PE-local dimension of control vector
      INTEGER(c_int), INTENT(in) :: dim_cv_ens_p
      ! PE-local ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! PE-local input vector
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: Vcv_p
      ! PE-local result vector
      REAL(c_double), DIMENSION(dim_cv_ens_p), INTENT(inout) :: cv_p
   END SUBROUTINE c__cvt_adj_ens_pdaf

   SUBROUTINE c__cvt_adj_pdaf(iter, dim_p, dim_cvec, Vcv_p, cv_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Iteration of optimization
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Dimension of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! PE-local result vector (state vector increment)
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: Vcv_p
      ! PE-local control vector
      REAL(c_double), DIMENSION(dim_cvec), INTENT(inout) :: cv_p
   END SUBROUTINE c__cvt_adj_pdaf

   SUBROUTINE c__cvt_pdaf(iter, dim_p, dim_cvec, cv_p, Vv_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Iteration of optimization
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Dimension of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec
      ! PE-local control vector
      REAL(c_double), DIMENSION(dim_cvec), INTENT(in) :: cv_p
      ! PE-local result vector (state vector increment)
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: Vv_p
   END SUBROUTINE c__cvt_pdaf

   SUBROUTINE c__cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, Vv_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      IMPLICIT NONE
      ! Iteration of optimization
      INTEGER(c_int), INTENT(in) :: iter
      ! PE-local dimension of state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Dimension of control vector
      INTEGER(c_int), INTENT(in) :: dim_cvec_ens
      ! PE-local ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! PE-local control vector
      REAL(c_double), DIMENSION(dim_cvec_ens), INTENT(in) :: v_p
      ! PE-local state increment
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: Vv_p
   END SUBROUTINE c__cvt_ens_pdaf

   SUBROUTINE c__obs_op_adj_pdaf(step, dim_p, dim_obs_p, m_state_p, state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Dimension of observed state
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local observed state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: m_state_p
      ! PE-local model state
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
   END SUBROUTINE c__obs_op_adj_pdaf

   SUBROUTINE c__likelihood_hyb_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, gamma, likely_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer(c_int), intent(in) :: dim_obs_l
      ! Local vector of observations
      real(c_double), dimension(dim_obs_l), intent(in) :: obs_l
      ! Hybrid weight provided by PDAF
      real(c_double), intent(in) :: gamma
      ! Input vector holding the local residual
      real(c_double), dimension(dim_obs_l), intent(inout) :: resid_l
      ! Output value of the local likelihood
      real(c_double), intent(out) :: likely_l
   END SUBROUTINE c__likelihood_hyb_l_pdaf

   SUBROUTINE c__prodRinvA_hyb_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, gamma, A_l, C_l) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! Index of current local analysis domain
      integer(c_int), intent(in) :: domain_p
      ! Current time step
      integer(c_int), intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer(c_int), intent(in) :: dim_obs_l
      ! Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
      integer(c_int), intent(in) :: rank
      ! Local vector of observations
      real(c_double), dimension(dim_obs_l), intent(in) :: obs_l
      ! Hybrid weight provided by PDAF
      real(c_double), intent(in) :: gamma
      ! Input matrix provided by PDAF
      real(c_double), dimension(dim_obs_l, rank), intent(inout) :: A_l
      ! Output matrix
      real(c_double), dimension(dim_obs_l, rank), intent(out) :: C_l
   END SUBROUTINE c__prodRinvA_hyb_l_pdaf
end interface

end module pdaf_c_cb_interface