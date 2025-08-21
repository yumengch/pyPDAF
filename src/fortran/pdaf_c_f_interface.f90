module pdaf_c_f_interface
use pdaf_c_cb_interface
implicit none

procedure(c__init_ens_pdaf), pointer :: init_ens_pdaf_c_ptr => null()
procedure(c__add_obs_err_pdaf), pointer :: add_obs_err_pdaf_c_ptr => null()
procedure(c__next_observation_pdaf), pointer :: next_observation_pdaf_c_ptr => null()
procedure(c__collect_state_pdaf), pointer :: collect_state_pdaf_c_ptr => null()
procedure(c__distribute_state_pdaf), pointer :: distribute_state_pdaf_c_ptr => null()
procedure(c__prepoststep_pdaf), pointer :: prepoststep_pdaf_c_ptr => null()
procedure(c__init_dim_obs_pdaf), pointer :: init_dim_obs_pdaf_c_ptr => null()
procedure(c__init_dim_obs_f_pdaf), pointer :: init_dim_obs_f_pdaf_c_ptr => null()
procedure(c__init_obs_pdaf), pointer :: init_obs_pdaf_c_ptr => null()
procedure(c__init_obs_covar_pdaf), pointer :: init_obs_covar_pdaf_c_ptr => null()
procedure(c__init_obsvar_pdaf), pointer :: init_obsvar_pdaf_c_ptr => null()
procedure(c__init_obsvars_pdaf), pointer :: init_obsvars_pdaf_c_ptr => null()
procedure(c__prodRinvA_pdaf), pointer :: prodRinvA_pdaf_c_ptr => null()
procedure(c__obs_op_pdaf), pointer :: obs_op_pdaf_c_ptr => null()
procedure(c__obs_op_f_pdaf), pointer :: obs_op_f_pdaf_c_ptr => null()
procedure(c__g2l_obs_pdaf), pointer :: g2l_obs_pdaf_c_ptr => null()
procedure(c__g2l_state_pdaf), pointer :: g2l_state_pdaf_c_ptr => null()
procedure(c__init_dim_l_pdaf), pointer :: init_dim_l_pdaf_c_ptr => null()
procedure(c__init_dim_obs_l_pdaf), pointer :: init_dim_obs_l_pdaf_c_ptr => null()
procedure(c__init_n_domains_p_pdaf), pointer :: init_n_domains_p_pdaf_c_ptr => null()
procedure(c__init_obs_f_pdaf), pointer :: init_obs_f_pdaf_c_ptr => null()
procedure(c__init_obs_l_pdaf), pointer :: init_obs_l_pdaf_c_ptr => null()
procedure(c__init_obsvar_l_pdaf), pointer :: init_obsvar_l_pdaf_c_ptr => null()
procedure(c__init_obserr_f_pdaf), pointer :: init_obserr_f_pdaf_c_ptr => null()
procedure(c__l2g_state_pdaf), pointer :: l2g_state_pdaf_c_ptr => null()
procedure(c__prodRinvA_l_pdaf), pointer :: prodRinvA_l_pdaf_c_ptr => null()
procedure(c__localize_covar_pdaf), pointer :: localize_covar_pdaf_c_ptr => null()
procedure(c__localize_covar_serial_pdaf), pointer :: localize_covar_serial_pdaf_c_ptr => null()
procedure(c__likelihood_pdaf), pointer :: likelihood_pdaf_c_ptr => null()
procedure(c__likelihood_l_pdaf), pointer :: likelihood_l_pdaf_c_ptr => null()
procedure(c__get_obs_f_pdaf), pointer :: get_obs_f_pdaf_c_ptr => null()
procedure(c__cvt_adj_ens_pdaf), pointer :: cvt_adj_ens_pdaf_c_ptr => null()
procedure(c__cvt_adj_pdaf), pointer :: cvt_adj_pdaf_c_ptr => null()
procedure(c__cvt_pdaf), pointer :: cvt_pdaf_c_ptr => null()
procedure(c__cvt_ens_pdaf), pointer :: cvt_ens_pdaf_c_ptr => null()
procedure(c__obs_op_adj_pdaf), pointer :: obs_op_adj_pdaf_c_ptr => null()
procedure(c__obs_op_lin_pdaf), pointer :: obs_op_lin_pdaf_c_ptr => null()
procedure(c__dist_stateinc_pdaf), pointer :: dist_stateinc_pdaf_c_ptr => null()
procedure(c__likelihood_hyb_l_pdaf), pointer :: likelihood_hyb_l_pdaf_c_ptr => null()
procedure(c__prodRinvA_hyb_l_pdaf), pointer :: prodRinvA_hyb_l_pdaf_c_ptr => null()



contains
   SUBROUTINE f__add_obs_err_pdaf(step, dim_obs_p, C_p)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: step
      INTEGER, INTENT(in) :: dim_obs_p
      REAL, DIMENSION(dim_obs_p,dim_obs_p), INTENT(inout) :: C_p
      call add_obs_err_pdaf_c_ptr(step, dim_obs_p, C_p)
   END SUBROUTINE f__add_obs_err_pdaf

   SUBROUTINE f__init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, flag)
      implicit none
      integer, intent(in) :: filtertype
      integer, intent(in) :: dim_p
      integer, intent(in) :: dim_ens
      real, DIMENSION(dim_p), intent(inout) :: state_p
      real, DIMENSION(dim_ens - 1, dim_ens-1), intent(inout) :: uinv
      real, DIMENSION(dim_p, dim_ens), intent(inout) :: ens_p
      integer, intent(inout) :: flag

      call init_ens_pdaf_c_ptr(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, flag)
   end subroutine f__init_ens_pdaf

   subroutine f__next_observation_pdaf(stepnow, nsteps, doexit, time)
      implicit none
      ! the current time step given by PDAF
      integer, intent(in)  :: stepnow
      ! number of forecast time steps until next assimilation;
      ! this can also be interpreted as
      ! number of assimilation function calls
      ! to perform a new assimilation
      integer, intent(out) :: nsteps
      ! whether to exit forecasting (1 for exit)
      integer, intent(out) :: doexit
      ! current model (physical) time
      real, intent(out) :: time

      call next_observation_pdaf_c_ptr(stepnow, nsteps, doexit, time)
   end subroutine f__next_observation_pdaf

   subroutine f__collect_state_pdaf(dim_p, state_p)
      implicit none
      ! pe-local state dimension
      integer, intent(in) :: dim_p
      ! local state vector
      real, DIMENSION(dim_p), intent(inout) :: state_p
      call collect_state_pdaf_c_ptr(dim_p, state_p)
   end subroutine f__collect_state_pdaf

   subroutine f__distribute_state_pdaf(dim_p, state_p)
      implicit none
      ! PE-local state dimension
      integer, intent(in) :: dim_p
      ! PE-local state vector
      real, DIMENSION(dim_p), intent(inout) :: state_p
      call distribute_state_pdaf_c_ptr(dim_p, state_p)
   end subroutine f__distribute_state_pdaf

   subroutine f__prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_l, &
      dim_obs_p, state_p, uinv, ens_p, flag)
      implicit none
      ! current time step
      ! (negative for call before analysis/preprocessing)
      integer, intent(in) :: step
      ! PE-local state vector dimension
      integer, intent(in) :: dim_p
      ! number of ensemble members
      integer, intent(in) :: dim_ens
      ! number of ensemble members run serially
      ! on each model task
      integer, intent(in) :: dim_ens_l
      ! PE-local dimension of observation vector
      integer, intent(in) :: dim_obs_p
      ! pe-local forecast/analysis state
      ! (the array 'state_p' is generally not
      ! initialised in the case of ESTKF/ETKF/EnKF/SEIK,
      ! so it can be used freely here.)
      real, DIMENSION(dim_p), intent(inout) :: state_p
      ! Inverse of the transformation matrix in ETKF and ESKTF;
      ! inverse of matrix formed by right singular vectors of error
      ! covariance matrix of ensemble perturbations in SEIK/SEEK.
      ! not used in EnKF.
      real, DIMENSION(dim_ens-1, dim_ens-1), intent(inout) :: uinv
      ! PE-local ensemble
      real, DIMENSION(dim_p, dim_ens), intent(inout) :: ens_p
      ! pdaf status flag
      integer, intent(in) :: flag
      call prepoststep_pdaf_c_ptr(step, dim_p, dim_ens, dim_ens_l, &
         dim_obs_p, state_p, uinv, ens_p, flag)
   end subroutine f__prepoststep_pdaf

   subroutine f__init_dim_obs_pdaf(step, dim_obs_p)
      implicit none
      ! current time step
      integer, intent(in)    :: step
      ! dimension of observation vector
      integer, intent(out) :: dim_obs_p

      call init_dim_obs_pdaf_c_ptr(step, dim_obs_p)
   end subroutine f__init_dim_obs_pdaf

   subroutine f__init_dim_obs_f_pdaf(step, dim_obs_p)
      implicit none
      ! current time step
      integer, intent(in)    :: step
      ! dimension of observation vector
      integer, intent(out) :: dim_obs_p
      call init_dim_obs_f_pdaf_c_ptr(step, dim_obs_p)
   end subroutine f__init_dim_obs_f_pdaf

   SUBROUTINE f__init_obs_pdaf(step, dim_obs_p, observation_p)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of the observation vector
      integer, intent(in) :: dim_obs_p
      ! Vector of observations
      real, DIMENSION(dim_obs_p), intent(out) :: observation_p
      call init_obs_pdaf_c_ptr(step, dim_obs_p, observation_p)
   END SUBROUTINE f__init_obs_pdaf

   SUBROUTINE f__init_obs_covar_pdaf(step, dim_obs, dim_obs_p, covar, obs_p, isdiag)
      use iso_c_binding, only: c_bool
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Global size of observation vector
      integer, intent(in) :: dim_obs
      ! Size of process-local observation vector
      integer, intent(in) :: dim_obs_p
      ! Observation error covariance matrix
      real, intent(out), dimension(dim_obs,dim_obs) :: covar
      ! Process-local vector of observations
      real, intent(in), dimension(dim_obs_p) :: obs_p
      logical, intent(out) :: isdiag

      logical(c_bool) :: isdiag_c
      call init_obs_covar_pdaf_c_ptr(step, dim_obs, dim_obs_p, covar, obs_p, isdiag_c)
      isdiag = isdiag_c
   END SUBROUTINE f__init_obs_covar_pdaf

   SUBROUTINE f__init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of observation vector
      integer, intent(in) :: dim_obs_p
      ! Vector of observations
      real, intent(in), dimension(dim_obs_p) :: obs_p
      ! Mean observation error variance
      real, intent(out) :: meanvar
      call init_obsvar_pdaf_c_ptr(step, dim_obs_p, obs_p, meanvar)
   END SUBROUTINE f__init_obsvar_pdaf

   SUBROUTINE f__init_obsvars_pdaf(step, dim_obs_f, var_f)
      IMPLICIT NONE
      ! Current time step
      INTEGER, INTENT(in) :: step
      ! Dimension of full observation vector
      INTEGER, INTENT(in) :: dim_obs_f
      ! vector of observation error variances
      REAL, DIMENSION(dim_obs_f), INTENT(out) :: var_f

      call init_obsvars_pdaf_c_ptr(step, dim_obs_f, var_f)
   END SUBROUTINE f__init_obsvars_pdaf

   SUBROUTINE f__prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Number of observations at current time step (i.e. the size of the observation vector)
      integer, intent(in) :: dim_obs_p
      ! Number of the columns in the matrix processes here.
      ! This is usually the ensemble size minus one
      ! (or the rank of the initial covariance matrix)
      integer, intent(in) :: rank
      ! Vector of observations
      real, intent(in), dimension(dim_obs_p) :: obs_p
      ! Input matrix provided by PDAF
      real, intent(in), dimension(dim_obs_p, rank) :: A_p
      ! Output matrix
      real, intent(out), dimension(dim_obs_p, rank) :: C_p
      call prodRinvA_pdaf_c_ptr(step, dim_obs_p, rank, obs_p, A_p, C_p)
   END SUBROUTINE f__prodRinvA_pdaf

   SUBROUTINE f__obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of state vector
      ! (local part in case of parallel decomposed state)
      integer, intent(in) :: dim_p
      ! Size of PE-local observation vector
      integer, intent(in) :: dim_obs_p
      ! Model state vector
      real, intent(in), dimension(dim_p) :: state_p
      ! Observed state vector
      ! (i.e. the result after applying the observation operator to state_p)
      real, intent(inout), dimension(dim_obs_p) :: m_state_p
      call obs_op_pdaf_c_ptr(step, dim_p, dim_obs_p, state_p, m_state_p)
   END SUBROUTINE f__obs_op_pdaf

   SUBROUTINE f__obs_op_f_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of state vector (local part in case of parallel decomposed state)
      integer, intent(in) :: dim_p
      ! Size of observation vector
      integer, intent(in) :: dim_obs_p
      ! Model state vector
      real, intent(in), dimension(dim_p) :: state_p
      ! Observed state vector (i.e. the result after applying the observation operator to state_p)
      real, intent(out), dimension(dim_obs_p) :: m_state_p
      call obs_op_f_pdaf_c_ptr(step, dim_p, dim_obs_p, state_p, m_state_p)
   END SUBROUTINE f__obs_op_f_pdaf

   SUBROUTINE f__g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, dim_p, mstate_l, dim_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Size of full observation vector for model sub-domain
      integer, intent(in) :: dim_obs_f
      ! Size of observation vector for local analysis domain
      integer, intent(in) :: dim_obs_l
      ! Size of full observation vector for model sub-domain
      integer, intent(in) :: dim_p
      ! Size of observation vector for local analysis domain
      integer, intent(in) :: dim_l
      ! Full observation vector for model sub-domain
      integer, intent(in), dimension(dim_p) :: mstate_f
      ! Observation vector for local analysis domain
      integer, intent(out), dimension(dim_l) :: mstate_l
      call g2l_obs_pdaf_c_ptr(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, dim_p, mstate_l, dim_l)
   END SUBROUTINE f__g2l_obs_pdaf

   subroutine f__g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

      implicit none
      ! current time step
      integer, intent(in) :: step
      ! current local analysis domain
      integer, intent(in) :: domain_p
      ! pe-local full state dimension
      integer, intent(in) :: dim_p
      ! local state dimension
      integer, intent(in) :: dim_l
      ! pe-local full state vector
      real, dimension(dim_p), intent(in)    :: state_p
      ! state vector on local analysis domain
      real, dimension(dim_l), intent(out)   :: state_l
      call g2l_state_pdaf_c_ptr(step, domain_p, dim_p, state_p, dim_l, state_l)
   end subroutine f__g2l_state_pdaf

   subroutine f__init_dim_l_pdaf(step, domain_p, dim_l)
      implicit none
      ! current time step
      integer, intent(in)  :: step
      ! current local analysis domain
      integer, intent(in)  :: domain_p
      ! local state dimension
      integer, intent(out) :: dim_l
      call init_dim_l_pdaf_c_ptr(step, domain_p, dim_l)
   end subroutine f__init_dim_l_pdaf

   subroutine f__init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)
      implicit none
      ! index of current local analysis domain
      integer, intent(in)  :: domain_p
      ! current time step
      integer, intent(in)  :: step
      ! full dimension of observation vector
      integer, intent(in)  :: dim_obs_f
      ! local dimension of observation vector
      integer, intent(out) :: dim_obs_l
      call init_dim_obs_l_pdaf_c_ptr(domain_p, step, dim_obs_f, dim_obs_l)
   end subroutine f__init_dim_obs_l_pdaf

   subroutine f__init_n_domains_p_pdaf(step, n_domains_p)
      implicit none
      ! current time step
      integer, intent(in)  :: step
      ! pe-local number of analysis domains
      integer, intent(out) :: n_domains_p
      call init_n_domains_p_pdaf_c_ptr(step, n_domains_p)
   end subroutine f__init_n_domains_p_pdaf

   SUBROUTINE f__init_obs_f_pdaf(step, dim_obs_f, observation_f)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of the full observation vector
      integer, intent(in) :: dim_obs_f
      ! Full vector of observations
      real, intent(out), dimension(dim_obs_f) :: observation_f
      call init_obs_f_pdaf_c_ptr(step, dim_obs_f, observation_f)
   END SUBROUTINE f__init_obs_f_pdaf

   SUBROUTINE f__init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Local size of the observation vector
      integer, intent(in) :: dim_obs_l
      ! Local vector of observations
      real, intent(out), dimension(dim_obs_l) :: observation_l
      call init_obs_l_pdaf_c_ptr(domain_p, step, dim_obs_l, observation_l)
   END SUBROUTINE f__init_obs_l_pdaf

   SUBROUTINE f__init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, dim_obs_p, meanvar_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Local dimension of observation vector
      integer, intent(in) :: dim_obs_l
      ! Dimension of local observation vector
      integer, intent(in) :: dim_obs_p
      ! Local observation vector
      real, intent(in), dimension(dim_obs_p) :: obs_l
      ! Mean local observation error variance
      real, intent(out) :: meanvar_l
      call init_obsvar_l_pdaf_c_ptr(domain_p, step, dim_obs_l, obs_l, dim_obs_p, meanvar_l)
   END SUBROUTINE f__init_obsvar_l_pdaf

   SUBROUTINE f__init_obserr_f_pdaf(step, dim_obs_f, obs_f, obserr_f)
      IMPLICIT NONE
      ! Current time step
      INTEGER, INTENT(in) :: step
      ! Full dimension of observation vector
      INTEGER, INTENT(in) :: dim_obs_f
      ! Full observation vector
      REAL, DIMENSION(dim_obs_f), INTENT(in)    :: obs_f
      ! Full observation error stddev
      REAL, DIMENSION(dim_obs_f), INTENT(out)   :: obserr_f
      call init_obserr_f_pdaf_c_ptr(step, dim_obs_f, obs_f, obserr_f)
   END SUBROUTINE f__init_obserr_f_pdaf

   subroutine f__l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)
      implicit none
      ! current time step
      integer, intent(in) :: step
      ! current local analysis domain
      integer, intent(in) :: domain_p
      ! local state dimension
      integer, intent(in) :: dim_l
      ! pe-local full state dimension
      integer, intent(in) :: dim_p
      ! state vector on local analysis domain
      real, DIMENSION(dim_l), intent(in)    :: state_l
      ! pe-local full state vector
      real, DIMENSION(dim_p), intent(inout) :: state_p
      call l2g_state_pdaf_c_ptr(step, domain_p, dim_l, state_l, dim_p, state_p)
   end subroutine f__l2g_state_pdaf

   SUBROUTINE f__prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer, intent(in) :: dim_obs_l
      ! Number of the columns in the matrix processes here.
      ! This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
      integer, intent(in) :: rank
      ! Local vector of observations
      real, intent(in), dimension(dim_obs_l) :: obs_l
      ! Input matrix provided by PDAF
      real, intent(inout), dimension(dim_obs_l, rank) :: A_l
      ! Output matrix
      real, intent(out), dimension(dim_obs_l, rank) :: C_l
      call prodRinvA_l_pdaf_c_ptr(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
   END SUBROUTINE f__prodRinvA_l_pdaf

   subroutine f__localize_covar_pdaf(dim_p, dim_obs, hp_p, hph)
      implicit none
      ! pe-local state dimension
      integer, intent(in) :: dim_p
      ! number of observations
      integer, intent(in) :: dim_obs
      ! pe local part of matrix hp
      real, DIMENSION(dim_obs, dim_p), intent(inout) :: hp_p
      ! matrix hph
      real, DIMENSION(dim_obs, dim_obs), intent(inout) :: hph
      call localize_covar_pdaf_c_ptr(dim_p, dim_obs, hp_p, hph)
   end subroutine f__localize_covar_pdaf

   SUBROUTINE f__localize_covar_serial_pdaf(iobs, dim_p, dim_obs, HP_p, HXY_p)
      IMPLICIT NONE
      ! Index of current observation
      INTEGER, INTENT(in) :: iobs
      ! Process-local state dimension
      INTEGER, INTENT(in) :: dim_p
      ! Number of observations
      INTEGER, INTENT(in) :: dim_obs
      ! Process-local part of matrix HP for observation iobs
      REAL, DIMENSION(dim_p), INTENT(inout) :: HP_p
      ! Process-local part of matrix HX(HX_all) for full observations
      REAL, DIMENSION(dim_obs), INTENT(inout) :: HXY_p
      call localize_covar_serial_pdaf_c_ptr(iobs, dim_p, dim_obs, HP_p, HXY_p)
   END SUBROUTINE f__localize_covar_serial_pdaf


   SUBROUTINE f__likelihood_pdaf(step, dim_obs_p, obs_p, resid, likely)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Number of observations at current time step (i.e. the size of the observation vector)
      integer, intent(in) :: dim_obs_p
      ! Vector of observations
      real, intent(in), dimension(dim_obs_p) :: obs_p
      ! Input vector holding the residual
      real, intent(in), dimension(dim_obs_p) :: resid
      ! Output value of the likelihood
      real, intent(out) :: likely
      call likelihood_pdaf_c_ptr(step, dim_obs_p, obs_p, resid, likely)
   END SUBROUTINE f__likelihood_pdaf

   SUBROUTINE f__likelihood_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer, intent(in) :: dim_obs_l
      ! Local vector of observations
      real, intent(in), dimension(dim_obs_l) :: obs_l
      ! nput vector holding the local residual
      real, intent(inout), dimension(dim_obs_l) :: resid_l
      ! Output value of the local likelihood
      real, intent(out) :: likely_l
      call likelihood_l_pdaf_c_ptr(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l)
   END SUBROUTINE f__likelihood_l_pdaf

   SUBROUTINE f__get_obs_f_pdaf(step, dim_obs_f, observation_f)
      implicit none
      ! Current time step
      integer, intent(in) :: step
      ! Size of the full observation vector
      integer, intent(in) :: dim_obs_f
      ! Full vector of synthetic observations (process-local)
      real, intent(out), dimension(dim_obs_f) :: observation_f
      call get_obs_f_pdaf_c_ptr(step, dim_obs_f, observation_f)
   END SUBROUTINE f__get_obs_f_pdaf

   SUBROUTINE f__cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vcv_p, cv_p)
      implicit none
      ! Iteration of optimization
      INTEGER, INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER, INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER, INTENT(in) :: dim_ens
      ! PE-local dimension of control vector
      INTEGER, INTENT(in) :: dim_cv_ens_p
      ! PE-local ensemble
      REAL, DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! PE-local input vector
      REAL, DIMENSION(dim_p), INTENT(in) :: Vcv_p
      ! PE-local result vector
      REAL, DIMENSION(dim_cv_ens_p), INTENT(inout) :: cv_p

      call cvt_adj_ens_pdaf_c_ptr(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vcv_p, cv_p)
   END SUBROUTINE f__cvt_adj_ens_pdaf

   SUBROUTINE f__cvt_adj_pdaf(iter, dim_p, dim_cvec, Vcv_p, cv_p)
      implicit none
      ! Iteration of optimization
      INTEGER, INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER, INTENT(in) :: dim_p
      ! Dimension of control vector
      INTEGER, INTENT(in) :: dim_cvec
      ! PE-local result vector (state vector increment)
      REAL, DIMENSION(dim_p), INTENT(in) :: Vcv_p
      ! PE-local control vector
      REAL, DIMENSION(dim_cvec), INTENT(inout) :: cv_p
      call cvt_adj_pdaf_c_ptr(iter, dim_p, dim_cvec, Vcv_p, cv_p)
   END SUBROUTINE f__cvt_adj_pdaf

   SUBROUTINE f__cvt_pdaf(iter, dim_p, dim_cvec, cv_p, Vv_p)
      implicit none
      ! Iteration of optimization
      INTEGER, INTENT(in) :: iter
      ! PE-local observation dimension
      INTEGER, INTENT(in) :: dim_p
      ! Dimension of control vector
      INTEGER, INTENT(in) :: dim_cvec
      ! PE-local control vector
      REAL, DIMENSION(dim_cvec), INTENT(in) :: cv_p
      ! PE-local result vector (state vector increment)
      REAL, DIMENSION(dim_p), INTENT(inout) :: Vv_p
      call cvt_pdaf_c_ptr(iter, dim_p, dim_cvec, cv_p, Vv_p)
   END SUBROUTINE f__cvt_pdaf

   SUBROUTINE f__cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, Vv_p)
      IMPLICIT NONE
      ! Iteration of optimization
      INTEGER, INTENT(in) :: iter
      ! PE-local dimension of state
      INTEGER, INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER, INTENT(in) :: dim_ens
      ! Dimension of control vector
      INTEGER, INTENT(in) :: dim_cvec_ens
      ! PE-local ensemble
      REAL, DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! PE-local control vector
      REAL, DIMENSION(dim_cvec_ens), INTENT(in) :: v_p
      ! PE-local state increment
      REAL, DIMENSION(dim_p), INTENT(inout) :: Vv_p
      call cvt_ens_pdaf_c_ptr(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, Vv_p)
   END SUBROUTINE f__cvt_ens_pdaf

   SUBROUTINE f__obs_op_adj_pdaf(step, dim_p, dim_obs_p, m_state_p, state_p)
      implicit none
      ! Current time step
      INTEGER, INTENT(in) :: step
      ! PE-local dimension of state
      INTEGER, INTENT(in) :: dim_p
      ! Dimension of observed state
      INTEGER, INTENT(in) :: dim_obs_p
      ! PE-local observed state
      REAL, DIMENSION(dim_obs_p), INTENT(in) :: m_state_p
      ! PE-local model state
      REAL, DIMENSION(dim_p), INTENT(inout) :: state_p
      call obs_op_adj_pdaf_c_ptr(step, dim_p, dim_obs_p, m_state_p, state_p)
   END SUBROUTINE f__obs_op_adj_pdaf

   SUBROUTINE f__obs_op_lin_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)
      implicit none
      ! Current time step
      INTEGER, INTENT(in) :: step
      ! PE-local dimension of state
      INTEGER, INTENT(in) :: dim_p
      ! Dimension of observed state
      INTEGER, INTENT(in) :: dim_obs_p
      ! PE-local model state
      REAL, DIMENSION(dim_p), INTENT(in) :: state_p
      ! PE-local observed state
      REAL, DIMENSION(dim_obs_p), INTENT(inout) :: m_state_p
      call obs_op_lin_pdaf_c_ptr(step, dim_p, dim_obs_p, state_p, m_state_p)
   END SUBROUTINE f__obs_op_lin_pdaf

   SUBROUTINE f__dist_stateinc_pdaf(dim_p, state_inc_p, first, steps)
      implicit none
      ! Dimension of PE-local state
      INTEGER, INTENT(in) :: dim_p
      ! PE-local increment of state vector
      REAL, DIMENSION(dim_p), INTENT(in) :: state_inc_p
      ! Flag for first call of each forecast
      INTEGER, INTENT(in) :: first
      ! number of time steps in forecast
      INTEGER, INTENT(in) :: steps
      call dist_stateinc_pdaf_c_ptr(dim_p, state_inc_p, first, steps)
   END SUBROUTINE f__dist_stateinc_pdaf

   SUBROUTINE f__likelihood_hyb_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, gamma, likely_l)
      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer, intent(in) :: dim_obs_l
      ! Local vector of observations
      real, dimension(dim_obs_l), intent(in) :: obs_l
      ! Hybrid weight provided by PDAF
      real, intent(in) :: gamma
      ! Input vector holding the local residual
      real, dimension(dim_obs_l), intent(inout) :: resid_l
      ! Output value of the local likelihood
      real, intent(out) :: likely_l
      call likelihood_hyb_l_pdaf_c_ptr(domain_p, step, dim_obs_l, obs_l, resid_l, gamma, likely_l)
   END SUBROUTINE f__likelihood_hyb_l_pdaf

   SUBROUTINE f__prodRinvA_hyb_l_pdaf(domain_p, step, dim_obs_l, dim_ens, obs_l, gamma, A_l, C_l)

      implicit none
      ! Index of current local analysis domain
      integer, intent(in) :: domain_p
      ! Current time step
      integer, intent(in) :: step
      ! Number of local observations at current time step (i.e. the size of the local observation vector)
      integer, intent(in) :: dim_obs_l
      ! Number of the columns in the matrix processes here. This is usually the ensemble size minus one (or the rank of the initial covariance matrix)
      integer, intent(in) :: dim_ens
      ! Local vector of observations
      real, dimension(dim_obs_l), intent(in) :: obs_l
      ! Hybrid weight provided by PDAF
      real, intent(in) :: gamma
      ! Input matrix provided by PDAF
      real, dimension(dim_obs_l, dim_ens), intent(inout) :: A_l
      ! Output matrix
      real, dimension(dim_obs_l, dim_ens), intent(out) :: C_l
      call prodRinvA_hyb_l_pdaf_c_ptr(domain_p, step, dim_obs_l, dim_ens, obs_l, gamma, A_l, C_l)
   END SUBROUTINE f__prodRinvA_hyb_l_pdaf

end module pdaf_c_f_interface