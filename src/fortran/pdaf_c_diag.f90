MODULE pdaf_c_diag
use PDAF

implicit none

contains
   SUBROUTINE c__PDAF_diag_ensmean(dim, dim_ens, state, ens, status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: ens
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_ensmean(dim, dim_ens, state, ens, status)

   END SUBROUTINE c__PDAF_diag_ensmean

   SUBROUTINE c__PDAF_diag_stddev_nompi(dim, dim_ens, state, ens, stddev,  &
      do_mean, status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: ens
      ! Standard deviation of ensemble
      REAL(c_double), INTENT(out) :: stddev
      ! Whether to compute ensemble mean
      INTEGER(c_int), INTENT(in) :: do_mean
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_stddev_nompi(dim, dim_ens, state, ens, stddev, do_mean,  &
         status)

   END SUBROUTINE c__PDAF_diag_stddev_nompi

   SUBROUTINE c__PDAF_diag_stddev(dim_p, dim_ens, state_p, ens_p, stddev_g,  &
      do_mean, comm_filter, status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! State ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Global mean standard deviation of ensemble
      REAL(c_double), INTENT(out) :: stddev_g
      ! Whether to compute ensemble mean
      INTEGER(c_int), INTENT(in) :: do_mean
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: comm_filter
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_stddev(dim_p, dim_ens, state_p, ens_p, stddev_g, do_mean,  &
         comm_filter, status)

   END SUBROUTINE c__PDAF_diag_stddev

   SUBROUTINE c__PDAF_diag_variance_nompi(dim, dim_ens, state, ens, variance,  &
      stddev, do_mean, do_stddev, status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: ens
      ! Variance state vector
      REAL(c_double), DIMENSION(dim), INTENT(out) :: variance
      ! Standard deviation of ensemble
      REAL(c_double), INTENT(out) :: stddev
      ! Whether to compute ensemble mean
      INTEGER(c_int), INTENT(in) :: do_mean
      ! Whether to compute the ensemble mean standard deviation
      INTEGER(c_int), INTENT(in) :: do_stddev
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_variance_nompi(dim, dim_ens, state, ens, variance, stddev,  &
         do_mean, do_stddev, status)

   END SUBROUTINE c__PDAF_diag_variance_nompi

   SUBROUTINE c__PDAF_diag_variance(dim_p, dim_ens, state_p, ens_p, variance_p,  &
      stddev_g, do_mean, do_stddev, comm_filter, status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! State vector
      REAL(c_double), DIMENSION(dim_p), INTENT(inout) :: state_p
      ! State ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: ens_p
      ! Variance state vector
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: variance_p
      ! Global standard deviation of ensemble
      REAL(c_double), INTENT(out) :: stddev_g
      ! Whether to compute ensemble mean
      INTEGER(c_int), INTENT(in) :: do_mean
      ! Whether to compute the ensemble mean standard deviation
      INTEGER(c_int), INTENT(in) :: do_stddev
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: comm_filter
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_variance(dim_p, dim_ens, state_p, ens_p, variance_p,  &
         stddev_g, do_mean, do_stddev, comm_filter, status)

   END SUBROUTINE c__PDAF_diag_variance

   SUBROUTINE c__PDAF_diag_rmsd_nompi(dim_p, statea_p, stateb_p, rmsd_p,  &
      status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! State vector A
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: statea_p
      ! State vector B
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: stateb_p
      ! RSMD
      REAL(c_double), INTENT(out) :: rmsd_p
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_rmsd_nompi(dim_p, statea_p, stateb_p, rmsd_p, status)

   END SUBROUTINE c__PDAF_diag_rmsd_nompi

   SUBROUTINE c__PDAF_diag_rmsd(dim_p, statea_p, stateb_p, rmsd_g, comm_filter,  &
      status) bind(c)
      ! state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! State vector A
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: statea_p
      ! State vector B
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: stateb_p
      ! Global RSMD
      REAL(c_double), INTENT(out) :: rmsd_g
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: comm_filter
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_rmsd(dim_p, statea_p, stateb_p, rmsd_g, comm_filter, status)

   END SUBROUTINE c__PDAF_diag_rmsd

   SUBROUTINE c__PDAF_diag_crps(dim_p, dim_ens, element, oens, obs, crps, reli,  &
      pot_crps, uncert, status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! index of element in full state vector
      INTEGER(c_int), INTENT(in) :: element
      ! State ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: oens
      ! Observation / truth
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: obs
      ! CRPS
      REAL(c_double), INTENT(out) :: crps
      ! Reliability
      REAL(c_double), INTENT(out) :: reli
      ! potential CRPS
      REAL(c_double), INTENT(out) :: pot_crps
      ! uncertainty
      REAL(c_double), INTENT(out) :: uncert
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_crps(dim_p, dim_ens, element, oens, obs, crps, reli,  &
         pot_crps, uncert, status)

   END SUBROUTINE c__PDAF_diag_crps

   SUBROUTINE c__PDAF_diag_crps_mpi(dim_p, dim_ens, element, oens, obs,  &
      comm_filter, mype_filter, npes_filter, crps, reli, pot_crps, uncert,  &
      status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! index of element in full state vector
      INTEGER(c_int), INTENT(in) :: element
      ! State ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens), INTENT(in) :: oens
      ! Observation / truth
      REAL(c_double), DIMENSION(dim_p), INTENT(in) :: obs
      ! MPI communicator for filter
      INTEGER(c_int), INTENT(in) :: comm_filter
      ! rank of MPI communicator
      INTEGER(c_int), INTENT(in) :: mype_filter
      ! size of MPI communicator
      INTEGER(c_int), INTENT(in) :: npes_filter
      ! CRPS
      REAL(c_double), INTENT(out) :: crps
      ! Reliability
      REAL(c_double), INTENT(out) :: reli
      ! potential CRPS
      REAL(c_double), INTENT(out) :: pot_crps
      ! uncertainty
      REAL(c_double), INTENT(out) :: uncert
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_crps_mpi(dim_p, dim_ens, element, oens, obs, comm_filter,  &
         mype_filter, npes_filter, crps, reli, pot_crps, uncert, status)

   END SUBROUTINE c__PDAF_diag_crps_mpi

   SUBROUTINE c__PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, crps,  &
      reli, resol, uncert, status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! ID of element to be used
      INTEGER(c_int), INTENT(in) :: element
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: oens
      ! State ensemble
      REAL(c_double), DIMENSION(dim), INTENT(in) :: obs
      ! CRPS
      REAL(c_double), INTENT(out) :: crps
      ! Reliability
      REAL(c_double), INTENT(out) :: reli
      ! resolution
      REAL(c_double), INTENT(out) :: resol
      ! uncertainty
      REAL(c_double), INTENT(out) :: uncert
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, crps, reli,  &
         resol, uncert, status)

   END SUBROUTINE c__PDAF_diag_CRPS_nompi

   SUBROUTINE c__PDAF_diag_effsample(dim_sample, weights, n_eff) bind(c)
      ! Sample size
      INTEGER(c_int), INTENT(in) :: dim_sample
      ! Weights of the samples
      REAL(c_double), DIMENSION(dim_sample), INTENT(in) :: weights
      ! Effecfive sample size
      REAL(c_double), INTENT(out) :: n_eff


      call PDAF_diag_effsample(dim_sample, weights, n_eff)

   END SUBROUTINE c__PDAF_diag_effsample

   SUBROUTINE c__PDAF_diag_ensstats(dim, dim_ens, element, state, ens,  &
      skewness, kurtosis, status) bind(c)
      ! PE-local state dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! ID of element to be used
      INTEGER(c_int), INTENT(in) :: element
      ! State vector
      REAL(c_double), DIMENSION(dim), INTENT(in) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: ens
      ! Skewness of ensemble
      REAL(c_double), INTENT(out) :: skewness
      ! Kurtosis of ensemble
      REAL(c_double), INTENT(out) :: kurtosis
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_ensstats(dim, dim_ens, element, state, ens, skewness,  &
         kurtosis, status)

   END SUBROUTINE c__PDAF_diag_ensstats

   SUBROUTINE c__PDAF_diag_compute_moments(dim_p, dim_ens, ens, kmax, moments,  &
      bias) bind(c)
      ! local size of the state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! number of ensemble members/samples
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! ensemble matrix
      REAL(c_double), DIMENSION(dim_p,dim_ens), INTENT(in) :: ens
      ! maximum order of central moment that is computed, maximum is 4
      INTEGER(c_int), INTENT(in) :: kmax
      ! The columns contain the moments of the ensemble
      REAL(c_double), DIMENSION(dim_p, kmax), INTENT(out) :: moments
      ! if 0 bias correction is applied (default)
      INTEGER(c_int), INTENT(in) :: bias


      call PDAF_diag_compute_moments(dim_p, dim_ens, ens, kmax, moments, bias)

   END SUBROUTINE c__PDAF_diag_compute_moments

   SUBROUTINE c__PDAF_diag_histogram(ncall, dim, dim_ens, element, state, ens,  &
      hist, delta, status) bind(c)
      ! Number of calls to routine
      INTEGER(c_int), INTENT(in) :: ncall
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Element of vector used for histogram
      INTEGER(c_int), INTENT(in) :: element
      ! State vector
      REAL(c_double), DIMENSION(dim), INTENT(in) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(in) :: ens
      ! Histogram about the state
      INTEGER(c_int), DIMENSION(dim_ens+1), INTENT(inout) :: hist
      ! deviation measure from flat histogram
      REAL(c_double), INTENT(out) :: delta
      ! Status flag (0=success)
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_diag_histogram(ncall, dim, dim_ens, element, state, ens, hist,  &
         delta, status)

   END SUBROUTINE c__PDAF_diag_histogram

   SUBROUTINE c__PDAF_diag_reliability_budget(n_times, dim_ens, dim_p, ens_p,  &
      obsvar, obs_p, budget, bias_2) bind(c)
      ! Number of time steps
      INTEGER(c_int), INTENT(in) :: n_times
      ! Number of ensemble members
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Dimension of the state vector
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble matrix over times
      REAL(c_double), DIMENSION(dim_p, dim_ens, n_times), INTENT(in) :: ens_p
      ! Squared observation error/variance at n_times
      REAL(c_double), DIMENSION(dim_p, dim_ens, n_times), INTENT(in) :: obsvar
      ! Observation vector
      REAL(c_double), DIMENSION(dim_p, n_times), INTENT(in) :: obs_p
      ! Budget term for a single time step
      REAL(c_double), DIMENSION(dim_p, n_times, 5), INTENT(out) :: budget
      ! bias^2 uses
      REAL(c_double), DIMENSION(dim_p), INTENT(out) :: bias_2


      call PDAF_diag_reliability_budget(n_times, dim_ens, dim_p, ens_p, obsvar,  &
         obs_p, budget, bias_2)

   END SUBROUTINE c__PDAF_diag_reliability_budget
END MODULE pdaf_c_diag
