module pdafomi_c_internal
use iso_c_binding, only: c_int, c_double, c_bool
use pdafomi_c, only: n_obs_omi, thisobs, thisobs_l
use pdafomi_obs_f
use pdafomi_obs_l
use PDAFomi_obs_op
use PDAFomi_dim_obs_l
use PDAFomi_obs_diag
implicit none
contains
   SUBROUTINE c__PDAFomi_set_globalobs(globalobs_in) bind(c)
      ! Input value of globalobs
      INTEGER(c_int), INTENT(in) :: globalobs_in


      call PDAFomi_set_globalobs(globalobs_in)

   END SUBROUTINE c__PDAFomi_set_globalobs

   SUBROUTINE c__PDAFomi_diag_omit_by_inno() bind(c)
      call PDAFomi_diag_omit_by_inno()

   END SUBROUTINE c__PDAFomi_diag_omit_by_inno

   SUBROUTINE c__PDAFomi_cnt_dim_obs_l(i_obs, coords_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (thisobs%ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l


      call PDAFomi_cnt_dim_obs_l(thisobs_l(i_obs), thisobs(i_obs), coords_l)

   END SUBROUTINE c__PDAFomi_cnt_dim_obs_l

   SUBROUTINE c__PDAFomi_cnt_dim_obs_l_noniso(i_obs, coords_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (thisobs%ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l


      call PDAFomi_cnt_dim_obs_l_noniso(thisobs_l(i_obs), thisobs(i_obs), coords_l)

   END SUBROUTINE c__PDAFomi_cnt_dim_obs_l_noniso

   SUBROUTINE c__PDAFomi_init_obsarrays_l(i_obs, coords_l, off_obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current water column (thisobs%ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! input: offset of current obs. in local obs. vector
      INTEGER(c_int), INTENT(inout) :: off_obs_l_all


      call PDAFomi_init_obsarrays_l(thisobs_l(i_obs), thisobs(i_obs), coords_l,  &
         off_obs_l_all)

   END SUBROUTINE c__PDAFomi_init_obsarrays_l

   SUBROUTINE c__PDAFomi_init_obsarrays_l_noniso(i_obs, coords_l, off_obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current water column (thisobs%ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! input: offset of current obs. in local obs. vector
      INTEGER(c_int), INTENT(inout) :: off_obs_l_all


      call PDAFomi_init_obsarrays_l_noniso(thisobs_l(i_obs), thisobs(i_obs),  &
         coords_l, off_obs_l_all)

   END SUBROUTINE c__PDAFomi_init_obsarrays_l_noniso

   SUBROUTINE c__PDAFomi_g2l_obs(i_obs, obs_f_all, obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Full obs. vector of current obs. for all variables
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_all
      ! Local observation vector for all variables
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_l_all


      call PDAFomi_g2l_obs(thisobs_l(i_obs), thisobs(i_obs), obs_f_all, obs_l_all)

   END SUBROUTINE c__PDAFomi_g2l_obs

   SUBROUTINE c__PDAFomi_init_obs_l(i_obs, obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Local observation vector for all variables
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_l_all


      call PDAFomi_init_obs_l(thisobs_l(i_obs), thisobs(i_obs), obs_l_all)

   END SUBROUTINE c__PDAFomi_init_obs_l

   SUBROUTINE c__PDAFomi_init_obsvar_l(i_obs, meanvar_l, cnt_obs_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Mean variance
      REAL(c_double), INTENT(inout) :: meanvar_l
      ! Observation counter
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l


      call PDAFomi_init_obsvar_l(thisobs_l(i_obs), thisobs(i_obs), meanvar_l,  &
         cnt_obs_l)

   END SUBROUTINE c__PDAFomi_init_obsvar_l

   SUBROUTINE c__PDAFomi_prodRinvA_l(i_obs, nobs_all, ncols, a_l, c_l, verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Dimension of local obs. vector (all obs. types)
      INTEGER(c_int), INTENT(in) :: nobs_all
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix (thisobs_l%dim_obs_l, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: a_l
      ! Output matrix (thisobs_l%dim_obs_l, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(out) :: c_l
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_prodRinvA_l(thisobs_l(i_obs), thisobs(i_obs), nobs_all,  &
         ncols, a_l, c_l, verbose)

   END SUBROUTINE c__PDAFomi_prodRinvA_l

   SUBROUTINE c__PDAFomi_prodRinvA_hyb_l(i_obs, nobs_all, ncols, gamma, a_l, c_l,  &
      verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Dimension of local obs. vector (all obs. types)
      INTEGER(c_int), INTENT(in) :: nobs_all
      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: ncols
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: gamma
      ! Input matrix (thisobs_l%dim_obs_l, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: a_l
      ! Output matrix (thisobs_l%dim_obs_l, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(out) :: c_l
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_prodRinvA_hyb_l(thisobs_l(i_obs), thisobs(i_obs), nobs_all,  &
         ncols, gamma, a_l, c_l, verbose)

   END SUBROUTINE c__PDAFomi_prodRinvA_hyb_l

   SUBROUTINE c__PDAFomi_likelihood_l(i_obs, resid_l_all, lhood_l, verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Input vector of residuum
      REAL(c_double), DIMENSION(:), INTENT(inout) :: resid_l_all
      ! Output vector - log likelihood
      REAL(c_double), INTENT(inout) :: lhood_l
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_likelihood_l(thisobs_l(i_obs), thisobs(i_obs), resid_l_all,  &
         lhood_l, verbose)

   END SUBROUTINE c__PDAFomi_likelihood_l

   SUBROUTINE c__PDAFomi_likelihood_hyb_l(i_obs, resid_l_all, gamma, lhood_l, verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Input vector of residuum
      REAL(c_double), DIMENSION(:), INTENT(inout) :: resid_l_all
      ! Hybrid weight
      REAL(c_double), INTENT(in) :: gamma
      ! Output vector - log likelihood
      REAL(c_double), INTENT(inout) :: lhood_l
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_likelihood_hyb_l(thisobs_l(i_obs), thisobs(i_obs),  &
         resid_l_all, gamma, lhood_l, verbose)

   END SUBROUTINE c__PDAFomi_likelihood_hyb_l

   SUBROUTINE c__PDAFomi_g2l_obs_internal(i_obs, obs_f_one, offset_obs_l_all,  &
      obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Full obs. vector of current obs. type (nobs_f_one)
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_one
      ! Offset of current observation in obs_l_all and ivar_l_all
      INTEGER(c_int), INTENT(in) :: offset_obs_l_all
      ! Local observation vector for all variables (nobs_l_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_l_all


      call PDAFomi_g2l_obs_internal(thisobs_l(i_obs), obs_f_one,  &
         offset_obs_l_all, obs_l_all)

   END SUBROUTINE c__PDAFomi_g2l_obs_internal

   SUBROUTINE c__PDAFomi_comp_dist2(i_obs, coordsa, coordsb, distance2,  &
      verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsa
      ! Coordinates of observation (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsb
      ! Squared distance
      REAL(c_double), INTENT(out) :: distance2
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_comp_dist2(thisobs(i_obs), coordsa, coordsb, distance2, verbose)

   END SUBROUTINE c__PDAFomi_comp_dist2

   SUBROUTINE c__PDAFomi_check_dist2(i_obs, coordsa, coordsb, distance2, checkdist,  &
      verbose, cnt_obs) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsa
      ! Coordinates of observation (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsb
      ! Squared distance
      REAL(c_double), INTENT(out) :: distance2
      ! Flag whether distance is within cut-off radius
      LOGICAL(c_bool), INTENT(out) :: checkdist
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Count number of local observations
      INTEGER(c_int), INTENT(inout) :: cnt_obs

      logical :: checkdist_out
      call PDAFomi_check_dist2(thisobs(i_obs), thisobs_l(i_obs), coordsa,  &
         coordsb, distance2, checkdist_out, verbose, cnt_obs)
      checkdist = checkdist_out
   END SUBROUTINE c__PDAFomi_check_dist2

   SUBROUTINE c__PDAFomi_check_dist2_noniso(i_obs, coordsa, coordsb, distance2, dists,  &
      cradius, sradius, checkdist, verbose, cnt_obs) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsa
      ! Coordinates of observation (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsb
      ! Squared distance
      REAL(c_double), INTENT(out) :: distance2
      ! Vector of distance in each coordinate direction
      REAL(c_double), DIMENSION(:), INTENT(inout) :: dists
      ! Directional cut-off radius
      REAL(c_double), INTENT(out) :: cradius
      ! Directional support radius
      REAL(c_double), INTENT(inout) :: sradius
      ! Flag whether distance is within cut-off radius
      LOGICAL(c_bool), INTENT(out) :: checkdist
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose
      ! Count number of local observations
      INTEGER(c_int), INTENT(inout) :: cnt_obs
      
      logical :: checkdist_out

      call PDAFomi_check_dist2_noniso(thisobs(i_obs), thisobs_l(i_obs),  &
         coordsa, coordsb, distance2, dists, cradius, sradius, checkdist_out,  &
         verbose, cnt_obs)

      checkdist = checkdist_out
   END SUBROUTINE c__PDAFomi_check_dist2_noniso

   SUBROUTINE c__PDAFomi_weights_l(verbose, nobs_l, ncols, locweight, cradius,  &
      sradius, mata, ivar_obs_l, dist_l, weight_l) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Number of local observations
      INTEGER(c_int), INTENT(in) :: nobs_l
      !
      INTEGER(c_int), INTENT(in) :: ncols
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! Localization cut-off radius
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! support radius for weight functions
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      !
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: mata
      ! Local vector of inverse obs. variances (nobs_l)
      REAL(c_double), DIMENSION(:), INTENT(in) :: ivar_obs_l
      ! Local vector of obs. distances (nobs_l)
      REAL(c_double), DIMENSION(:), INTENT(in) :: dist_l
      ! Output: vector of weights
      REAL(c_double), DIMENSION(:), INTENT(out) :: weight_l


      call PDAFomi_weights_l(verbose, nobs_l, ncols, locweight, cradius,  &
         sradius, mata, ivar_obs_l, dist_l, weight_l)

   END SUBROUTINE c__PDAFomi_weights_l

   SUBROUTINE c__PDAFomi_weights_l_sgnl(verbose, nobs_l, ncols, locweight,  &
      cradius, sradius, mata, ivar_obs_l, dist_l, weight_l) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Number of local observations
      INTEGER(c_int), INTENT(in) :: nobs_l
      !
      INTEGER(c_int), INTENT(in) :: ncols
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! Localization cut-off radius
      REAL(c_double), INTENT(in) :: cradius
      ! support radius for weight functions
      REAL(c_double), INTENT(in) :: sradius
      !
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: mata
      ! Local vector of inverse obs. variances (nobs_l)
      REAL(c_double), DIMENSION(:), INTENT(in) :: ivar_obs_l
      ! Local vector of obs. distances (nobs_l)
      REAL(c_double), DIMENSION(:), INTENT(in) :: dist_l
      ! Output: vector of weights
      REAL(c_double), DIMENSION(:), INTENT(out) :: weight_l


      call PDAFomi_weights_l_sgnl(verbose, nobs_l, ncols, locweight, cradius,  &
         sradius, mata, ivar_obs_l, dist_l, weight_l)

   END SUBROUTINE c__PDAFomi_weights_l_sgnl

   SUBROUTINE c__PDAFomi_omit_by_inno_l(i_obs, inno_l, obs_l_all, obsid, cnt_all,  &
      verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Input vector of observation innovation
      REAL(c_double), DIMENSION(:), INTENT(in) :: inno_l
      ! Input vector of local observations
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_l_all
      ! ID of observation type
      INTEGER(c_int), INTENT(in) :: obsid
      ! Count of omitted observation over all types
      INTEGER(c_int), INTENT(inout) :: cnt_all
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_omit_by_inno_l(thisobs_l(i_obs), thisobs(i_obs), inno_l,  &
         obs_l_all, obsid, cnt_all, verbose)

   END SUBROUTINE c__PDAFomi_omit_by_inno_l

   SUBROUTINE c__PDAFomi_obsstats_l(screen) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAFomi_obsstats_l(screen)

   END SUBROUTINE c__PDAFomi_obsstats_l

   SUBROUTINE c__PDAFomi_dealloc() bind(c)
      call PDAFomi_dealloc()

   END SUBROUTINE c__PDAFomi_dealloc

   SUBROUTINE c__PDAFomi_ocoord_all(ncoord, oc_all) bind(c)
      ! Number of coordinate directions
      INTEGER(c_int), INTENT(in) :: ncoord
      ! Array of observation coordinates size(ncoord, dim_obs)
      REAL(c_double), DIMENSION(:,:), INTENT(out) :: oc_all


      call PDAFomi_ocoord_all(ncoord, oc_all)

   END SUBROUTINE c__PDAFomi_ocoord_all

   SUBROUTINE c__PDAFomi_local_weight(wtype, rtype, cradius, sradius, distance,  &
      nrows, ncols, a, var_obs, weight, verbose) bind(c)
      ! Type of weight function
      INTEGER(c_int), INTENT(in) :: wtype
      ! Type of regulated weighting
      INTEGER(c_int), INTENT(in) :: rtype
      ! Cut-off radius
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius
      REAL(c_double), INTENT(in) :: sradius
      ! Distance to observation
      REAL(c_double), INTENT(in) :: distance
      ! Number of rows in matrix A
      INTEGER(c_int), INTENT(in) :: nrows
      ! Number of columns in matrix A
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix
      REAL(c_double), DIMENSION(nrows, ncols), INTENT(in) :: a
      ! Observation variance
      REAL(c_double), INTENT(in) :: var_obs
      ! Weights
      REAL(c_double), INTENT(out) :: weight
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_local_weight(wtype, rtype, cradius, sradius, distance,  &
         nrows, ncols, a, var_obs, weight, verbose)

   END SUBROUTINE c__PDAFomi_local_weight

   SUBROUTINE c__PDAFomi_check_dist2_loop(i_obs, coordsa, cnt_obs, mode) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsa
      ! Count number of local observations
      INTEGER(c_int), INTENT(inout) :: cnt_obs
      ! 1: count local observations
      INTEGER(c_int), INTENT(in) :: mode


      call PDAFomi_check_dist2_loop(thisobs_l(i_obs), thisobs(i_obs), coordsa,  &
         cnt_obs, mode)

   END SUBROUTINE c__PDAFomi_check_dist2_loop

   SUBROUTINE c__PDAFomi_check_dist2_noniso_loop(i_obs, coordsa, cnt_obs, mode) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain (ncoord)
      REAL(c_double), DIMENSION(:), INTENT(in) :: coordsa
      ! Count number of local observations
      INTEGER(c_int), INTENT(inout) :: cnt_obs
      ! 1: count local observations
      INTEGER(c_int), INTENT(in) :: mode


      call PDAFomi_check_dist2_noniso_loop(thisobs_l(i_obs), thisobs(i_obs),  &
         coordsa, cnt_obs, mode)

   END SUBROUTINE c__PDAFomi_check_dist2_noniso_loop

   SUBROUTINE c__PDAFomi_obs_op_gatheronly(i_obs, state_p, obs_f_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(in) :: state_p
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_f_all


      call PDAFomi_obs_op_gatheronly(thisobs(i_obs), state_p, obs_f_all)

   END SUBROUTINE c__PDAFomi_obs_op_gatheronly

   SUBROUTINE c__PDAFomi_obs_op_adj_gatheronly(i_obs, obs_f_all, state_p) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(in) :: state_p


      call PDAFomi_obs_op_adj_gatheronly(thisobs(i_obs), obs_f_all, state_p)

   END SUBROUTINE c__PDAFomi_obs_op_adj_gatheronly

   SUBROUTINE c__PDAFomi_init_obs_f(i_obs, dim_obs_f, obsstate_f, offset) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Dimension of full observed state (all observed fields)
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector (dim_obs_f)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obsstate_f
      ! input: offset of module-type observations in obsstate_f
      INTEGER(c_int), INTENT(inout) :: offset


      call PDAFomi_init_obs_f(thisobs(i_obs), dim_obs_f, obsstate_f, offset)

   END SUBROUTINE c__PDAFomi_init_obs_f

   SUBROUTINE c__PDAFomi_init_obsvars_f(i_obs, dim_obs_f, var_f, offset) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Dimension of full observed state (all observed fields)
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! Full vector of observation variances (dim_obs_f)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: var_f
      ! input: offset of module-type observations in obsstate_f
      INTEGER(c_int), INTENT(inout) :: offset


      call PDAFomi_init_obsvars_f(thisobs(i_obs), dim_obs_f, var_f, offset)

   END SUBROUTINE c__PDAFomi_init_obsvars_f

   SUBROUTINE c__PDAFomi_init_obsvar_f(i_obs, meanvar, cnt_obs) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Mean variance
      REAL(c_double), INTENT(inout) :: meanvar
      ! Observation counter
      INTEGER(c_int), INTENT(inout) :: cnt_obs


      call PDAFomi_init_obsvar_f(thisobs(i_obs), meanvar, cnt_obs)

   END SUBROUTINE c__PDAFomi_init_obsvar_f

   SUBROUTINE c__PDAFomi_prodRinvA(i_obs, ncols, a_p, c_p) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of columns in A_p and C_p
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix (nobs_f, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(in) :: a_p
      ! Output matrix (nobs_f, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(out) :: c_p


      call PDAFomi_prodRinvA(thisobs(i_obs), ncols, a_p, c_p)

   END SUBROUTINE c__PDAFomi_prodRinvA

   SUBROUTINE c__PDAFomi_likelihood(i_obs, resid, lhood) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Input vector of residuum
      REAL(c_double), DIMENSION(:), INTENT(in) :: resid
      ! Output vector - log likelihood
      REAL(c_double), INTENT(inout) :: lhood


      call PDAFomi_likelihood(thisobs(i_obs), resid, lhood)

   END SUBROUTINE c__PDAFomi_likelihood

   SUBROUTINE c__PDAFomi_add_obs_error(i_obs, nobs_all, matc) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of observations
      INTEGER(c_int), INTENT(in) :: nobs_all
      ! Input/Output matrix (nobs_f, rank)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: matc


      call PDAFomi_add_obs_error(thisobs(i_obs), nobs_all, matc)

   END SUBROUTINE c__PDAFomi_add_obs_error

   SUBROUTINE c__PDAFomi_init_obscovar(i_obs, nobs_all, covar, isdiag) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of observations
      INTEGER(c_int), INTENT(in) :: nobs_all
      ! Input/Output matrix (nobs_all, nobs_all)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: covar
      ! Whether matrix R is diagonal
      LOGICAL(c_bool), INTENT(out) :: isdiag

      logical :: isdiag_out

      call PDAFomi_init_obscovar(thisobs(i_obs), nobs_all, covar, isdiag_out)
      isdiag = isdiag_out
   END SUBROUTINE c__PDAFomi_init_obscovar

   SUBROUTINE c__PDAFomi_init_obserr_f(i_obs, obserr_f) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Full vector of observation errors
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obserr_f


      call PDAFomi_init_obserr_f(thisobs(i_obs), obserr_f)

   END SUBROUTINE c__PDAFomi_init_obserr_f

   SUBROUTINE c__PDAFomi_get_local_ids_obs_f(dim_obs_g, lradius, oc_f, cnt_lim,  &
      id_lim, disttype, domainsize) bind(c)
      ! Global full number of observations
      INTEGER(c_int), INTENT(in) :: dim_obs_g
      ! Localization radius (used is a constant one here)
      REAL(c_double), INTENT(in) :: lradius
      ! observation coordinates (radians), row 1: lon, 2: lat
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: oc_f
      ! Number of full observation for local process domain
      INTEGER(c_int), INTENT(out) :: cnt_lim
      ! Indices of process-local full obs. in global full vector
      INTEGER(c_int), DIMENSION(:), INTENT(out) :: id_lim
      ! type of distance computation
      INTEGER(c_int), INTENT(in) :: disttype
      ! Global size of model domain
      REAL(c_double), DIMENSION(:), INTENT(in) :: domainsize


      call PDAFomi_get_local_ids_obs_f(dim_obs_g, lradius, oc_f, cnt_lim,  &
         id_lim, disttype, domainsize)

   END SUBROUTINE c__PDAFomi_get_local_ids_obs_f

   SUBROUTINE c__PDAFomi_limit_obs_f(i_obs, offset, obs_f_one, obs_f_lim) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! offset of this observation in obs_f_lim
      INTEGER(c_int), INTENT(in) :: offset
      ! Global full observation vector (nobs_f)
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_one
      ! full observation vector for process domains (nobs_lim)
      REAL(c_double), DIMENSION(:), INTENT(out) :: obs_f_lim


      call PDAFomi_limit_obs_f(thisobs(i_obs), offset, obs_f_one, obs_f_lim)

   END SUBROUTINE c__PDAFomi_limit_obs_f

   SUBROUTINE c__PDAFomi_gather_dim_obs_f(dim_obs_p, dim_obs_f) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_f


      call PDAFomi_gather_dim_obs_f(dim_obs_p, dim_obs_f)

   END SUBROUTINE c__PDAFomi_gather_dim_obs_f

   SUBROUTINE c__PDAFomi_gather_obs_f_flex(dim_obs_p, obs_p, obs_f, status) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local vector
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_p
      ! Full gathered vector
      REAL(c_double), DIMENSION(:), INTENT(out) :: obs_f
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status


      call PDAFomi_gather_obs_f_flex(dim_obs_p, obs_p, obs_f, status)

   END SUBROUTINE c__PDAFomi_gather_obs_f_flex

   SUBROUTINE c__PDAFomi_gather_obs_f2_flex(dim_obs_p, coords_p, coords_f,  &
      nrows, status) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! PE-local array
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords_p
      ! Full gathered array
      REAL(c_double), DIMENSION(:,:), INTENT(out) :: coords_f
      ! Number of rows in array
      INTEGER(c_int), INTENT(in) :: nrows
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status


      call PDAFomi_gather_obs_f2_flex(dim_obs_p, coords_p, coords_f, nrows, status)

   END SUBROUTINE c__PDAFomi_gather_obs_f2_flex

   SUBROUTINE c__PDAFomi_omit_by_inno(i_obs, inno_f, obs_f_all, obsid, cnt_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Input vector of observation innovation
      REAL(c_double), DIMENSION(:), INTENT(in) :: inno_f
      ! Input vector of local observations
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_all
      ! ID of observation type
      INTEGER(c_int), INTENT(in) :: obsid
      ! Count of omitted observation over all types
      INTEGER(c_int), INTENT(inout) :: cnt_all


      call PDAFomi_omit_by_inno(thisobs(i_obs), inno_f, obs_f_all, obsid, cnt_all)

   END SUBROUTINE c__PDAFomi_omit_by_inno

   SUBROUTINE c__PDAFomi_obsstats(screen) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAFomi_obsstats(screen)

   END SUBROUTINE c__PDAFomi_obsstats

   SUBROUTINE c__PDAFomi_gather_obsdims() bind(c)
      call PDAFomi_gather_obsdims()

   END SUBROUTINE c__PDAFomi_gather_obsdims
end module pdafomi_c_internal
