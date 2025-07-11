module pdafomi_c
implicit none
contains
   SUBROUTINE c__PDAFomi_check_error(flag) bind(c)
      ! Error flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAFomi_check_error(flag)

   END SUBROUTINE c__PDAFomi_check_error

   SUBROUTINE c__PDAFomi_gather_obs(i_obs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p,  &
      ncoord, lradius, dim_obs_f) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of process-local observation
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Vector of process-local observations
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_p
      ! Vector of process-local inverse observation error variance
      REAL(c_double), DIMENSION(:), INTENT(in) :: ivar_obs_p
      ! Array of process-local observation coordinates
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: ocoord_p
      ! Number of rows of coordinate array
      INTEGER(c_int), INTENT(in) :: ncoord
      ! Localization radius (the maximum radius used in this process domain)
      REAL(c_double), INTENT(in) :: lradius
      ! Full number of observations
      INTEGER(c_int), INTENT(out) :: dim_obs_f


      call PDAFomi_gather_obs(thisobs(i_obs), dim_obs_p, obs_p, ivar_obs_p,  &
         ocoord_p, ncoord, lradius, dim_obs_f)

   END SUBROUTINE c__PDAFomi_gather_obs

   SUBROUTINE c__PDAFomi_gather_obsstate(i_obs, obsstate_p, obsstate_f) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Vector of process-local observed state
      REAL(c_double), DIMENSION(:), INTENT(in) :: obsstate_p
      ! Full observed vector for all types
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obsstate_f


      call PDAFomi_gather_obsstate(thisobs(i_obs), obsstate_p, obsstate_f)

   END SUBROUTINE c__PDAFomi_gather_obsstate

   SUBROUTINE c__PDAFomi_get_interp_coeff_tri(gpc, oc, icoeff) bind(c)
      ! Coordinates of grid points; dim(3,2)
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: gpc
      ! Coordinates of observation; dim(2)
      REAL(c_double), DIMENSION(:), INTENT(in) :: oc
      ! Interpolation coefficients; dim(3)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: icoeff


      call PDAFomi_get_interp_coeff_tri(gpc, oc, icoeff)

   END SUBROUTINE c__PDAFomi_get_interp_coeff_tri

   SUBROUTINE c__PDAFomi_get_interp_coeff_lin1D(gpc, oc, icoeff) bind(c)
      ! Coordinates of grid points (dim=2)
      REAL(c_double), DIMENSION(:), INTENT(in) :: gpc
      ! Coordinates of observation
      REAL(c_double), INTENT(in) :: oc
      ! Interpolation coefficients (dim=2)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: icoeff


      call PDAFomi_get_interp_coeff_lin1D(gpc, oc, icoeff)

   END SUBROUTINE c__PDAFomi_get_interp_coeff_lin1D

   SUBROUTINE c__PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc,  &
      icoeff) bind(c)
      ! Length of icoeff
      INTEGER(c_int), INTENT(in) :: num_gp
      ! Number of dimensions in interpolation
      INTEGER(c_int), INTENT(in) :: n_dim
      ! Coordinates of grid points
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: gpc
      ! Coordinates of observation
      REAL(c_double), DIMENSION(:), INTENT(in) :: oc
      ! Interpolation coefficients (num_gp)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: icoeff


      call PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)

   END SUBROUTINE c__PDAFomi_get_interp_coeff_lin

   SUBROUTINE c__PDAFomi_init_dim_obs_l_iso(i_obs, coords_l, locweight, cradius,  &
      sradius, cnt_obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight
      ! Localization cut-off radius
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius of localization function
      REAL(c_double), INTENT(in) :: sradius
      ! Local dimension of current observation vector
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l_all


      call PDAFomi_init_dim_obs_l_iso(thisobs_l(i_obs), thisobs(i_obs),  &
         coords_l, locweight, cradius, sradius, cnt_obs_l_all)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_iso

   SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso(i_obs, coords_l, locweight, cradius,  &
      sradius, cnt_obs_l_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Local dimension of current observation vector
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l_all


      call PDAFomi_init_dim_obs_l_noniso(thisobs_l(i_obs), thisobs(i_obs),  &
         coords_l, locweight, cradius, sradius, cnt_obs_l_all)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso

   SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_locweights(i_obs, coords_l, locweights,  &
      cradius, sradius, cnt_obs_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! Types of localization function
      INTEGER(c_int), DIMENSION(:), INTENT(in) :: locweights
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Local dimension of current observation vector
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l


      call PDAFomi_init_dim_obs_l_noniso_locweights(thisobs_l(i_obs),  &
         thisobs(i_obs), coords_l, locweights, cradius, sradius, cnt_obs_l)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_locweights

   SUBROUTINE c__PDAFomi_obs_op_gridpoint(i_obs, state_p, obs_f_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(in) :: state_p
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_f_all


      call PDAFomi_obs_op_gridpoint(thisobs(i_obs), state_p, obs_f_all)

   END SUBROUTINE c__PDAFomi_obs_op_gridpoint

   SUBROUTINE c__PDAFomi_obs_op_gridavg(i_obs, nrows, state_p, obs_f_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(in) :: state_p
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_f_all


      call PDAFomi_obs_op_gridavg(thisobs(i_obs), nrows, state_p, obs_f_all)

   END SUBROUTINE c__PDAFomi_obs_op_gridavg

   SUBROUTINE c__PDAFomi_obs_op_interp_lin(i_obs, nrows, state_p, obs_f_all) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(in) :: state_p
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: obs_f_all


      call PDAFomi_obs_op_interp_lin(thisobs(i_obs), nrows, state_p, obs_f_all)

   END SUBROUTINE c__PDAFomi_obs_op_interp_lin

   SUBROUTINE c__PDAFomi_obs_op_adj_gridpoint(i_obs, obs_f_all, state_p) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state_p


      call PDAFomi_obs_op_adj_gridpoint(thisobs(i_obs), obs_f_all, state_p)

   END SUBROUTINE c__PDAFomi_obs_op_adj_gridpoint

   SUBROUTINE c__PDAFomi_obs_op_adj_gridavg(i_obs, nrows, obs_f_all, state_p) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state_p


      call PDAFomi_obs_op_adj_gridavg(thisobs(i_obs), nrows, obs_f_all, state_p)

   END SUBROUTINE c__PDAFomi_obs_op_adj_gridavg

   SUBROUTINE c__PDAFomi_obs_op_adj_interp_lin(i_obs, nrows, obs_f_all, state_p) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), DIMENSION(:), INTENT(in) :: obs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state_p


      call PDAFomi_obs_op_adj_interp_lin(thisobs(i_obs), nrows, obs_f_all, state_p)

   END SUBROUTINE c__PDAFomi_obs_op_adj_interp_lin

   SUBROUTINE c__PDAFomi_observation_localization_weights(i_obs, ncols, a_l, weight,  &
      verbose) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Rank of initial covariance matrix
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix (thisobs_l%dim_obs_l, ncols)
      REAL(c_double), DIMENSION(:, :), INTENT(in) :: a_l
      ! > Localization weights
      REAL(c_double), DIMENSION(thisobs_l%dim_obs_l), INTENT(out) :: weight
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_observation_localization_weights(thisobs_l(i_obs),  &
         thisobs(i_obs), ncols, a_l, weight, verbose)

   END SUBROUTINE c__PDAFomi_observation_localization_weights

   SUBROUTINE c__PDAFomi_set_debug_flag(debugval) bind(c)
      ! Value for debugging flag
      INTEGER(c_int), INTENT(in) :: debugval


      call PDAFomi_set_debug_flag(debugval)

   END SUBROUTINE c__PDAFomi_set_debug_flag

   SUBROUTINE c__PDAFomi_set_dim_obs_l(i_obs, cnt_obs_l_all, cnt_obs_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Local dimension of observation vector over all obs. types
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l_all
      ! Local dimension of single observation type vector
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l


      call PDAFomi_set_dim_obs_l(thisobs_l(i_obs), thisobs(i_obs),  &
         cnt_obs_l_all, cnt_obs_l)

   END SUBROUTINE c__PDAFomi_set_dim_obs_l

   SUBROUTINE c__PDAFomi_set_localization(i_obs, cradius, sradius, locweight) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Localization cut-off radius
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius of localization function
      REAL(c_double), INTENT(in) :: sradius
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight


      call PDAFomi_set_localization(thisobs_l(i_obs), cradius, sradius, locweight)

   END SUBROUTINE c__PDAFomi_set_localization

   SUBROUTINE c__PDAFomi_set_localization_noniso(i_obs, nradii, cradius, sradius,  &
      locweight, locweight_v) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of radii to consider for localization
      INTEGER(c_int), INTENT(in) :: nradii
      ! Localization cut-off radius
      REAL(c_double), DIMENSION(nradii), INTENT(in) :: cradius
      ! Support radius of localization function
      REAL(c_double), DIMENSION(nradii), INTENT(in) :: sradius
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight
      ! Type of localization function in vertical direction (only for nradii=3)
      INTEGER(c_int), INTENT(in) :: locweight_v


      call PDAFomi_set_localization_noniso(thisobs_l(i_obs), nradii, cradius,  &
         sradius, locweight, locweight_v)

   END SUBROUTINE c__PDAFomi_set_localization_noniso

   SUBROUTINE c__PDAFomi_set_localize_covar_iso(i_obs, dim, ncoords, coords,  &
      locweight, cradius, sradius) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! number of coordinate directions
      INTEGER(c_int), INTENT(in) :: ncoords
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! localization radius
      REAL(c_double), INTENT(in) :: cradius
      ! support radius for weight functions
      REAL(c_double), INTENT(in) :: sradius


      call PDAFomi_set_localize_covar_iso(thisobs(i_obs), dim, ncoords, coords,  &
         locweight, cradius, sradius)

   END SUBROUTINE c__PDAFomi_set_localize_covar_iso

   SUBROUTINE c__PDAFomi_set_localize_covar_noniso(i_obs, dim, ncoords, coords,  &
      locweight, cradius, sradius) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! number of coordinate directions
      INTEGER(c_int), INTENT(in) :: ncoords
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius


      call PDAFomi_set_localize_covar_noniso(thisobs(i_obs), dim, ncoords,  &
         coords, locweight, cradius, sradius)

   END SUBROUTINE c__PDAFomi_set_localize_covar_noniso

   SUBROUTINE c__PDAFomi_set_localize_covar_noniso_locweights(i_obs, dim, ncoords,  &
      coords, locweights, cradius, sradius) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! number of coordinate directions
      INTEGER(c_int), INTENT(in) :: ncoords
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Types of localization function
      INTEGER(c_int), DIMENSION(:), INTENT(in) :: locweights
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius


      call PDAFomi_set_localize_covar_noniso_locweights(thisobs(i_obs), dim,  &
         ncoords, coords, locweights, cradius, sradius)

   END SUBROUTINE c__PDAFomi_set_localize_covar_noniso_locweights

   SUBROUTINE c__PDAFomi_set_obs_diag(diag) bind(c)
      ! Value for observation diagnostics mode
      INTEGER(c_int), INTENT(in) :: diag


      call PDAFomi_set_obs_diag(diag)

   END SUBROUTINE c__PDAFomi_set_obs_diag

   SUBROUTINE c__PDAFomi_set_domain_limits(lim_coords) bind(c)
      ! geographic coordinate array (1: longitude, 2: latitude)
      REAL(c_double), DIMENSION(2,2), INTENT(in) :: lim_coords


      call PDAFomi_set_domain_limits(lim_coords)

   END SUBROUTINE c__PDAFomi_set_domain_limits

   SUBROUTINE c__PDAFomi_get_domain_limits_unstr(npoints_p, coords_p) bind(c)
      ! number of process-local grid points
      INTEGER(c_int), INTENT(in) :: npoints_p
      ! geographic coordinate array (row 1: longitude, 2: latitude)
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords_p


      call PDAFomi_get_domain_limits_unstr(npoints_p, coords_p)

   END SUBROUTINE c__PDAFomi_get_domain_limits_unstr

   SUBROUTINE c__PDAFomi_store_obs_l_index(i_obs, idx, id_obs_l, distance, cradius_l,  &
      sradius_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Element of local observation array to be filled
      INTEGER(c_int), INTENT(in) :: idx
      ! Index of local observation in full observation array
      INTEGER(c_int), INTENT(in) :: id_obs_l
      ! Distance between local analysis domain and observation
      REAL(c_double), INTENT(in) :: distance
      ! cut-off radius for this local observation
      REAL(c_double), INTENT(in) :: cradius_l
      ! support radius for this local observation
      REAL(c_double), INTENT(in) :: sradius_l


      call PDAFomi_store_obs_l_index(thisobs_l(i_obs), idx, id_obs_l, distance,  &
         cradius_l, sradius_l)

   END SUBROUTINE c__PDAFomi_store_obs_l_index

   SUBROUTINE c__PDAFomi_store_obs_l_index_vdist(i_obs, idx, id_obs_l, distance,  &
      cradius_l, sradius_l, vdist) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Element of local observation array to be filled
      INTEGER(c_int), INTENT(in) :: idx
      ! Index of local observation in full observation array
      INTEGER(c_int), INTENT(in) :: id_obs_l
      ! Distance between local analysis domain and observation
      REAL(c_double), INTENT(in) :: distance
      ! cut-off radius for this local observation
      REAL(c_double), INTENT(in) :: cradius_l
      ! support radius for this local observation
      REAL(c_double), INTENT(in) :: sradius_l
      ! support radius in vertical direction for 2+1D factorized localization
      REAL(c_double), INTENT(in) :: vdist


      call PDAFomi_store_obs_l_index_vdist(thisobs_l(i_obs), idx, id_obs_l,  &
         distance, cradius_l, sradius_l, vdist)

   END SUBROUTINE c__PDAFomi_store_obs_l_index_vdist
end module pdafomi_c