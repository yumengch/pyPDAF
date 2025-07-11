module pdafomi_c_diag
implicit none

contains
   SUBROUTINE c__PDAFomi_diag_dimobs(dim_obs_ptr) bind(c)
      ! Pointer to observation dimensions
      INTEGER(c_int), POINTER, DIMENSION(:), INTENT(inout) :: dim_obs_ptr


      call PDAFomi_diag_dimobs(dim_obs_ptr)
   END SUBROUTINE c__PDAFomi_diag_dimobs

   SUBROUTINE c__PDAFomi_diag_get_HX(id_obs, dim_obs_diag, hx_p_ptr) bind(c)
      ! Index of observation type to return
      INTEGER(c_int), INTENT(in) :: id_obs
      ! Observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_diag
      ! Pointer to observed ensemble mean
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: hx_p_ptr


      call PDAFomi_diag_get_HX(id_obs, dim_obs_diag, hx_p_ptr)
   END SUBROUTINE c__PDAFomi_diag_get_HX

   SUBROUTINE c__PDAFomi_diag_get_HXmean(id_obs, dim_obs_diag,  &
      hxmean_p_ptr) bind(c)
      ! Index of observation type to return
      INTEGER(c_int), INTENT(in) :: id_obs
      ! Observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_diag
      ! Pointer to observed ensemble mean
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: hxmean_p_ptr


      call PDAFomi_diag_get_HXmean(id_obs, dim_obs_diag, hxmean_p_ptr)
   END SUBROUTINE c__PDAFomi_diag_get_HXmean

   SUBROUTINE c__PDAFomi_diag_get_ivar(id_obs, dim_obs_diag, ivar_ptr) bind(c)
      ! Index of observation type to return
      INTEGER(c_int), INTENT(in) :: id_obs
      ! Observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_diag
      ! Pointer to inverse observation error variances
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: ivar_ptr


      call PDAFomi_diag_get_ivar(id_obs, dim_obs_diag, ivar_ptr)
   END SUBROUTINE c__PDAFomi_diag_get_ivar

   SUBROUTINE c__PDAFomi_diag_get_obs(id_obs, dim_obs_diag, ncoord, obs_p_ptr,  &
      ocoord_p_ptr) bind(c)
      ! Index of observation type to return
      INTEGER(c_int), INTENT(in) :: id_obs
      ! Observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_diag
      ! Number of observation dimensions
      INTEGER(c_int), INTENT(out) :: ncoord
      ! Pointer to observation vector
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: obs_p_ptr
      ! Pointer to Coordinate array
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: ocoord_p_ptr


      call PDAFomi_diag_get_obs(id_obs, dim_obs_diag, ncoord, obs_p_ptr,  &
         ocoord_p_ptr)
   END SUBROUTINE c__PDAFomi_diag_get_obs

   SUBROUTINE c__PDAFomi_diag_nobstypes(nobs) bind(c)
      ! Number of observation types
      INTEGER(c_int), INTENT(inout) :: nobs


      call PDAFomi_diag_nobstypes(nobs)

   END SUBROUTINE c__PDAFomi_diag_nobstypes

   SUBROUTINE c__PDAFomi_diag_obs_rmsd(nobs, rmsd_pointer, verbose) bind(c)
      ! Number of observation types
      INTEGER(c_int), INTENT(inout) :: nobs
      ! Vector of RMSD values
      REAL(c_double), POINTER, DIMENSION(:), INTENT(inout) :: rmsd_pointer
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_diag_obs_rmsd(nobs, rmsd_pointer, verbose)
   END SUBROUTINE c__PDAFomi_diag_obs_rmsd

   SUBROUTINE c__PDAFomi_diag_stats(nobs, obsstats_ptr, verbose) bind(c)
      ! Number of observation types
      INTEGER(c_int), INTENT(inout) :: nobs
      ! Array of observation statistics
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(inout) :: obsstats_ptr
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAFomi_diag_stats(nobs, obsstats_ptr, verbose)
   END SUBROUTINE c__PDAFomi_diag_stats
end module pdafomi_c_diag