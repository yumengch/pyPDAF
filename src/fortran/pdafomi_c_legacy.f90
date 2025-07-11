module pdafomi_c_legacy
implicit none
contains
   SUBROUTINE c__PDAFomi_localize_covar_iso(i_obs, dim, locweight, cradius, sradius,  &
      coords, hp, hph) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! localization radius
      REAL(c_double), INTENT(in) :: cradius
      ! support radius for weight functions
      REAL(c_double), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Matrix HP, dimension (nobs, dim)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hp
      ! Matrix HPH, dimension (nobs, nobs)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hph


      call PDAFomi_localize_covar_iso(thisobs(i_obs), dim, locweight, cradius,  &
         sradius, coords, hp, hph)

   END SUBROUTINE c__PDAFomi_localize_covar_iso

   SUBROUTINE c__PDAFomi_localize_covar_noniso_locweights(i_obs, dim, locweights,  &
      cradius, sradius, coords, hp, hph) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Types of localization function
      INTEGER(c_int), DIMENSION(:), INTENT(in) :: locweights
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Matrix HP, dimension (nobs, dim)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hp
      ! Matrix HPH, dimension (nobs, nobs)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hph


      call PDAFomi_localize_covar_noniso_locweights(thisobs(i_obs), dim,  &
         locweights, cradius, sradius, coords, hp, hph)

   END SUBROUTINE c__PDAFomi_localize_covar_noniso_locweights

   SUBROUTINE c__PDAFomi_localize_covar_noniso(i_obs, dim, locweight, cradius,  &
      sradius, coords, hp, hph) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Matrix HP, dimension (nobs, dim)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hp
      ! Matrix HPH, dimension (nobs, nobs)
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: hph


      call PDAFomi_localize_covar_noniso(thisobs(i_obs), dim, locweight,  &
         cradius, sradius, coords, hp, hph)

   END SUBROUTINE c__PDAFomi_localize_covar_noniso

   SUBROUTINE c__PDAFomi_localize_covar_serial_iso(i_obs, iobs_all, dim, dim_obs,  &
      locweight, cradius, sradius, coords, hp, hxy) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Index of current observation
      INTEGER(c_int), INTENT(in) :: iobs_all
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Overall full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! localization radius
      REAL(c_double), INTENT(in) :: cradius
      ! support radius for weight functions
      REAL(c_double), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Vector HP, dimension (dim)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hp
      ! Matrix HXY, dimension (nobs)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hxy


      call PDAFomi_localize_covar_serial_iso(thisobs(i_obs), iobs_all, dim,  &
         dim_obs, locweight, cradius, sradius, coords, hp, hxy)

   END SUBROUTINE c__PDAFomi_localize_covar_serial_iso

   SUBROUTINE c__PDAFomi_localize_covar_serial_noniso_locweights(i_obs, iobs_all, dim,  &
      dim_obs, locweights, cradius, sradius, coords, hp, hxy) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Index of current observation
      INTEGER(c_int), INTENT(in) :: iobs_all
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Overall full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Types of localization function
      INTEGER(c_int), DIMENSION(:), INTENT(in) :: locweights
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Vector HP, dimension (dim)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hp
      ! Matrix HXY, dimension (nobs)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hxy


      call PDAFomi_localize_covar_serial_noniso_locweights(thisobs(i_obs),  &
         iobs_all, dim, dim_obs, locweights, cradius, sradius, coords, hp, hxy)

   END SUBROUTINE c__PDAFomi_localize_covar_serial_noniso_locweights

   SUBROUTINE c__PDAFomi_localize_covar_serial_noniso(i_obs, iobs_all, dim, dim_obs,  &
      locweight, cradius, sradius, coords, hp, hxy) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Index of current observation
      INTEGER(c_int), INTENT(in) :: iobs_all
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim
      ! Overall full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! Vector of localization cut-off radii
      REAL(c_double), DIMENSION(:), INTENT(in) :: cradius
      ! Vector of support radii of localization function
      REAL(c_double), DIMENSION(:), INTENT(in) :: sradius
      ! Coordinates of state vector elements
      REAL(c_double), DIMENSION(:,:), INTENT(in) :: coords
      ! Vector HP, dimension (dim)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hp
      ! Matrix HXY, dimension (nobs)
      REAL(c_double), DIMENSION(:), INTENT(inout) :: hxy


      call PDAFomi_localize_covar_serial_noniso(thisobs(i_obs), iobs_all, dim,  &
         dim_obs, locweight, cradius, sradius, coords, hp, hxy)

   END SUBROUTINE c__PDAFomi_localize_covar_serial_noniso

   SUBROUTINE c__PDAFomi_init_dim_obs_l_iso_old(i_obs, coords_l, locweight, cradius,  &
      sradius, cnt_obs_l) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Coordinates of current analysis domain
      REAL(c_double), DIMENSION(:), INTENT(in) :: coords_l
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight
      ! Localization cut-off radius (single or vector)
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius of localization function (single or vector)
      REAL(c_double), INTENT(in) :: sradius
      ! Local dimension of current observation vector
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l


      call PDAFomi_init_dim_obs_l_iso_old(thisobs_l(i_obs), thisobs(i_obs),  &
         coords_l, locweight, cradius, sradius, cnt_obs_l)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_iso_old

   SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_old(i_obs, coords_l, locweight,  &
      cradius, sradius, cnt_obs_l) bind(c)
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
      INTEGER(c_int), INTENT(inout) :: cnt_obs_l


      call PDAFomi_init_dim_obs_l_noniso_old(thisobs_l(i_obs), thisobs(i_obs),  &
         coords_l, locweight, cradius, sradius, cnt_obs_l)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_old

   SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_locweights_old(i_obs, coords_l,  &
      locweights, cradius, sradius, cnt_obs_l) bind(c)
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


      call PDAFomi_init_dim_obs_l_noniso_locweights_old(thisobs_l(i_obs),  &
         thisobs(i_obs), coords_l, locweights, cradius, sradius, cnt_obs_l)

   END SUBROUTINE c__PDAFomi_init_dim_obs_l_noniso_locweights_old

   SUBROUTINE c__PDAFomi_deallocate_obs(i_obs) bind(c)
      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs



      call PDAFomi_deallocate_obs(thisobs(i_obs))

   END SUBROUTINE c__PDAFomi_deallocate_obs
end module pdafomi_c_legacy
