MODULE pdafomi_c_setter
use PDAF
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAFomi_set_doassim(i_obs, doassim) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Flag whether to assimilate this observation type
      INTEGER(c_int), INTENT(in) :: doassim


      call PDAFomi_set_doassim(thisobs(i_obs), doassim)

   END SUBROUTINE c__PDAFomi_set_doassim

   SUBROUTINE c__PDAFomi_set_disttype(i_obs, disttype) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Index of distance type
      INTEGER(c_int), INTENT(in) :: disttype


      call PDAFomi_set_disttype(thisobs(i_obs), disttype)

   END SUBROUTINE c__PDAFomi_set_disttype

   SUBROUTINE c__PDAFomi_set_ncoord(i_obs, ncoord) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Number of coordinates
      INTEGER(c_int), INTENT(in) :: ncoord


      call PDAFomi_set_ncoord(thisobs(i_obs), ncoord)

   END SUBROUTINE c__PDAFomi_set_ncoord

   SUBROUTINE c__PDAFomi_set_obs_err_type(i_obs, obs_err_type) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Type of observation error
      INTEGER(c_int), INTENT(in) :: obs_err_type


      call PDAFomi_set_obs_err_type(thisobs(i_obs), obs_err_type)

   END SUBROUTINE c__PDAFomi_set_obs_err_type

   SUBROUTINE c__PDAFomi_set_use_global_obs(i_obs, use_global_obs) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Set whether to use global full observations
      INTEGER(c_int), INTENT(in) :: use_global_obs


      call PDAFomi_set_use_global_obs(thisobs(i_obs), use_global_obs)

   END SUBROUTINE c__PDAFomi_set_use_global_obs

   SUBROUTINE c__PDAFomi_set_inno_omit(i_obs, inno_omit) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Set observation omission error variance level
      REAL(c_double), INTENT(in) :: inno_omit


      call PDAFomi_set_inno_omit(thisobs(i_obs), inno_omit)

   END SUBROUTINE c__PDAFomi_set_inno_omit

   SUBROUTINE c__PDAFomi_set_inno_omit_ivar(i_obs, inno_omit_ivar) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! Value of inverse variance to omit observation
      REAL(c_double), INTENT(in) :: inno_omit_ivar


      call PDAFomi_set_inno_omit_ivar(thisobs(i_obs), inno_omit_ivar)

   END SUBROUTINE c__PDAFomi_set_inno_omit_ivar

   SUBROUTINE c__PDAFomi_set_id_obs_p(i_obs, nrows, dim_obs_p, id_obs_p) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! number of rows required in observation operator
      INTEGER(c_int), INTENT(in) :: nrows
      ! number of process local observations
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Indices of process-local observed field in state vector
      INTEGER(c_int), DIMENSION(nrows, dim_obs_p), INTENT(in) :: id_obs_p


      call PDAFomi_set_id_obs_p(thisobs(i_obs), nrows, dim_obs_p, id_obs_p)

   END SUBROUTINE c__PDAFomi_set_id_obs_p

   SUBROUTINE c__PDAFomi_set_icoeff_p(i_obs, nrows, dim_obs_p, icoeff_p) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! number of rows required in observation operator
      INTEGER(c_int), INTENT(in) :: nrows
      ! number of process local observations
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Interpolation coeffs. for obs. operator
      REAL(c_double), DIMENSION(nrows, dim_obs_p), INTENT(in) :: icoeff_p


      call PDAFomi_set_icoeff_p(thisobs(i_obs), nrows, dim_obs_p, icoeff_p)

   END SUBROUTINE c__PDAFomi_set_icoeff_p

   SUBROUTINE c__PDAFomi_set_domainsize(i_obs, ncoord, domainsize) bind(c)
      use iso_c_binding

      ! index into observation arrays
      INTEGER(c_int), INTENT(in) :: i_obs

      ! number of coordinates considered for localizations
      INTEGER(c_int), INTENT(in) :: ncoord
      ! Size of domain for periodicity (<=0 for no periodicity)
      REAL(c_double), DIMENSION(ncoord), INTENT(in) :: domainsize


      call PDAFomi_set_domainsize(thisobs(i_obs), ncoord, domainsize)

   END SUBROUTINE c__PDAFomi_set_domainsize

   SUBROUTINE c__PDAFomi_set_globalobs(globalobs_in) bind(c)
      use iso_c_binding

      ! Input value of globalobs
      INTEGER(c_int), INTENT(in) :: globalobs_in


      call PDAFomi_set_globalobs(globalobs_in)

   END SUBROUTINE c__PDAFomi_set_globalobs
END MODULE pdafomi_c_setter
