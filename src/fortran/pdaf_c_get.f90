MODULE pdaf_c_get
use iso_c_binding, only: c_int, c_double, c_bool
use PDAF
implicit none

contains
   SUBROUTINE c__PDAF_get_assim_flag(did_assim) bind(c)
      ! Flag: (1) for assimilation; (0) else
      INTEGER(c_int), INTENT(out) :: did_assim


      call PDAF_get_assim_flag(did_assim)

   END SUBROUTINE c__PDAF_get_assim_flag

   SUBROUTINE c__PDAF_get_localfilter(localfilter_out) bind(c)
      ! Whether the filter is domain-localized
      INTEGER(c_int), INTENT(out) :: localfilter_out


      call PDAF_get_localfilter(localfilter_out)

   END SUBROUTINE c__PDAF_get_localfilter

   SUBROUTINE c__PDAF_get_local_type(localtype) bind(c)
      ! Localization type of the filter
      INTEGER(c_int), INTENT(out) :: localtype


      call PDAF_get_local_type(localtype)

   END SUBROUTINE c__PDAF_get_local_type

   SUBROUTINE c__PDAF_get_memberid(memberid) bind(c)
      ! Index in the local ensemble
      INTEGER(c_int), INTENT(inout) :: memberid


      call PDAF_get_memberid(memberid)

   END SUBROUTINE c__PDAF_get_memberid

   SUBROUTINE c__PDAF_get_obsmemberid(memberid) bind(c)
      ! Index in the local ensemble
      INTEGER(c_int), INTENT(inout) :: memberid


      call PDAF_get_obsmemberid(memberid)

   END SUBROUTINE c__PDAF_get_obsmemberid

   SUBROUTINE c__PDAF_get_seed(seedvec) bind(c)
      ! Random seed vector
      INTEGER(c_int), DIMENSION(4), INTENT(out) :: seedvec

      call PDAF_get_seed(seedvec)

   END SUBROUTINE c__PDAF_get_seed

   SUBROUTINE c__PDAF_get_seedvec(seedvec) bind(c)
      ! Random seed vector
      INTEGER(c_int), DIMENSION(4), INTENT(out) :: seedvec

      call PDAF_get_seedvec(seedvec)

   END SUBROUTINE c__PDAF_get_seedvec

   SUBROUTINE c__PDAF_get_rndcount(rndcount) bind(c)
      ! Number of random-number generation calls
      INTEGER(c_int), INTENT(out) :: rndcount

      call PDAF_get_rndcount(rndcount)

   END SUBROUTINE c__PDAF_get_rndcount

   SUBROUTINE c__PDAF_reset_fcst_flag(reset_fcst_flag) bind(c)
      ! Whether model time should be reset in the forecast
      INTEGER(c_int), INTENT(out) :: reset_fcst_flag

      reset_fcst_flag = PDAF_reset_fcst_flag()

   END SUBROUTINE c__PDAF_reset_fcst_flag

   SUBROUTINE c__PDAF_get_smootherens(sens_point, maxlag, status) bind(c)
      ! Pointer to smoother array
      REAL(c_double), POINTER, DIMENSION(:,:,:), INTENT(out) :: sens_point
      ! Number of past timesteps processed in sens
      INTEGER(c_int), INTENT(out) :: maxlag
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_get_smootherens(sens_point, maxlag, status)
   END SUBROUTINE c__PDAF_get_smootherens
END MODULE pdaf_c_get
