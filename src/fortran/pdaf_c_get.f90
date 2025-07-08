MODULE pdaf_c_get
use PDAF
implicit none

contains
   SUBROUTINE c__PDAF_get_assim_flag(did_assim) bind(c)
      use iso_c_binding

      ! Flag: (1) for assimilation; (0) else
      INTEGER(c_int), INTENT(out) :: did_assim


      call PDAF_get_assim_flag(did_assim)

   END SUBROUTINE c__PDAF_get_assim_flag

   SUBROUTINE c__PDAF_get_localfilter(localfilter_out) bind(c)
      use iso_c_binding

      ! Whether the filter is domain-localized
      INTEGER(c_int), INTENT(out) :: localfilter_out


      call PDAF_get_localfilter(localfilter_out)

   END SUBROUTINE c__PDAF_get_localfilter

   SUBROUTINE c__PDAF_get_local_type(localtype) bind(c)
      use iso_c_binding

      ! Localization type of the filter
      INTEGER(c_int), INTENT(out) :: localtype


      call PDAF_get_local_type(localtype)

   END SUBROUTINE c__PDAF_get_local_type

   SUBROUTINE c__PDAF_get_memberid(memberid) bind(c)
      use iso_c_binding

      ! Index in the local ensemble
      INTEGER(c_int), INTENT(inout) :: memberid


      call PDAF_get_memberid(memberid)

   END SUBROUTINE c__PDAF_get_memberid

   SUBROUTINE c__PDAF_get_obsmemberid(memberid) bind(c)
      use iso_c_binding

      ! Index in the local ensemble
      INTEGER(c_int), INTENT(inout) :: memberid


      call PDAF_get_obsmemberid(memberid)

   END SUBROUTINE c__PDAF_get_obsmemberid

   SUBROUTINE c__PDAF_get_smootherens(sens_point, maxlag, status) bind(c)
      use iso_c_binding

      ! Pointer to smoother array
      REAL(c_double), POINTER, DIMENSION(:,:,:), INTENT(out) :: sens_point
      ! Number of past timesteps processed in sens
      INTEGER(c_int), INTENT(out) :: maxlag
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_get_smootherens(sens_point, maxlag, status)
   END SUBROUTINE c__PDAF_get_smootherens
END MODULE pdaf_c_get
