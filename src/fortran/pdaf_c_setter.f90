MODULE pdaf_c_setter
use iso_c_binding, only: c_int, c_double, c_bool
use PDAF
use pdaf_c_cb_interface
implicit none

contains
   SUBROUTINE c__PDAF_set_comm_pdaf(in_comm_pdaf) bind(c)
      ! MPI communicator for PDAF
      INTEGER(c_int), INTENT(in) :: in_comm_pdaf


      call PDAF_set_comm_pdaf(in_comm_pdaf)

   END SUBROUTINE c__PDAF_set_comm_pdaf

   SUBROUTINE c__PDAF_set_debug_flag(debugval) bind(c)
      ! Value for debugging flag
      INTEGER(c_int), INTENT(in) :: debugval


      call PDAF_set_debug_flag(debugval)

   END SUBROUTINE c__PDAF_set_debug_flag

   SUBROUTINE c__PDAF_set_ens_pointer(ens_ptr, status) bind(c)
      ! Pointer to ensemble array
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: ens_ptr
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_set_ens_pointer(ens_ptr, status)
   END SUBROUTINE c__PDAF_set_ens_pointer

   SUBROUTINE c__PDAF_set_iparam(id, value, flag) bind(c)
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_set_iparam(id, value, flag)

   END SUBROUTINE c__PDAF_set_iparam

   SUBROUTINE c__PDAF_set_memberid(memberid) bind(c)
      ! Index in the local ensemble
      INTEGER(c_int), INTENT(inout) :: memberid


      call PDAF_set_memberid(memberid)

   END SUBROUTINE c__PDAF_set_memberid

   SUBROUTINE c__PDAF_set_offline_mode(screen) bind(c)
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen


      call PDAF_set_offline_mode(screen)

   END SUBROUTINE c__PDAF_set_offline_mode

   SUBROUTINE c__PDAF_set_rparam(id, value, flag) bind(c)
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_set_rparam(id, value, flag)

   END SUBROUTINE c__PDAF_set_rparam

   SUBROUTINE c__PDAF_set_seedset(seedset_in) bind(c)
      ! Seedset index (1-20)
      INTEGER(c_int), INTENT(in) :: seedset_in


      call PDAF_set_seedset(seedset_in)

   END SUBROUTINE c__PDAF_set_seedset

   SUBROUTINE c__PDAF_set_smootherens(sens_point, maxlag, status) bind(c)
      ! Pointer to smoother array
      REAL(c_double), POINTER, DIMENSION(:,:,:), INTENT(out) :: sens_point
      ! Number of past timesteps in sens
      INTEGER(c_int), INTENT(in) :: maxlag
      ! Status flag,
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_set_smootherens(sens_point, maxlag, status)
   END SUBROUTINE c__PDAF_set_smootherens
END MODULE pdaf_c_setter
