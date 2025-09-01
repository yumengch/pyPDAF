MODULE pdaf_c_iau
use iso_c_binding, only: c_int, c_double, c_bool
use PDAF
use pdaf_c_cb_interface

implicit none

contains
   SUBROUTINE c__PDAF_iau_init(type_iau_in, nsteps_iau_in, flag) bind(c)
      ! Type of IAU, (0) no IAU
      INTEGER(c_int), INTENT(in) :: type_iau_in
      ! number of time steps in IAU
      INTEGER(c_int), INTENT(in) :: nsteps_iau_in
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_init(type_iau_in, nsteps_iau_in, flag)

   END SUBROUTINE c__PDAF_iau_init

   SUBROUTINE c__PDAF_iau_reset(type_iau_in, nsteps_iau_in, flag) bind(c)
      ! Type of IAU, (0) no IAU
      INTEGER(c_int), INTENT(in) :: type_iau_in
      ! number of time steps in IAU
      INTEGER(c_int), INTENT(in) :: nsteps_iau_in
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_reset(type_iau_in, nsteps_iau_in, flag)

   END SUBROUTINE c__PDAF_iau_reset

   SUBROUTINE c__PDAF_iau_set_weights(iweights, weights) bind(c)
      ! Length of weights input vector
      INTEGER(c_int), INTENT(in) :: iweights
      ! Input weight vector
      REAL(c_double), DIMENSION(iweights), INTENT(in) :: weights


      call PDAF_iau_set_weights(iweights, weights)

   END SUBROUTINE c__PDAF_iau_set_weights

   SUBROUTINE c__PDAF_iau_set_pointer(iau_ptr, flag) bind(c)
      ! Pointer to IAU ensemble array
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: iau_ptr
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_set_pointer(iau_ptr, flag)
   END SUBROUTINE c__PDAF_iau_set_pointer

   SUBROUTINE c__PDAF_iau_init_inc(dim_p, dim_ens_l, ens_inc, flag) bind(c)
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Task-local size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens_l
      ! PE-local increment ensemble
      REAL(c_double), DIMENSION(dim_p, dim_ens_l), INTENT(in) :: ens_inc
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_init_inc(dim_p, dim_ens_l, ens_inc, flag)

   END SUBROUTINE c__PDAF_iau_init_inc

   SUBROUTINE c__PDAF_iau_add_inc(u_collect_state, u_distribute_state) bind(c)
      use pdaf_c_f_interface, only: collect_state_pdaf_c_ptr, &
                                    distribute_state_pdaf_c_ptr, &
                                    f__collect_state_pdaf, f__distribute_state_pdaf
      implicit none
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state

      collect_state_pdaf_c_ptr => u_collect_state
      distribute_state_pdaf_c_ptr => u_distribute_state

      call PDAF_iau_add_inc(f__collect_state_pdaf, f__distribute_state_pdaf)

   END SUBROUTINE c__PDAF_iau_add_inc

   SUBROUTINE c__PDAF_iau_set_ens_pointer(iau_ptr, flag) bind(c)
      use PDAF_IAU, only: PDAF_iau_set_ens_pointer
      IMPLICIT NONE
      !< Pointer to IAU ensemble array
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: iau_ptr
      !< Status flag
      INTEGER(c_int), INTENT(out)       :: flag

      call PDAF_iau_set_ens_pointer(iau_ptr, flag)

   END SUBROUTINE c__PDAF_iau_set_ens_pointer

   SUBROUTINE c__PDAF_iau_set_state_pointer(iau_x_ptr, flag)
      use PDAF_IAU, only: PDAF_iau_set_state_pointer
      IMPLICIT NONE
      !< Pointer to IAU state vector
      REAL(c_double), POINTER, DIMENSION(:), INTENT(out) :: iau_x_ptr
      !< Status flag
      INTEGER(c_int), INTENT(out)       :: flag

      call PDAF_iau_set_state_pointer(iau_x_ptr, flag)

   END SUBROUTINE c__PDAF_iau_set_state_pointer
END MODULE pdaf_c_iau