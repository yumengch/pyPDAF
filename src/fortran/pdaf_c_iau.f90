MODULE pdaf_c_iau
use PDAF
use U_PDAF_interface_c_binding

implicit none

contains
   SUBROUTINE c__PDAF_iau_init(type_iau_in, nsteps_iau_in, flag) bind(c)
      use iso_c_binding

      ! Type of IAU, (0) no IAU
      INTEGER(c_int), INTENT(in) :: type_iau_in
      ! number of time steps in IAU
      INTEGER(c_int), INTENT(in) :: nsteps_iau_in
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_init(type_iau_in, nsteps_iau_in, flag)

   END SUBROUTINE c__PDAF_iau_init

   SUBROUTINE c__PDAF_iau_reset(type_iau_in, nsteps_iau_in, flag) bind(c)
      use iso_c_binding

      ! Type of IAU, (0) no IAU
      INTEGER(c_int), INTENT(in) :: type_iau_in
      ! number of time steps in IAU
      INTEGER(c_int), INTENT(in) :: nsteps_iau_in
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_reset(type_iau_in, nsteps_iau_in, flag)

   END SUBROUTINE c__PDAF_iau_reset

   SUBROUTINE c__PDAF_iau_set_weights(iweights, weights) bind(c)
      use iso_c_binding

      ! Length of weights input vector
      INTEGER(c_int), INTENT(in) :: iweights
      ! Input weight vector
      REAL(c_double), DIMENSION(iweights), INTENT(in) :: weights


      call PDAF_iau_set_weights(iweights, weights)

   END SUBROUTINE c__PDAF_iau_set_weights

   SUBROUTINE c__PDAF_iau_set_pointer(iau_ptr, flag) bind(c)
      use iso_c_binding

      ! Pointer to IAU ensemble array
      REAL(c_double), POINTER, DIMENSION(:,:), INTENT(out) :: iau_ptr
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_iau_set_pointer(iau_ptr, flag)
   END SUBROUTINE c__PDAF_iau_set_pointer

   SUBROUTINE c__PDAF_iau_init_inc(dim_p, dim_ens_l, ens_inc, flag) bind(c)
      use iso_c_binding

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
      use iso_c_binding


      ! Routine to collect a state vector
      procedure(c__u_collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__u_distribute_state_pdaf) :: u_distribute_state

      call PDAF_iau_add_inc(u_collect_state, u_distribute_state)

   END SUBROUTINE c__PDAF_iau_add_inc

END MODULE pdaf_c_iau