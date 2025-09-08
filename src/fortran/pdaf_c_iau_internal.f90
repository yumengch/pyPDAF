MODULE pdaf_c_iau_internal
use iso_c_binding, only: c_double, c_int, c_bool
use PDAF
use PDAF_iau
use pdaf_c_cb_interface
implicit none

contains
   SUBROUTINE c__PDAF_iau_init_weights(type_iau, nsteps_iau) bind(c)
      ! Type of IAU, (0) no IAU
      INTEGER(c_int), INTENT(in) :: type_iau
      ! number of time steps in IAU
      INTEGER(c_int), INTENT(in) :: nsteps_iau


      call PDAF_iau_init_weights(type_iau, nsteps_iau)

   END SUBROUTINE c__PDAF_iau_init_weights

   SUBROUTINE c__PDAF_iau_update_inc(ens_ana, state_ana) bind(c)
      ! PE-local analysis ensemble
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: ens_ana
      ! PE-local state vector
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state_ana

      call PDAF_iau_update_inc(ens_ana, state_ana)

   END SUBROUTINE c__PDAF_iau_update_inc

   SUBROUTINE c__PDAF_iau_add_inc_ens(step, dim_p, dim_ens_task, ens,  &
      state, u_collect_state, u_distribute_state) bind(c)
      ! Time step
      INTEGER(c_int), INTENT(in) :: step
      ! PE-local dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! Ensemble size of model task
      INTEGER(c_int), INTENT(in) :: dim_ens_task
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: ens
      ! PE-local state vector
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state

      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: u_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state

      call PDAF_iau_add_inc_ens(step, dim_p, dim_ens_task, ens,  state, &
         u_collect_state, u_distribute_state)

   END SUBROUTINE c__PDAF_iau_add_inc_ens

   SUBROUTINE c__PDAF_iau_update_ens(ens, state) bind(c)
      ! PE-local state ensemble
      REAL(c_double), DIMENSION(:, :), INTENT(inout) :: ens
      ! PE-local state vector
      REAL(c_double), DIMENSION(:), INTENT(inout) :: state

      call PDAF_iau_update_ens(ens, state)

   END SUBROUTINE c__PDAF_iau_update_ens

   SUBROUTINE c__PDAF_iau_dealloc() bind(c)
      call PDAF_iau_dealloc()
   END SUBROUTINE c__PDAF_iau_dealloc
END MODULE pdaf_c_iau_internal
