module pdaf_c
use iso_c_binding, only: c_int, c_double, c_bool
use pdaf
implicit none

contains
   SUBROUTINE c__PDAF3_init(filtertype, subtype, stepnull, param_int, dim_pint, &
      param_real, dim_preal, U_init_ens, in_screen, outflag) bind(c)
      use pdaf_c_f_interface, only: init_ens_pdaf_c_ptr, &
                                    f__init_ens_pdaf
      IMPLICIT NONE

      ! *** Arguments ***
      ! For valid and default values see PDAF_mod_core.F90
      !< Type of filter
      INTEGER(c_int), INTENT(in) :: filtertype
      !< Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      !< Initial time step of assimilation
      INTEGER(c_int), INTENT(in) :: stepnull
      !< Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      !< Integer parameter array
      INTEGER(c_int), dimension(dim_pint), INTENT(inout) :: param_int
      !< Number of real parameter
      INTEGER(c_int), INTENT(in) :: dim_preal
      !< Real parameter array
      REAL(c_double), dimension(dim_preal), INTENT(inout) :: param_real
      !< Control screen output:
      !< (0) none, (1) some, default, (2) extensive
      INTEGER(c_int), INTENT(in) :: in_screen
      !< Status flag, 0: no error, error codes:
      !< -1: Call with subtype=-1 for info display
      !<  1: No valid filter type
      !<  2: No valid sub type
      !<  3: Invalid dim_pint
      !<  4: Invalid dim_preal
      !<  5: Invalid state dimension
      !<  6: Invalid ensemble size
      !<  7: Invalid value for forgetting factor
      !<  8: Invalid other integer parameter value
      !<  9: Invalid other real parameter value
      !< 10: MPI information not initialized
      !< 20: error in allocation of array at PDAF init
      INTEGER(c_int), INTENT(out):: outflag
      ! *** External subroutines ***
      ! (PDAF-internal names, real names are defined in the call to PDAF)
      ! User-supplied routine for ensemble initialization
      procedure(c__init_ens_pdaf) :: u_init_ens

      init_ens_pdaf_c_ptr => u_init_ens

      call PDAF3_init(filtertype, subtype, stepnull, param_int, dim_pint, &
      param_real, dim_preal, f__init_ens_pdaf, in_screen, outflag)

   END SUBROUTINE c__PDAF3_init

   SUBROUTINE c__PDAF3_init_forecast(U_next_observation, U_distribute_state, &
       U_prepoststep, outflag) bind(c)
      use pdaf_c_f_interface, only: next_observation_pdaf_c_ptr, &
                                    f__next_observation_pdaf, &
                                    distribute_state_pdaf_c_ptr, &
                                    f__distribute_state_pdaf, &
                                    prepoststep_pdaf_c_ptr, &
                                    f__prepoststep_pdaf
      implicit none
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Provide information on next forecast
      procedure(c__next_observation_pdaf) :: u_next_observation
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: u_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: u_prepoststep

      next_observation_pdaf_c_ptr => u_next_observation
      distribute_state_pdaf_c_ptr => u_distribute_state
      prepoststep_pdaf_c_ptr => u_prepoststep

      call PDAF3_init_forecast(f__next_observation_pdaf, f__distribute_state_pdaf,  &
         f__prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAF3_init_forecast

   SUBROUTINE c__PDAF3_set_parallel(in_COMM_pdaf, in_COMM_model, in_COMM_filter, in_COMM_couple, &
      in_task_id, in_n_modeltasks, in_filterpe, flag) bind(c)
      IMPLICIT NONE
      !< MPI communicator for all PEs involved in PDAF
      INTEGER(c_int), INTENT(in) :: in_COMM_pdaf
      !< Model communicator
      INTEGER(c_int), INTENT(in) :: in_COMM_model
      !< Filter communicator
      INTEGER(c_int), INTENT(in) :: in_COMM_filter
      !< Coupling communicator
      INTEGER(c_int), INTENT(in) :: in_COMM_couple
      !< Task ID of current PE
      INTEGER(c_int), INTENT(in) :: in_task_id
      !< Number of model tasks
      INTEGER(c_int), INTENT(in) :: in_n_modeltasks
      !< Is my PE a filter-PE?
      LOGICAL(c_bool), INTENT(in) :: in_filterpe
      !< Status flag
      INTEGER(c_int), INTENT(inout):: flag

      call PDAF3_set_parallel(in_COMM_pdaf, in_COMM_model, in_COMM_filter, in_COMM_couple, &
         in_task_id, in_n_modeltasks, logical(in_filterpe), flag)
   END SUBROUTINE c__PDAF3_set_parallel