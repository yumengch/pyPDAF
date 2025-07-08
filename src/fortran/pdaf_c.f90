module pdaf_c
implicit none

contains
   SUBROUTINE c__PDAF_correlation_function(ctype, length, distance, value) bind(c)
      use iso_c_binding

      ! Type of correlation function
      INTEGER(c_int), INTENT(in) :: ctype
      ! Length scale of function
      REAL(c_double), INTENT(in) :: length
      ! Distance at which the function is evaluated
      REAL(c_double), INTENT(in) :: distance
      ! Value of the function
      REAL(c_double), INTENT(out) :: value


      call PDAF_correlation_function(ctype, length, distance, value)

   END SUBROUTINE c__PDAF_correlation_function

   SUBROUTINE c__PDAF_deallocate() bind(c)
      call PDAF_deallocate()

   END SUBROUTINE c__PDAF_deallocate

   SUBROUTINE c__PDAF_eofcovar(dim, nstates, nfields, dim_fields, offsets,  &
      remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose,  &
      status) bind(c)
      use iso_c_binding

      ! Dimension of state vector
      INTEGER(c_int), INTENT(in) :: dim
      ! Number of state vectors
      INTEGER(c_int), INTENT(in) :: nstates
      ! Number of fields in state vector
      INTEGER(c_int), INTENT(in) :: nfields
      ! Size of each field
      INTEGER(c_int), DIMENSION(nfields), INTENT(in) :: dim_fields
      ! Start position of each field
      INTEGER(c_int), DIMENSION(nfields), INTENT(in) :: offsets
      ! 1: subtract mean state from states
      INTEGER(c_int), INTENT(in) :: remove_mstate
      ! 1: Do multivariate scaling; 0: no scaling
      INTEGER(c_int), INTENT(in) :: do_mv
      ! State perturbations
      REAL(c_double), DIMENSION(dim, nstates), INTENT(inout) :: states
      ! Standard deviation of field variability
      REAL(c_double), DIMENSION(nfields), INTENT(out) :: stddev
      ! Singular values divided by sqrt(nstates-1)
      REAL(c_double), DIMENSION(nstates), INTENT(out) :: svals
      ! Singular vectors
      REAL(c_double), DIMENSION(dim, nstates), INTENT(out) :: svec
      ! Mean state (only changed if remove_mstate=1)
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: meanstate
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_eofcovar(dim, nstates, nfields, dim_fields, offsets,  &
         remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose,  &
         status)

   END SUBROUTINE c__PDAF_eofcovar

   SUBROUTINE c__PDAF_force_analysis() bind(c)
      call PDAF_force_analysis()

   END SUBROUTINE c__PDAF_force_analysis

   SUBROUTINE c__PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f) bind(c)
      use iso_c_binding

      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_f


      call PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

   END SUBROUTINE c__PDAF_gather_dim_obs_f

   SUBROUTINE c__PDAF_gather_obs_f(obs_p, obs_f, status) bind(c)
      use iso_c_binding

      ! PE-local vector
      REAL(c_double), DIMENSION(dimobs_p), INTENT(in) :: obs_p
      ! Full gathered vector
      REAL(c_double), DIMENSION(dimobs_f), INTENT(out) :: obs_f
      ! Status flag:
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_gather_obs_f(obs_p, obs_f, status)

   END SUBROUTINE c__PDAF_gather_obs_f

   SUBROUTINE c__PDAF_gather_obs_f2(coords_p, coords_f, nrows, status) bind(c)
      use iso_c_binding

      ! PE-local array
      REAL(c_double), DIMENSION(nrows, dimobs_p), INTENT(in) :: coords_p
      ! Full gathered array
      REAL(c_double), DIMENSION(nrows, dimobs_f), INTENT(out) :: coords_f
      ! Number of rows in array
      INTEGER(c_int), INTENT(in) :: nrows
      ! Status flag:
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_gather_obs_f2(coords_p, coords_f, nrows, status)

   END SUBROUTINE c__PDAF_gather_obs_f2

   SUBROUTINE c__PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, obs_f,  &
      status) bind(c)
      use iso_c_binding

      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! PE-local vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Full gathered vector
      REAL(c_double), DIMENSION(dim_obs_f), INTENT(out) :: obs_f
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, obs_f, status)

   END SUBROUTINE c__PDAF_gather_obs_f_flex

   SUBROUTINE c__PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p,  &
      coords_f, nrows, status) bind(c)
      use iso_c_binding

      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_f
      ! PE-local array
      REAL(c_double), DIMENSION(nrows, dim_obs_p), INTENT(in) :: coords_p
      ! Full gathered array
      REAL(c_double), DIMENSION(nrows, dim_obs_f), INTENT(out) :: coords_f
      ! Number of rows in array
      INTEGER(c_int), INTENT(in) :: nrows
      ! Status flag: (0) no error
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p, coords_f,  &
         nrows, status)

   END SUBROUTINE c__PDAF_gather_obs_f2_flex

   SUBROUTINE c__PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint,  &
      param_real, dim_preal, comm_model, comm_filter, comm_couple, task_id,  &
      n_modeltasks, in_filterpe, u_init_ens, in_screen, outflag) bind(c)
      use iso_c_binding

      ! Type of filter
      INTEGER(c_int), INTENT(in) :: filtertype
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Initial time step of assimilation
      INTEGER(c_int), INTENT(in) :: stepnull
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
      ! Number of real parameter
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Model communicator
      INTEGER(c_int), INTENT(in) :: comm_model
      ! Filter communicator
      INTEGER(c_int), INTENT(in) :: comm_filter
      ! Coupling communicator
      INTEGER(c_int), INTENT(in) :: comm_couple
      ! Id of my ensemble task
      INTEGER(c_int), INTENT(in) :: task_id
      ! Number of parallel model tasks
      INTEGER(c_int), INTENT(in) :: n_modeltasks
      ! Is my PE a filter-PE?
      LOGICAL(c_bool), INTENT(in) :: in_filterpe
      ! Control screen output:
      INTEGER(c_int), INTENT(in) :: in_screen
      ! Status flag, 0: no error, error codes:
      INTEGER(c_int), INTENT(out) :: outflag

      ! User-supplied routine for ensemble initialization
      procedure(c__u_init_ens_pdaf) :: u_init_ens

      call PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint,  &
         param_real, dim_preal, comm_model, comm_filter, comm_couple, task_id,  &
         n_modeltasks, in_filterpe, u_init_ens, in_screen, outflag)

   END SUBROUTINE c__PDAF_init

   SUBROUTINE c__PDAF_init_forecast(u_next_observation, u_distribute_state,  &
      u_prepoststep, outflag) bind(c)
      use iso_c_binding

      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag

      ! Provide information on next forecast
      procedure(c__u_next_observation_pdaf) :: u_next_observation
      ! Routine to distribute a state vector
      procedure(c__u_distribute_state_pdaf) :: u_distribute_state
      ! User supplied pre/poststep routine
      procedure(c__u_prepoststep_pdaf) :: u_prepoststep

      call PDAF_init_forecast(u_next_observation, u_distribute_state,  &
         u_prepoststep, outflag)

   END SUBROUTINE c__PDAF_init_forecast

   SUBROUTINE c__PDAF_local_weight(wtype, rtype, cradius, sradius, distance,  &
      nrows, ncols, a, var_obs, weight, verbose) bind(c)
      use iso_c_binding

      ! Type of weight function
      INTEGER(c_int), INTENT(in) :: wtype
      ! Type of regulated weighting
      INTEGER(c_int), INTENT(in) :: rtype
      ! Cut-off radius
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius
      REAL(c_double), INTENT(in) :: sradius
      ! Distance to observation
      REAL(c_double), INTENT(in) :: distance
      ! Number of rows in matrix A
      INTEGER(c_int), INTENT(in) :: nrows
      ! Number of columns in matrix A
      INTEGER(c_int), INTENT(in) :: ncols
      ! Input matrix
      REAL(c_double), DIMENSION(nrows, ncols), INTENT(in) :: a
      ! Observation variance
      REAL(c_double), INTENT(in) :: var_obs
      ! Weights
      REAL(c_double), INTENT(out) :: weight
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_local_weight(wtype, rtype, cradius, sradius, distance, nrows,  &
         ncols, a, var_obs, weight, verbose)

   END SUBROUTINE c__PDAF_local_weight

   SUBROUTINE c__PDAF_local_weights(wtype, cradius, sradius, dim, distance,  &
      weight, verbose) bind(c)
      use iso_c_binding

      ! Type of weight function
      INTEGER(c_int), INTENT(in) :: wtype
      ! Parameter for cut-off
      REAL(c_double), INTENT(in) :: cradius
      ! Support radius
      REAL(c_double), INTENT(in) :: sradius
      ! Size of distance and weight arrays
      INTEGER(c_int), INTENT(in) :: dim
      ! Array holding distances
      REAL(c_double), DIMENSION(dim), INTENT(in) :: distance
      ! Array for weights
      REAL(c_double), DIMENSION(dim), INTENT(out) :: weight
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_local_weights(wtype, cradius, sradius, dim, distance, weight,  &
         verbose)

   END SUBROUTINE c__PDAF_local_weights

   SUBROUTINE c__PDAF_print_filter_types(verbose) bind(c)
      use iso_c_binding

      !
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_print_filter_types(verbose)

   END SUBROUTINE c__PDAF_print_filter_types

   SUBROUTINE c__PDAF_print_DA_types(verbose) bind(c)
      use iso_c_binding

      !
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_print_DA_types(verbose)

   END SUBROUTINE c__PDAF_print_DA_types

   SUBROUTINE c__PDAF_print_info(printtype) bind(c)
      use iso_c_binding

      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_print_info(printtype)

   END SUBROUTINE c__PDAF_print_info

   SUBROUTINE c__PDAF_reset_forget(forget_in) bind(c)
      use iso_c_binding

      ! New value of forgetting factor
      REAL(c_double), INTENT(in) :: forget_in


      call PDAF_reset_forget(forget_in)

   END SUBROUTINE c__PDAF_reset_forget

   SUBROUTINE c__PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens,  &
      verbose, flag) bind(c)
      use iso_c_binding

      ! Size of state vector
      INTEGER(c_int), INTENT(in) :: dim
      ! Size of ensemble
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Array of EOF modes
      REAL(c_double), DIMENSION(dim, dim_ens-1), INTENT(inout) :: modes
      ! Vector of singular values
      REAL(c_double), DIMENSION(dim_ens-1), INTENT(in) :: svals
      ! PE-local model state
      REAL(c_double), DIMENSION(dim), INTENT(inout) :: state
      ! State ensemble
      REAL(c_double), DIMENSION(dim, dim_ens), INTENT(out) :: ens
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: verbose
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag


      call PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens, verbose, flag)

   END SUBROUTINE c__PDAF_SampleEns
end module pdaf_c
