module pdaf_c
use iso_c_binding, only: c_int, c_double, c_bool, c_char, c_null_char
use pdaf
use pdaf_c_cb_interface
implicit none

contains
   subroutine c__pdaf_flush_fortran_stdout() bind(C)
      use iso_fortran_env, only: output_unit, error_unit
      flush(output_unit)
      flush(error_unit)
   end subroutine c__pdaf_flush_fortran_stdout

   SUBROUTINE c__PDAF_print_version() bind(c)
      use PDAF_info
      implicit none
      call PDAF_print_version()

   END SUBROUTINE c__PDAF_print_version

   SUBROUTINE c__PDAF_configinfo_filters(subtype, verbose) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Sub-type of filter
      INTEGER(c_int), INTENT(inout) :: subtype
      ! Control screen output
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_configinfo_filters(subtype, verbose)

   END SUBROUTINE c__PDAF_configinfo_filters

   SUBROUTINE c__PDAF_options_filters(type_filter) bind(c)
      use PDAF_utils_filters
      implicit none
      ! Type of filter
      INTEGER(c_int), INTENT(in) :: type_filter


      call PDAF_options_filters(type_filter)

   END SUBROUTINE c__PDAF_options_filters

   SUBROUTINE c__PDAF_get_fcst_info(steps, time, doexit) bind(c)
      ! Flag and number of time steps
      INTEGER(c_int), INTENT(inout) :: steps
      ! current model time
      REAL(c_double), INTENT(inout) :: time
      ! Whether to exit from forecasts
      INTEGER(c_int), INTENT(inout) :: doexit


      call PDAF_get_fcst_info(steps, time, doexit)

   END SUBROUTINE c__PDAF_get_fcst_info

   SUBROUTINE c__PDAF_correlation_function(ctype, length, distance, value) bind(c)
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

   SUBROUTINE c__PDAF_finalize() bind(c)
      call PDAF_finalize()

   END SUBROUTINE c__PDAF_finalize

   SUBROUTINE c__PDAF_abort(err) bind(c)
      INTEGER(c_int), INTENT(in) :: err

      call PDAF_abort(err)

   END SUBROUTINE c__PDAF_abort

   SUBROUTINE c__PDAF_eofcovar(dim, nstates, nfields, dim_fields, offsets,  &
      remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose,  &
      status) bind(c)
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

   SUBROUTINE c__PDAF_generate_rndvec(len, vec, stddev, dist, iseed) bind(c)
      ! Length of vector to process
      INTEGER(c_int), INTENT(in) :: len
      ! Values to be perturbed
      REAL(c_double), DIMENSION(len), INTENT(inout) :: vec
      ! Standard deviation of random perturbation
      REAL(c_double), INTENT(in) :: stddev
      ! Distribution type
      INTEGER(c_int), INTENT(in) :: dist
      ! Seed vector for LAPACK dlarnv
      INTEGER(c_int), DIMENSION(4), INTENT(inout) :: iseed

      call PDAF_generate_rndvec(len, vec, stddev, dist, iseed)

   END SUBROUTINE c__PDAF_generate_rndvec

   SUBROUTINE c__PDAF_parse_int(handle, intvalue) bind(c)
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: handle
      INTEGER(c_int), INTENT(inout) :: intvalue
      CHARACTER(len=32) :: clean_handle
      INTEGER :: i

      clean_handle = ""
      i = 1
      DO WHILE (i <= LEN(clean_handle))
         IF (handle(i) == c_null_char) EXIT
         clean_handle(i:i) = handle(i)
         i = i + 1
      END DO

      call PDAF_parse(clean_handle, intvalue)

   END SUBROUTINE c__PDAF_parse_int

   SUBROUTINE c__PDAF_parse_real(handle, realvalue) bind(c)
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: handle
      REAL(c_double), INTENT(inout) :: realvalue
      CHARACTER(len=32) :: clean_handle
      INTEGER :: i

      clean_handle = ""
      i = 1
      DO WHILE (i <= LEN(clean_handle))
         IF (handle(i) == c_null_char) EXIT
         clean_handle(i:i) = handle(i)
         i = i + 1
      END DO

      call PDAF_parse(clean_handle, realvalue)

   END SUBROUTINE c__PDAF_parse_real

   SUBROUTINE c__PDAF_parse_string(handle, charvalue) bind(c)
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: handle
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(inout) :: charvalue
      CHARACTER(len=32) :: clean_handle
      CHARACTER(len=100) :: clean_charvalue
      INTEGER :: i

      clean_handle = ""
      i = 1
      DO WHILE (i <= LEN(clean_handle))
         IF (handle(i) == c_null_char) EXIT
         clean_handle(i:i) = handle(i)
         i = i + 1
      END DO

      clean_charvalue = ""
      i = 1
      DO WHILE (i <= LEN(clean_charvalue))
         IF (charvalue(i) == c_null_char) EXIT
         clean_charvalue(i:i) = charvalue(i)
         i = i + 1
      END DO

      call PDAF_parse(clean_handle, clean_charvalue)

      DO i = 1, LEN(clean_charvalue)
         charvalue(i) = clean_charvalue(i:i)
      END DO
      charvalue(LEN(clean_charvalue) + 1) = c_null_char

   END SUBROUTINE c__PDAF_parse_string

   SUBROUTINE c__PDAF_parse_logical(handle, logvalue) bind(c)
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: handle
      LOGICAL(c_bool), INTENT(inout) :: logvalue
      CHARACTER(len=32) :: clean_handle
      LOGICAL :: f_logvalue
      INTEGER :: i

      clean_handle = ""
      i = 1
      DO WHILE (i <= LEN(clean_handle))
         IF (handle(i) == c_null_char) EXIT
         clean_handle(i:i) = handle(i)
         i = i + 1
      END DO

      f_logvalue = LOGICAL(logvalue)
      call PDAF_parse(clean_handle, f_logvalue)
      logvalue = f_logvalue

   END SUBROUTINE c__PDAF_parse_logical

   SUBROUTINE c__PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f) bind(c)
      ! PE-local observation dimension
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Full observation dimension
      INTEGER(c_int), INTENT(out) :: dim_obs_f


      call PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

   END SUBROUTINE c__PDAF_gather_dim_obs_f

   SUBROUTINE c__PDAF_gather_obs_f(obs_p, obs_f, status) bind(c)
      USE PDAF_mod_parallel, ONLY: dimobs_p, dimobs_f
      implicit none
      ! PE-local vector
      REAL(c_double), DIMENSION(dimobs_p), INTENT(in) :: obs_p
      ! Full gathered vector
      REAL(c_double), DIMENSION(dimobs_f), INTENT(out) :: obs_f
      ! Status flag:
      INTEGER(c_int), INTENT(out) :: status


      call PDAF_gather_obs_f(obs_p, obs_f, status)

   END SUBROUTINE c__PDAF_gather_obs_f

   SUBROUTINE c__PDAF_gather_obs_f2(coords_p, coords_f, nrows, status) bind(c)
      USE PDAF_mod_parallel, ONLY: dimobs_p, dimobs_f

      implicit none
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
      use pdaf_c_f_interface, only: init_ens_pdaf_c_ptr, &
                                    f__init_ens_pdaf
      implicit none
      ! Type of filter
      INTEGER(c_int), INTENT(in) :: filtertype
      ! Sub-type of filter
      INTEGER(c_int), INTENT(in) :: subtype
      ! Initial time step of assimilation
      INTEGER(c_int), INTENT(in) :: stepnull
      ! Number of integer parameters
      INTEGER(c_int), INTENT(in) :: dim_pint
      ! Integer parameter array
      INTEGER(c_int), DIMENSION(dim_pint), INTENT(inout) :: param_int
      ! Number of real parameter
      INTEGER(c_int), INTENT(in) :: dim_preal
      ! Real parameter array
      REAL(c_double), DIMENSION(dim_preal), INTENT(inout) :: param_real
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
      procedure(c__init_ens_pdaf) :: u_init_ens

      init_ens_pdaf_c_ptr => u_init_ens
      call PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint,  &
         param_real, dim_preal, comm_model, comm_filter, comm_couple, task_id,  &
         n_modeltasks, logical(in_filterpe), f__init_ens_pdaf, in_screen, outflag)

   END SUBROUTINE c__PDAF_init

   SUBROUTINE c__PDAF_init_forecast(u_next_observation, u_distribute_state,  &
      u_prepoststep, outflag) bind(c)
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

      call PDAF_init_forecast(f__next_observation_pdaf, f__distribute_state_pdaf,  &
         f__prepoststep_pdaf, outflag)

   END SUBROUTINE c__PDAF_init_forecast

   SUBROUTINE c__PDAF_local_weight(wtype, rtype, cradius, sradius, distance,  &
      nrows, ncols, a, var_obs, weight, verbose) bind(c)
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
      !
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_print_filter_types(verbose)

   END SUBROUTINE c__PDAF_print_filter_types

   SUBROUTINE c__PDAF_print_DA_types(verbose) bind(c)
      !
      INTEGER(c_int), INTENT(in) :: verbose


      call PDAF_print_DA_types(verbose)

   END SUBROUTINE c__PDAF_print_DA_types

   SUBROUTINE c__PDAF_print_info(printtype) bind(c)
      ! Type of screen output:
      INTEGER(c_int), INTENT(in) :: printtype


      call PDAF_print_info(printtype)

   END SUBROUTINE c__PDAF_print_info

   SUBROUTINE c__PDAF_reset_forget(forget_in) bind(c)
      ! New value of forgetting factor
      REAL(c_double), INTENT(in) :: forget_in


      call PDAF_reset_forget(forget_in)

   END SUBROUTINE c__PDAF_reset_forget

   SUBROUTINE c__PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens,  &
      verbose, flag) bind(c)
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
