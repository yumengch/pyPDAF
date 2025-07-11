MODULE pdaf_c_setter
use PDAF

implicit none

contains
   SUBROUTINE c__PDAF_set_iparam_filters(id, value, flag) bind(c)
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      INTEGER(c_int), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_set_iparam_filters(id, value, flag)

   END SUBROUTINE c__PDAF_set_iparam_filters

   SUBROUTINE c__PDAF_set_rparam_filters(id, value, flag) bind(c)
      ! Index of parameter
      INTEGER(c_int), INTENT(in) :: id
      ! Parameter value
      REAL(c_double), INTENT(in) :: value
      ! Status flag: 0 for no error
      INTEGER(c_int), INTENT(out) :: flag


      call PDAF_set_rparam_filters(id, value, flag)

   END SUBROUTINE c__PDAF_set_rparam_filters

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

   SUBROUTINE c__PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, mens_p,  &
      mstate_p, obs_p, u_init_obsvar, forget_in, forget_out, screen) bind(c)
      ! Current time step
      INTEGER(c_int), INTENT(in) :: step
      ! Whether filter is domain-local
      INTEGER(c_int), INTENT(in) :: localfilter
      ! Dimension of observation vector
      INTEGER(c_int), INTENT(in) :: dim_obs_p
      ! Ensemble size
      INTEGER(c_int), INTENT(in) :: dim_ens
      ! Observed PE-local ensemble
      REAL(c_double), DIMENSION(dim_obs_p, dim_ens), INTENT(in) :: mens_p
      ! Observed PE-local mean state
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: mstate_p
      ! Observation vector
      REAL(c_double), DIMENSION(dim_obs_p), INTENT(in) :: obs_p
      ! Prescribed forgetting factor
      REAL(c_double), INTENT(in) :: forget_in
      ! Adaptively estimated forgetting factor
      REAL(c_double), INTENT(out) :: forget_out
      ! Verbosity flag
      INTEGER(c_int), INTENT(in) :: screen

      ! Initialize mean obs. error variance
      procedure(c__init_obsvar_pdaf) :: u_init_obsvar

      call PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, mens_p,  &
         mstate_p, obs_p, u_init_obsvar, forget_in, forget_out, screen)

   END SUBROUTINE c__PDAF_set_forget



END MODULE pdaf_c_setter
