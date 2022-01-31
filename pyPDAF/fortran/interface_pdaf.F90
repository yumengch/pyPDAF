module interface_pdaf
use pdafomi, only: obs_f, obs_l
implicit none

abstract interface
   subroutine c__init_ens_pdaf(filtertype, dim_p, dim_ens, &
         state_p, uinv, ens_p, flag) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! *** arguments ***
      !< type of filter to initialize
      integer(c_int), intent(in) :: filtertype
      !< pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      !< size of ensemble
      integer(c_int), intent(in) :: dim_ens
      !< pe-local model state
      real(c_double), intent(inout) :: state_p(dim_p)
      !< array not referenced for ensemble filters
      real(c_double), intent(inout) :: uinv(dim_ens - 1,dim_ens - 1)
      !< pe-local state ensemble
      real(c_double), intent(inout) :: ens_p(dim_p, dim_ens)
      !< pdaf status flag   
      integer(c_int), intent(inout) :: flag
   end subroutine c__init_ens_pdaf

   subroutine c__collect_state_pdaf(dim_p, state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      !< pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      !< local state vector
      real(c_double), intent(inout) :: state_p(dim_p)
   end subroutine c__collect_state_pdaf

   subroutine c__next_observation_pdaf(stepnow, nsteps, &
      doexit, time) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! *** arguments ***
      !< number of the current time step
      integer(c_int), intent(in)  :: stepnow
      !< number of time steps until next obs
      integer(c_int), intent(out) :: nsteps
      !< whether to exit forecasting (1 for exit)
      integer(c_int), intent(out) :: doexit
      !< current model (physical) time
      real(c_double), intent(out) :: time
   end subroutine c__next_observation_pdaf

   subroutine c__distribute_state_pdaf(dim_p, state_p) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      integer(c_int), intent(in) :: dim_p
      real(c_double), intent(inout) :: state_p(dim_p)
   end subroutine c__distribute_state_pdaf

   subroutine c__prepoststep_ens_pdaf(step, dim_p, dim_ens, &
      dim_ens_p, dim_obs_p, &
      state_p, uinv, ens_p, flag) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! *** arguments ***
      !< current time step (negative for call after forecast)
      integer(c_int), intent(in) :: step
      !< pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      !< size of state ensemble
      integer(c_int), intent(in) :: dim_ens
      !< pe-local size of ensemble
      integer(c_int), intent(in) :: dim_ens_p
      !< pe-local dimension of observation vector
      integer(c_int), intent(in) :: dim_obs_p
      !< pe-local forecast/analysis state
      !< (the array 'state_p' is not generally not
      !< initialized in the case of seik.
      !< it can be used freely here.)
      real(c_double), intent(inout) :: state_p(dim_p) 
      !< inverse of matrix u
      real(c_double), intent(inout) :: uinv(dim_ens-1, dim_ens-1)
      !< pe-local state ensemble
      real(c_double), intent(inout) :: ens_p(dim_p, dim_ens)
      !< pdaf status flag
      integer(c_int), intent(in) :: flag
   end subroutine c__prepoststep_ens_pdaf

   subroutine c__init_dim_obs_pdafomi(step, dim_obs) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      ! *** arguments ***
      !< current time step
      integer(c_int), intent(in)    :: step
      !< dimension of full observation vector
      integer(c_int), intent(inout) :: dim_obs
   end subroutine c__init_dim_obs_pdafomi

   subroutine c__obs_op_pdafomi(step, dim_p, dim_obs, &
                                state_p, ostate) bind(c)
      use iso_c_binding, only: c_double, c_int
      implicit none
      !< current time step
      integer(c_int), intent(in) :: step
      !< pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      !< dimension of full observed state
      integer(c_int), intent(in) :: dim_obs
      !< pe-local model state
      real(c_double), intent(in)    :: state_p(dim_p)
      !< pe-local full observed state
      real(c_double), intent(inout) :: ostate(dim_obs)
   end subroutine c__obs_op_pdafomi

   subroutine c__init_dim_obs_l_pdafomi(domain_p, step, &
                                        dim_obs, dim_obs_l) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      !< index of current local analysis domain
      integer(c_int), intent(in)  :: domain_p
      !< current time step
      integer(c_int), intent(in)  :: step
      !< full dimension of observation vector
      integer(c_int), intent(in)  :: dim_obs
      !< local dimension of observation vector
      integer(c_int), intent(out) :: dim_obs_l
   end subroutine c__init_dim_obs_l_pdafomi

   subroutine c__localize_covar_pdafomi(dim_p, dim_obs, &
                                       hp_p, hph) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none
      !< pe-local state dimension
      integer(c_int), intent(in) :: dim_p
      !< number of observations
      integer(c_int), intent(in) :: dim_obs
      !< pe local part of matrix hp
      real(c_double), intent(inout) :: hp_p(dim_obs, dim_p)
      !< matrix hph
      real(c_double), intent(inout) :: hph(dim_obs, dim_obs)
   end subroutine c__localize_covar_pdafomi

   subroutine c__init_n_domains_pdaf(step, n_domains_p) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      !< current time step
      integer(c_int), intent(in)  :: step
      !< pe-local number of analysis domains        
      integer(c_int), intent(out) :: n_domains_p
   end subroutine c__init_n_domains_pdaf

   subroutine c__init_dim_l_pdaf(step, domain_p, dim_l) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      ! *** arguments ***
      !< current time step
      integer(c_int), intent(in)  :: step     
      !< current local analysis domain
      integer(c_int), intent(in)  :: domain_p 
      !< local state dimension 
      integer(c_int), intent(out) :: dim_l
   end subroutine c__init_dim_l_pdaf

   subroutine c__g2l_state_pdaf(step, domain_p, dim_p, &
                                state_p, dim_l, state_l) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none
      ! *** arguments ***
      !< current time step
      integer(c_int), intent(in) :: step           
      !< current local analysis domain
      integer(c_int), intent(in) :: domain_p       
      !< pe-local full state dimension
      integer(c_int), intent(in) :: dim_p          
      !< local state dimension
      integer(c_int), intent(in) :: dim_l          
      !< pe-local full state vector 
      real(c_double), intent(in)    :: state_p(dim_p) 
      !< state vector on local analysis domain
      real(c_double), intent(out)   :: state_l(dim_l) 
   end subroutine c__g2l_state_pdaf

   subroutine c__l2g_state_pdaf(step, domain_p, dim_l, &
      state_l, dim_p, state_p) bind(c)
      use iso_c_binding, only: c_int, c_double
      implicit none

      ! *** arguments ***
      !< current time step
      integer(c_int), intent(in) :: step
      !< current local analysis domain
      integer(c_int), intent(in) :: domain_p
      !< local state dimension
      integer(c_int), intent(in) :: dim_l
      !< pe-local full state dimension
      integer(c_int), intent(in) :: dim_p
      !< state vector on local analysis domain
      real(c_double), intent(in)    :: state_l(dim_l)
      !< pe-local full state vector
      real(c_double), intent(inout) :: state_p(dim_p)
   end subroutine c__l2g_state_pdaf
end interface

type(obs_f), allocatable, target :: thisobs(:)
type(obs_l), allocatable, target :: thisobs_l(:)

contains
   subroutine c__pdaf_init(filtertype, subtype, &
                           filter_param_i, dim_pint, & 
                           filter_param_r, dim_preal, &
                           comm_model, comm_filter, comm_couple, &
                           task_id, n_modeltasks, &
                           filterpe, init_ens_pdaf, &
                           screen, status_pdaf) bind(c)
      use pdaf_interfaces_module, only: pdaf_init
      use iso_c_binding, only: c_int, c_double, c_bool
      implicit none
      integer(c_int), intent(in) :: filtertype
      integer(c_int), intent(in) :: subtype
      integer(c_int), intent(inout) :: filter_param_i(dim_pint)
      integer(c_int), intent(in) :: dim_pint
      real(c_double), intent(inout) :: filter_param_r(dim_preal)
      integer(c_int), intent(in) :: dim_preal
      integer(c_int), intent(in) :: comm_model
      integer(c_int), intent(in) :: comm_couple
      integer(c_int), intent(in) :: comm_filter
      integer(c_int), intent(in) :: task_id
      integer(c_int), intent(in) :: n_modeltasks
      logical(c_bool), intent(in) :: filterpe
      integer(c_int), intent(in) :: screen
      procedure(c__init_ens_pdaf) :: init_ens_pdaf
      integer(c_int), intent(out):: status_pdaf

      logical :: filterpe_local 
      filterpe_local = filterpe

      call pdaf_init(filtertype, subtype, 0, &
            filter_param_i, dim_pint,&
            filter_param_r, dim_preal, &
            comm_model, comm_filter, comm_couple, &
            task_id, n_modeltasks, filterpe_local, init_ens_pdaf, &
            screen, status_pdaf)
   end subroutine c__pdaf_init

   subroutine c__pdaf_get_state(steps, timenow, doexit, &
                                next_observation_pdaf, &
                                distribute_state_pdaf, &
                                prepoststep_ens_pdaf, &
                                status_pdaf) bind(c)
      use pdaf_interfaces_module, only: pdaf_get_state
      use iso_c_binding, only: c_int, c_double
      implicit none
      ! flag and number of time steps
      integer(c_int), intent(inout) :: steps      
      ! current model time
      real(c_double), intent(out)   :: timenow
      ! whether to exit from forecasts
      integer(c_int), intent(inout) :: doexit
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      procedure(c__prepoststep_ens_pdaf) :: prepoststep_ens_pdaf
      ! status flag
      integer(c_int), intent(inout) :: status_pdaf

      call pdaf_get_state(steps, timenow, doexit, &
            next_observation_pdaf, &
            distribute_state_pdaf, &
            prepoststep_ens_pdaf, status_pdaf)
   end subroutine c__pdaf_get_state

   subroutine c__set_pdafomi_doassim(i_obs, doassim) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in) :: doassim
      thisobs(i_obs)%doassim = doassim
   end subroutine c__set_pdafomi_doassim

   subroutine c__set_pdafomi_disttype(i_obs, disttype) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: disttype
      thisobs(i_obs)%disttype = disttype
   end subroutine c__set_pdafomi_disttype

   subroutine c__set_pdafomi_ncoord(i_obs, ncoord) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: ncoord
      thisobs(i_obs)%ncoord = ncoord
   end subroutine c__set_pdafomi_ncoord

   subroutine c__set_pdafomi_id_obs_p(i_obs, nrows, dim_obs_p, id_obs_p) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: nrows, dim_obs_p
      integer(c_int),intent(in)  :: id_obs_p(nrows, dim_obs_p)
      allocate(thisobs(i_obs)%id_obs_p(nrows, dim_obs_p))
      thisobs(i_obs)%id_obs_p(:, :) = id_obs_p(:, :)
   end subroutine c__set_pdafomi_id_obs_p

   subroutine c__set_pdafomi_icoeff_p(i_obs, nrows, dim_obs_p, icoeff_p) bind(c)
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: nrows, dim_obs_p
      real(c_double),intent(in)  :: icoeff_p(nrows, dim_obs_p)
      allocate(thisobs(i_obs)%icoeff_p(nrows, dim_obs_p))
      thisobs(i_obs)%icoeff_p = icoeff_p
   end subroutine c__set_pdafomi_icoeff_p

   subroutine c__set_pdafomi_domainsize(i_obs, ncoord, domainsize) bind(c)
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in) :: i_obs      
      integer(c_int), intent(in)  :: ncoord
      real(c_double),intent(in)  :: domainsize(ncoord)
      allocate(thisobs(i_obs)%domainsize(ncoord))
      thisobs(i_obs)%domainsize(:) = domainsize(:)
   end subroutine c__set_pdafomi_domainsize

   subroutine c__set_pdafomi_obs_err_type(i_obs, obs_err_type) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: obs_err_type
      thisobs(i_obs)%obs_err_type = obs_err_type
   end subroutine c__set_pdafomi_obs_err_type

   subroutine c__set_pdafomi_use_global_obs(i_obs, use_global_obs) bind(c)
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in)  :: use_global_obs
      thisobs(i_obs)%use_global_obs = use_global_obs
   end subroutine c__set_pdafomi_use_global_obs

   subroutine c__pdafomi_gather_obs(i_obs, nrows, dim_obs_p, &
                                    obs_p, ivar_obs_p, ocoord_p, &
                                    local_range, dim_obs) bind(c)
      use pdafomi, only: pdafomi_gather_obs
      use iso_c_binding, only: c_int, c_double
      implicit none
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in) :: nrows, dim_obs_p
      ! pe-local observation vector
      real(c_double), intent(in) :: obs_p(dim_obs_p)
      ! pe-local inverse observation error variance
      real(c_double), intent(in) :: ivar_obs_p(dim_obs_p)
      ! pe-local observation coordinates
      real(c_double), intent(in) :: ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p)
      real(c_double), intent(in) :: local_range
      ! Full number of observations
      integer(c_int), intent(out) :: dim_obs

      call pdafomi_gather_obs(thisobs(i_obs), dim_obs_p, &
         obs_p, ivar_obs_p, ocoord_p, &
         thisobs(i_obs)%ncoord, local_range, dim_obs)
   end subroutine c__pdafomi_gather_obs

   subroutine c__pdafomi_obs_op_gridpoint(i_obs, dim_p, dim_obs, & 
                                          state_p, ostate) bind(c)
      use pdafomi, only: pdafomi_obs_op_gridpoint
      use iso_c_binding, only: c_int, c_double
      implicit none
      !< pe-local state dimension
      integer(c_int), intent(in) :: i_obs, dim_p
      !< dimension of full observed state (all observed fields)
      integer(c_int), intent(in) :: dim_obs
      !< pe-local model state
      real(c_double), intent(in) :: state_p(dim_p)
      !< full observed state
      real(c_double), intent(inout) :: ostate(dim_obs)
      ! observation operator for observed grid point values
      call pdafomi_obs_op_gridpoint(thisobs(i_obs), state_p, ostate)
   end subroutine c__pdafomi_obs_op_gridpoint

   subroutine c__pdafomi_set_domain_limits(lim_coords) bind(c)
      use pdafomi, only: pdafomi_set_domain_limits
      use iso_c_binding, only: c_double
      real(c_double), intent(in) :: lim_coords(2,2)
      call pdafomi_set_domain_limits(lim_coords)
   end subroutine c__pdafomi_set_domain_limits

   subroutine c__PDAFomi_init_dim_obs_l(i_obs, coords_l, locweight, &
                                        local_range, srange, &
                                        dim_obs_l) bind(c)
      use pdafomi, only: pdafomi_init_dim_obs_l
      use iso_c_binding, only: c_int, c_double
      implicit none
      INTEGER(c_int), intent(in) :: i_obs
      ! Coordinates of current local analysis domain
      REAL(c_double), INTENT(in) :: coords_l(:)
      ! Type of localization function
      INTEGER(c_int), INTENT(in) :: locweight
      ! Localization radius
      REAL(c_double), INTENT(in) :: local_range
      ! Support radius of localization function
      REAL(c_double), INTENT(in) :: srange              
      ! Local dimension of current observation vector
      INTEGER(c_int), INTENT(inout) :: dim_obs_l      

      call pdafomi_init_dim_obs_l(thisobs_l(i_obs), thisobs(i_obs), &
            coords_l, &
            locweight, local_range, srange, dim_obs_l)
   end subroutine c__PDAFomi_init_dim_obs_l

   subroutine c__PDAFomi_localize_covar(i_obs, dim_p, dim_obs, &
                                        dim_coords, &
                                        locweight, local_range, & 
                                        srange, coords_p, hp_p, &
                                        hph) bind(c)
      use pdafomi, only: pdafomi_localize_covar
      use iso_c_binding, only: c_int, c_double
      implicit none
      INTEGER(c_int), INTENT(in) :: i_obs
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim_p, dim_obs, dim_coords
      ! Localization weight type
      INTEGER(c_int), INTENT(in) :: locweight
      ! localization radius
      REAL(c_double), INTENT(in)    :: local_range
      ! support radius for weight functions
      REAL(c_double), INTENT(in)    :: srange
      ! Coordinates of state vector elements
      REAL(c_double), INTENT(in)    :: coords_p(dim_coords, dim_p)
      ! Matrix HP, dimension (nobs, dim)
      REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
      ! Matrix HPH, dimension (nobs, nobs)
      REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
      call pdafomi_localize_covar(thisobs(i_obs), dim_p, &
            locweight, local_range, srange, coords_p, hp_p, hph)
   end subroutine c__PDAFomi_localize_covar

   subroutine c__PDAFomi_deallocate_obs(i_obs, step) bind(c)
      ! include pdafomi function
      use pdafomi, only: pdafomi_deallocate_obs
      use iso_c_binding, only: c_int
      implicit none
      ! *** arguments ***
      integer(c_int), intent(in) :: i_obs
      integer(c_int), intent(in) :: step   !< current time step 

      call pdafomi_deallocate_obs(thisobs(i_obs))
   end subroutine c__PDAFomi_deallocate_obs

   subroutine c__pdaf_get_localfilter(localfilter) bind(c)
      use iso_c_binding, only: c_int
      use pdaf_interfaces_module, only: pdaf_get_localfilter
      implicit none
      integer(c_int), intent(out) :: localfilter
      call pdaf_get_localfilter(localfilter)
   end subroutine c__pdaf_get_localfilter

   subroutine c__pdafomi_assimilate_local(collect_state_pdaf, &
                                          distribute_state_pdaf, &
                                          init_dim_obs_pdafomi, &
                                          obs_op_pdafomi, &
                                          prepoststep_ens_pdaf, &
                                          init_n_domains_pdaf, &
                                          init_dim_l_pdaf, &
                                          init_dim_obs_l_pdafomi, &
                                          g2l_state_pdaf, &
                                          l2g_state_pdaf, &
                                          next_observation_pdaf, &
                                          status_pdaf) bind(c)
      use iso_c_binding, only: c_int
      use pdaf_interfaces_module, only: pdafomi_assimilate_local
      implicit none
      integer(c_int), intent(out) :: status_pdaf
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      procedure(c__init_dim_obs_pdafomi) :: init_dim_obs_pdafomi
      procedure(c__obs_op_pdafomi) :: obs_op_pdafomi
      procedure(c__prepoststep_ens_pdaf) :: prepoststep_ens_pdaf
      procedure(c__init_n_domains_pdaf) :: init_n_domains_pdaf
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      procedure(c__init_dim_obs_l_pdafomi) :: init_dim_obs_l_pdafomi
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      call pdafomi_assimilate_local(collect_state_pdaf, &
                                    distribute_state_pdaf, &
                                    init_dim_obs_pdafomi, &
                                    obs_op_pdafomi, &
                                    prepoststep_ens_pdaf, &
                                    init_n_domains_pdaf, &
                                    init_dim_l_pdaf, &
                                    init_dim_obs_l_pdafomi, &
                                    g2l_state_pdaf, l2g_state_pdaf, &
                                    next_observation_pdaf, &
                                    status_pdaf)
   end subroutine c__pdafomi_assimilate_local

   subroutine c__pdafomi_assimilate_global(collect_state_pdaf, &
                                          distribute_state_pdaf, &
                                          init_dim_obs_pdafomi, &
                                          obs_op_pdafomi, &
                                          prepoststep_ens_pdaf, &
                                          next_observation_pdaf, &
                                          status_pdaf) bind(c)
      use iso_c_binding, only: c_int
      use pdaf_interfaces_module, only: pdafomi_assimilate_global
      implicit none
      integer(c_int), intent(out) :: status_pdaf
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      procedure(c__init_dim_obs_pdafomi) :: init_dim_obs_pdafomi
      procedure(c__obs_op_pdafomi) :: obs_op_pdafomi
      procedure(c__prepoststep_ens_pdaf) :: prepoststep_ens_pdaf
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      call pdafomi_assimilate_global(collect_state_pdaf, &
                                     distribute_state_pdaf, &
                                     init_dim_obs_pdafomi, &
                                     obs_op_pdafomi, &
                                     prepoststep_ens_pdaf, &
                                     next_observation_pdaf, status_pdaf)
   end subroutine c__pdafomi_assimilate_global

   subroutine c__pdafomi_assimilate_lenkf(collect_state_pdaf, &
                                          distribute_state_pdaf, &
                                          init_dim_obs_pdafomi, &
                                          obs_op_pdafomi, &
                                          prepoststep_ens_pdaf, &
                                          localize_covar_pdafomi, &
                                          next_observation_pdaf, &
                                          status_pdaf) bind(c)
      use pdaf_interfaces_module, only: pdafomi_assimilate_lenkf
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(out) :: status_pdaf
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      procedure(c__init_dim_obs_pdafomi) :: init_dim_obs_pdafomi
      procedure(c__obs_op_pdafomi) :: obs_op_pdafomi
      procedure(c__prepoststep_ens_pdaf) :: prepoststep_ens_pdaf
      procedure(c__localize_covar_pdafomi) :: localize_covar_pdafomi
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! lenkf has its own omi interface routine
      call pdafomi_assimilate_lenkf(collect_state_pdaf, &
                                    distribute_state_pdaf, &
                                    init_dim_obs_pdafomi, &
                                    obs_op_pdafomi, &
                                    prepoststep_ens_pdaf, &
                                    localize_covar_pdafomi, &
                                    next_observation_pdaf, status_pdaf)
   end subroutine c__pdafomi_assimilate_lenkf

   subroutine c__init_pdafomi(n_obs) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(in) :: n_obs
      if (.not. allocated(thisobs)) allocate(thisobs(n_obs))
      if (.not. allocated(thisobs_l)) allocate(thisobs_l(n_obs))
   end subroutine c__init_pdafomi

end module interface_pdaf