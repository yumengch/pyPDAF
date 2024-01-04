module PDAFomi_obs_c_binding
use pdafomi
use U_PDAF_interface_c_binding
use iso_c_binding, only: c_double, c_int, c_bool
implicit none

type(obs_f), allocatable, target :: thisobs(:)
type(obs_l), allocatable, target :: thisobs_l(:)
contains
   subroutine c__PDAFomi_init(n_obs) bind(c)
      ! number of observations
      integer(c_int), intent(in) :: n_obs
      if (.not. allocated(thisobs)) allocate(thisobs(n_obs))
      if (.not. allocated(thisobs_l)) allocate(thisobs_l(n_obs))
   end subroutine c__PDAFomi_init

   subroutine c__PDAFomi_set_doassim(i_obs, doassim) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! setter value
      integer(c_int), intent(in) :: doassim
      thisobs(i_obs)%doassim = doassim
   end subroutine c__PDAFomi_set_doassim

   subroutine c__PDAFomi_set_disttype(i_obs, disttype) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! setter value
      integer(c_int), intent(in)  :: disttype
      thisobs(i_obs)%disttype = disttype
   end subroutine c__PDAFomi_set_disttype

   subroutine c__PDAFomi_set_ncoord(i_obs, ncoord) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! setter value
      integer(c_int), intent(in)  :: ncoord
      thisobs(i_obs)%ncoord = ncoord
   end subroutine c__PDAFomi_set_ncoord

   subroutine c__PDAFomi_set_id_obs_p(i_obs, nrows, dim_obs_p, id_obs_p) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! Number of values to be averaged
      integer(c_int), intent(in)  :: nrows
      ! dimension of PE local obs
      integer(c_int), intent(in)  :: dim_obs_p
      ! setter value
      integer(c_int),intent(in)  :: id_obs_p(nrows, dim_obs_p)
      if (.not. allocated(thisobs(i_obs)%id_obs_p)) allocate(thisobs(i_obs)%id_obs_p(nrows, dim_obs_p))
      thisobs(i_obs)%id_obs_p(:, :) = id_obs_p(:, :)
   end subroutine c__PDAFomi_set_id_obs_p

   subroutine c__PDAFomi_set_icoeff_p(i_obs, nrows, dim_obs_p, icoeff_p) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! Number of values to be averaged
      integer(c_int), intent(in)  :: nrows
      ! dimension of PE local obs
      integer(c_int), intent(in)  :: dim_obs_p
      ! setter value
      real(c_double),intent(in)  :: icoeff_p(nrows, dim_obs_p)
      if (.not. allocated(thisobs(i_obs)%icoeff_p)) allocate(thisobs(i_obs)%icoeff_p(nrows, dim_obs_p))
      thisobs(i_obs)%icoeff_p = icoeff_p
   end subroutine c__PDAFomi_set_icoeff_p

   subroutine c__PDAFomi_set_domainsize(i_obs, ncoord, domainsize) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! state dimension
      integer(c_int), intent(in)  :: ncoord
      ! setter value
      real(c_double),intent(in)  :: domainsize(ncoord)
      if (.not. allocated(thisobs(i_obs)%domainsize)) allocate(thisobs(i_obs)%domainsize(ncoord))
      thisobs(i_obs)%domainsize(:) = domainsize(:)
   end subroutine c__PDAFomi_set_domainsize

   subroutine c__PDAFomi_set_obs_err_type(i_obs, obs_err_type) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! setter value
      integer(c_int), intent(in)  :: obs_err_type
      thisobs(i_obs)%obs_err_type = obs_err_type
   end subroutine c__PDAFomi_set_obs_err_type

   subroutine c__PDAFomi_set_use_global_obs(i_obs, use_global_obs) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! setter value
      integer(c_int), intent(in)  :: use_global_obs
      thisobs(i_obs)%use_global_obs = use_global_obs
   end subroutine c__PDAFomi_set_use_global_obs

   subroutine c__pdafomi_gather_obs(i_obs, dim_obs_p, &
                                    obs_p, ivar_obs_p, ocoord_p, &
                                    local_range, dim_obs) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs
      ! State dimension
      integer(c_int), intent(in) :: dim_obs_p
      ! pe-local observation vector
      real(c_double), intent(in) :: obs_p(dim_obs_p)
      ! pe-local inverse observation error variance
      real(c_double), intent(in) :: ivar_obs_p(dim_obs_p)
      ! pe-local observation coordinates
      real(c_double), intent(in) :: ocoord_p(thisobs(i_obs)%ncoord, dim_obs_p)
      ! localization radius
      real(c_double), intent(in) :: local_range
      ! Full number of observations
      integer(c_int), intent(out) :: dim_obs

      call pdafomi_gather_obs(thisobs(i_obs), dim_obs_p, &
         obs_p, ivar_obs_p, ocoord_p, &
         thisobs(i_obs)%ncoord, local_range, dim_obs)
   end subroutine c__pdafomi_gather_obs


   SUBROUTINE c__PDAFomi_gather_obsstate(i_obs, obsstate_p, obsstate_f, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), intent(in) :: i_obs
      ! Vector of process-local observed state
      REAL(c_double), INTENT(in) :: obsstate_p(thisobs(i_obs)%dim_obs_p)
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! Full observed vector for all types
      REAL(c_double), INTENT(inout) :: obsstate_f(nobs_f_all)

      call PDAFomi_gather_obsstate(thisobs(i_obs), obsstate_p, obsstate_f)
   END SUBROUTINE c__PDAFomi_gather_obsstate

   subroutine c__PDAFomi_localize_covar(i_obs, dim_p, dim_obs, &
                                        dim_coords, &
                                        locweight, local_range, &
                                        srange, coords_p, hp_p, &
                                        hph) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! State dimension
      INTEGER(c_int), INTENT(in) :: dim_p
      ! obs dimension
      INTEGER(c_int), INTENT(in) :: dim_obs
      ! number of coords
      INTEGER(c_int), INTENT(in) :: dim_coords
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

   SUBROUTINE c__PDAFomi_set_domain_limits(lim_coords) bind(c)
      ! geographic coordinate array (1: longitude, 2: latitude)
      REAL(c_double), INTENT(in) :: lim_coords(2,2)

      call PDAFomi_set_domain_limits(lim_coords)
   END SUBROUTINE c__PDAFomi_set_domain_limits

   SUBROUTINE c__PDAFomi_set_debug_flag(debugval)  bind(c)
      ! Value for debugging flag
      INTEGER(c_int), INTENT(in) :: debugval
      call PDAFomi_set_debug_flag(debugval)
   END SUBROUTINE c__PDAFomi_set_debug_flag

   subroutine c__PDAFomi_deallocate_obs(i_obs) bind(c)
      ! index of observations
      integer(c_int), intent(in) :: i_obs

      call pdafomi_deallocate_obs(thisobs(i_obs))
   end subroutine c__PDAFomi_deallocate_obs

   subroutine c__PDAFomi_init_dim_obs_l(i_obs, coords_l, locweight, &
                                        local_range, srange, &
                                        dim_obs_l) bind(c)
      ! index of observations
      INTEGER(c_int), intent(in) :: i_obs
      ! Coordinates of current local analysis domain
      REAL(c_double), INTENT(in) :: coords_l(2)
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

   SUBROUTINE c__PDAFomi_obs_op_gridpoint(i_obs, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(in)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(inout) :: obs_f_all(nobs_f_all)

      call PDAFomi_obs_op_gridpoint(thisobs(i_obs), state_p, obs_f_all)
   END SUBROUTINE c__PDAFomi_obs_op_gridpoint

   SUBROUTINE c__PDAFomi_obs_op_gridavg(i_obs, nrows, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(in)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(inout) :: obs_f_all(nobs_f_all)
      call PDAFomi_obs_op_gridavg(thisobs(i_obs), nrows, state_p, obs_f_all)
   END SUBROUTINE c__PDAFomi_obs_op_gridavg

   SUBROUTINE c__PDAFomi_obs_op_interp_lin(i_obs, nrows, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(in)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(inout) :: obs_f_all(nobs_f_all)
      call PDAFomi_obs_op_interp_lin(thisobs(i_obs), nrows, state_p, obs_f_all)
   END SUBROUTINE c__PDAFomi_obs_op_interp_lin

   SUBROUTINE c__PDAFomi_obs_op_adj_gridavg(i_obs, nrows, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(inout)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(in) :: obs_f_all(nobs_f_all)

      call PDAFomi_obs_op_adj_gridavg(thisobs(i_obs), nrows, obs_f_all, state_p)
   END SUBROUTINE c__PDAFomi_obs_op_adj_gridavg

   SUBROUTINE c__PDAFomi_obs_op_adj_gridpoint(i_obs, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(inout)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(in) :: obs_f_all(nobs_f_all)

      call PDAFomi_obs_op_adj_gridpoint(thisobs(i_obs), obs_f_all, state_p)
   END SUBROUTINE c__PDAFomi_obs_op_adj_gridpoint

   SUBROUTINE c__PDAFomi_obs_op_adj_interp_lin(i_obs, nrows, state_p, dim_p, obs_f_all, nobs_f_all) bind(c)
      ! index of observations
      INTEGER(c_int), INTENT(in) :: i_obs
      ! Number of values to be averaged
      INTEGER(c_int), INTENT(in) :: nrows
      ! dimension of model state
      INTEGER(c_int), INTENT(in) :: dim_p
      ! dimension of the observation
      INTEGER(c_int), INTENT(in) :: nobs_f_all
      ! PE-local model state (dim_p)
      REAL(c_double), INTENT(inout)    :: state_p(dim_p)
      ! Full observed state for all observation types (nobs_f_all)
      REAL(c_double), INTENT(in) :: obs_f_all(nobs_f_all)

      call PDAFomi_obs_op_adj_interp_lin(thisobs(i_obs), nrows, obs_f_all, state_p)
   END SUBROUTINE c__PDAFomi_obs_op_adj_interp_lin


   SUBROUTINE c__PDAFomi_get_interp_coeff_tri(gpc, oc, icoeff) bind(c)
      ! Coordinates of grid points; dim(3,2)
      REAL(c_double), INTENT(in)    :: gpc(3,2)
      ! 3 rows; each containing lon and lat coordinates
      ! Coordinates of observation; dim(2)
      REAL(c_double), INTENT(in)    :: oc(2)
      ! Interpolation coefficients; dim(3)
      REAL(c_double), INTENT(inout) :: icoeff(3)
      call PDAFomi_get_interp_coeff_tri(gpc, oc, icoeff)
   END SUBROUTINE c__PDAFomi_get_interp_coeff_tri

   SUBROUTINE c__PDAFomi_get_interp_coeff_lin1D(gpc, oc, icoeff) bind(c)
      ! Coordinates of grid points (dim=2)
      REAL(c_double), INTENT(in)    :: gpc(2)
      ! Coordinates of observation
      REAL(c_double), INTENT(in)    :: oc
      ! Interpolation coefficients (dim=2)
      REAL(c_double), INTENT(inout) :: icoeff(2)
      call PDAFomi_get_interp_coeff_lin1D(gpc, oc, icoeff)
   END SUBROUTINE c__PDAFomi_get_interp_coeff_lin1D

   SUBROUTINE c__PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff) bind(c)
      ! Length of icoeff
      INTEGER(c_int), INTENT(in) :: num_gp
      ! Number of dimensions in interpolation
      INTEGER(c_int), INTENT(in) :: n_dim
      ! Coordinates of grid points
      REAL(c_double), INTENT(in)    :: gpc(num_gp,n_dim)
      ! Coordinates of observation
      REAL(c_double), INTENT(in)    :: oc(n_dim)
      ! Interpolation coefficients (num_gp)
      REAL(c_double), INTENT(inout) :: icoeff(num_gp)
      call PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)
   END SUBROUTINE c__PDAFomi_get_interp_coeff_lin


   SUBROUTINE c__PDAFomi_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
                                          init_dim_obs_pdaf, obs_op_pdaf, &
                                          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                                          prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf

      call PDAFomi_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
                                             init_dim_obs_pdaf, obs_op_pdaf, &
                                             cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                                             prepoststep_pdaf, next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_assimilate_3dvar

   SUBROUTINE c__PDAFomi_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
                                                   init_dim_obs_pdaf, obs_op_pdaf, &
                                                   cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                                                   prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf

      call PDAFomi_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
                                            init_dim_obs_pdaf, obs_op_pdaf, &
                                            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                                            prepoststep_pdaf, next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_assimilate_en3dvar_estkf

   SUBROUTINE c__PDAFomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf

      call PDAFomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
         g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_assimilate_en3dvar_lestkf

   SUBROUTINE c__PDAFomi_assimilate_global(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAFomi_assimilate_global(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_next_observation, flag)
   END SUBROUTINE c__PDAFomi_assimilate_global

   SUBROUTINE c__PDAFomi_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_pdaf, obs_op_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Apply ensemble control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint ensemble control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf

      call PDAFomi_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
         init_dim_obs_pdaf, obs_op_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
         obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_assimilate_hyb3dvar_estkf

   SUBROUTINE c__PDAFomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: distribute_state_pdaf
      ! Provide time step, time and dimension of next observation
      procedure(c__next_observation_pdaf) :: next_observation_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf

      call PDAFomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
         g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_assimilate_hyb3dvar_lestkf

   SUBROUTINE c__PDAFomi_assimilate_lenkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_localize, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: U_localize
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAFomi_assimilate_lenkf(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_localize, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAFomi_assimilate_lenkf

   SUBROUTINE c__PDAFomi_assimilate_local(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
         U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAFomi_assimilate_local(U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
         U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_next_observation, flag)
   END SUBROUTINE c__PDAFomi_assimilate_local


   SUBROUTINE c__PDAFomi_generate_obs(U_collect_state, U_distribute_state, &
         U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_prepoststep, &
         U_next_observation, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: U_get_obs_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation

      CALL PDAFomi_generate_obs(U_collect_state, U_distribute_state, &
         U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_prepoststep, &
         U_next_observation, flag)
   END SUBROUTINE c__PDAFomi_generate_obs


   SUBROUTINE c__PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
      cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf

      call PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
         cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_put_state_3dvar

   SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf(collect_state_pdaf, &
      init_dim_obs_pdaf, obs_op_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      call PDAFomi_put_state_en3dvar_estkf(collect_state_pdaf, &
         init_dim_obs_pdaf, obs_op_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
         prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_put_state_en3dvar_estkf

   SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf
      call PDAFomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
         g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_put_state_en3dvar_lestkf

   SUBROUTINE c__PDAFomi_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
         U_get_obs_f, U_prepoststep, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_f_pdaf) :: U_init_dim_obs_f
      ! Observation operator
      procedure(c__obs_op_f_pdaf) :: U_obs_op_f
      ! Provide observation vector to user
      procedure(c__get_obs_f_pdaf) :: U_get_obs_f
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

       CALL PDAFomi_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
         U_get_obs_f, U_prepoststep, flag)
   END SUBROUTINE c__PDAFomi_put_state_generate_obs

   SUBROUTINE c__PDAFomi_put_state_global(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_prepoststep, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAFomi_put_state_global(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_prepoststep, flag)
   END SUBROUTINE c__PDAFomi_put_state_global

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf(collect_state_pdaf, &
      init_dim_obs_pdaf, obs_op_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
      obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: init_dim_obs_pdaf
      ! Observation operator
      procedure(c__obs_op_pdaf) :: obs_op_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Apply ensemble control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint ensemble control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf

      call PDAFomi_put_state_hyb3dvar_estkf(collect_state_pdaf, &
         init_dim_obs_pdaf, obs_op_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
         obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_estkf

   SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
      init_dim_obs_f_pdaf, obs_op_f_pdaf, &
      cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
      init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: outflag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: collect_state_pdaf
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: prepoststep_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_ens_pdaf) :: cvt_ens_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_ens_pdaf) :: cvt_adj_ens_pdaf
      ! Apply control vector transform matrix to control vector
      procedure(c__cvt_pdaf) :: cvt_pdaf
      ! Apply adjoint control vector transform matrix
      procedure(c__cvt_adj_pdaf) :: cvt_adj_pdaf
      ! Linearized observation operator
      procedure(c__obs_op_lin_pdaf) :: obs_op_lin_pdaf
      ! Adjoint observation operator
      procedure(c__obs_op_adj_pdaf) :: obs_op_adj_pdaf
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: init_n_domains_pdaf
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: init_dim_l_pdaf
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: g2l_state_pdaf
      ! Init full state from local state
      procedure(c__l2g_state_pdaf) :: l2g_state_pdaf
      ! Initialize dimension of full observation vector
      procedure(c__init_dim_obs_f_pdaf) :: init_dim_obs_f_pdaf
      ! Full observation operator
      procedure(c__obs_op_f_pdaf) :: obs_op_f_pdaf
      ! Initialize local dimimension of obs. vector
      procedure(c__init_dim_obs_l_pdaf) :: init_dim_obs_l_pdaf

      call PDAFomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
         init_dim_obs_f_pdaf, obs_op_f_pdaf, &
         cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
         init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
         g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
   END SUBROUTINE c__PDAFomi_put_state_hyb3dvar_lestkf

   SUBROUTINE c__PDAFomi_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_prepoststep, U_localize, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
      ! Apply localization to HP and HPH^T
      procedure(c__localize_covar_pdaf) :: U_localize

      CALL PDAFomi_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
         U_prepoststep, U_localize, flag)
   END SUBROUTINE c__PDAFomi_put_state_lenkf

   SUBROUTINE c__PDAFomi_put_state_local(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, flag) bind(c)
      ! Status flag
      INTEGER(c_int), INTENT(out) :: flag
      ! Routine to collect a state vector
      procedure(c__collect_state_pdaf) :: U_collect_state
      ! Observation operator
      procedure(c__obs_op_pdaf) :: U_obs_op
      ! Provide number of local analysis domains
      procedure(c__init_n_domains_p_pdaf) :: U_init_n_domains_p
      ! Init state dimension for local ana. domain
      procedure(c__init_dim_l_pdaf) :: U_init_dim_l
      ! Initialize dimension of observation vector
      procedure(c__init_dim_obs_pdaf) :: U_init_dim_obs
      ! Initialize dim. of obs. vector for local ana. domain
      procedure(c__init_dim_obs_l_pdaf) :: U_init_dim_obs_l
      ! Get state on local ana. domain from full state
      procedure(c__g2l_state_pdaf) :: U_g2l_state
      ! Init full state from state on local analysis domain
      procedure(c__l2g_state_pdaf) :: U_l2g_state
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep

      CALL PDAFomi_put_state_local(U_collect_state, U_init_dim_obs, U_obs_op, &
         U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
         U_g2l_state, U_l2g_state, flag)
   END SUBROUTINE c__PDAFomi_put_state_local

   SUBROUTINE c__PDAFomi_init_obs_f_cb(step, dim_obs_f, observation_f) bind(c)
      external PDAFomi_init_obs_f_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of full observation vector
      integer(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector
      real(c_double), INTENT(out)   :: observation_f(dim_obs_f)
      call PDAFomi_init_obs_f_cb(step, dim_obs_f, observation_f)
   END SUBROUTINE c__PDAFomi_init_obs_f_cb

   SUBROUTINE c__PDAFomi_init_obsvar_cb(step, dim_obs_p, obs_p, meanvar) bind(c)
      external PDAFomi_init_obsvar_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! PE-local dimension of observation vector
      integer(c_int), INTENT(in) :: dim_obs_p
      ! PE-local observation vector
      real(c_double), INTENT(in) :: obs_p(dim_obs_p)
      ! Mean observation error variance
      real(c_double), INTENT(out)   :: meanvar
      call PDAFomi_init_obsvar_cb(step, dim_obs_p, obs_p, meanvar)
   END SUBROUTINE c__PDAFomi_init_obsvar_cb

   SUBROUTINE c__PDAFomi_g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
         ostate_l) bind(c)
      external PDAFomi_g2l_obs_cb
      ! Index of current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of full PE-local observation vector
      integer(c_int), INTENT(in) :: dim_obs_f
      ! Dimension of local observation vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! Full PE-local obs.ervation vector
      real(c_double), INTENT(in)    :: ostate_f(dim_obs_f)
      ! Observation vector on local domain
      real(c_double), INTENT(out)   :: ostate_l(dim_obs_l)
      call PDAFomi_g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
      ostate_l)
   END SUBROUTINE c__PDAFomi_g2l_obs_cb

   SUBROUTINE c__PDAFomi_init_obs_l_cb(domain_p, step, dim_obs_l, observation_l) bind(c)
      external PDAFomi_init_obs_l_cb
      ! Index of current local analysis domain index
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Local dimension of observation vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! Local observation vector
      real(c_double), INTENT(out)   :: observation_l(dim_obs_l)
      call PDAFomi_init_obs_l_cb(domain_p, step, dim_obs_l, observation_l)
   END SUBROUTINE c__PDAFomi_init_obs_l_cb

   SUBROUTINE c__PDAFomi_init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l, meanvar_l) bind(c)
      external PDAFomi_init_obs_l_cb
      ! Index of current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Local dimension of observation vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! Local observation vector
      real(c_double), INTENT(in) :: obs_l(dim_obs_l)
      ! Mean local observation error variance
      real(c_double), INTENT(out)   :: meanvar_l
      call PDAFomi_init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l, meanvar_l)
   END SUBROUTINE c__PDAFomi_init_obsvar_l_cb

   SUBROUTINE c__PDAFomi_prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l) bind(c)
      external PDAFomi_prodRinvA_l_cb
      ! Index of current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of local observation vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! Rank of initial covariance matrix
      integer(c_int), INTENT(in) :: rank
      ! Local vector of observations
      real(c_double), INTENT(in)    :: obs_l(dim_obs_l)
      ! Input matrix
      real(c_double), INTENT(inout) :: A_l(dim_obs_l, rank)
      ! Output matrix
      real(c_double), INTENT(out)   :: C_l(dim_obs_l, rank)
      call PDAFomi_prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
   END SUBROUTINE c__PDAFomi_prodRinvA_l_cb

   SUBROUTINE c__PDAFomi_likelihood_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l) bind(c)
      external PDAFomi_likelihood_l_cb
      ! Current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! PE-local vector of observations
      real(c_double), INTENT(in)    :: obs_l(dim_obs_l)
      ! Input vector of residuum
      real(c_double), INTENT(inout) :: resid_l(dim_obs_l)
      ! Output vector - log likelihood
      real(c_double), INTENT(out)   :: lhood_l
      call PDAFomi_likelihood_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l)
   END SUBROUTINE c__PDAFomi_likelihood_l_cb

   SUBROUTINE c__PDAFomi_prodRinvA_cb(step, dim_obs_p, ncol, obs_p, A_p, C_p) bind(c)
      external PDAFomi_prodRinvA_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of PE-local observation vector
      integer(c_int), INTENT(in) :: dim_obs_p
      ! Number of columns in A_p and C_p
      integer(c_int), INTENT(in) :: ncol
      ! PE-local vector of observations
      real(c_double), INTENT(in)    :: obs_p(dim_obs_p)
      ! Input matrix
      real(c_double), INTENT(in)    :: A_p(dim_obs_p, ncol)
      ! Output matrix
      real(c_double), INTENT(out)   :: C_p(dim_obs_p, ncol)
      call PDAFomi_prodRinvA_cb(step, dim_obs_p, ncol, obs_p, A_p, C_p)
   END SUBROUTINE c__PDAFomi_prodRinvA_cb

   SUBROUTINE c__PDAFomi_likelihood_cb(step, dim_obs, obs, resid, lhood) bind(c)
      external PDAFomi_prodRinvA_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      integer(c_int), INTENT(in) :: dim_obs
      ! PE-local vector of observations
      real(c_double), INTENT(in)    :: obs(dim_obs)
      ! Input vector of residuum
      real(c_double), INTENT(in)    :: resid(dim_obs)
      ! Output vector - log likelihood
      real(c_double), INTENT(out)   :: lhood
      call PDAFomi_likelihood_cb(step, dim_obs, obs, resid, lhood)
   END SUBROUTINE c__PDAFomi_likelihood_cb

   SUBROUTINE c__PDAFomi_add_obs_error_cb(step, dim_obs_p, C_p) bind(c)
      external PDAFomi_add_obs_error_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of PE-local observation vector
      integer(c_int), INTENT(in) :: dim_obs_p
      ! Matrix to which R is added
      real(c_double), INTENT(inout) :: C_p(dim_obs_p,dim_obs_p)
      call PDAFomi_add_obs_error_cb(step, dim_obs_p, C_p)
   END SUBROUTINE c__PDAFomi_add_obs_error_cb


   SUBROUTINE c__PDAFomi_init_obscovar_cb(step, dim_obs, dim_obs_p, covar, m_state_p, &
         isdiag) bind(c)
      external PDAFomi_init_obscovar_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of observation vector
      integer(c_int), INTENT(in) :: dim_obs
      ! PE-local dimension of obs. vector
      integer(c_int), INTENT(in) :: dim_obs_p
      ! Observation error covar. matrix
      real(c_double), INTENT(out) :: covar(dim_obs,dim_obs)
      ! Observation vector
      real(c_double), INTENT(in) :: m_state_p(dim_obs_p)
      ! Whether matrix R is diagonal
      LOGICAL(c_bool), INTENT(out) :: isdiag
      call PDAFomi_init_obscovar_cb(step, dim_obs, dim_obs_p, covar, m_state_p, &
         isdiag)
   END SUBROUTINE c__PDAFomi_init_obscovar_cb

   SUBROUTINE c__PDAFomi_init_obserr_f_cb(step, dim_obs_f, obs_f, obserr_f) bind(c)
      external PDAFomi_init_obserr_f_cb
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Full dimension of observation vector
      integer(c_int), INTENT(in) :: dim_obs_f
      ! Full observation vector
      real(c_double), INTENT(in) :: obs_f(dim_obs_f)
      ! Full observation error stddev
      real(c_double), INTENT(out)   :: obserr_f(dim_obs_f)
      call PDAFomi_init_obserr_f_cb(step, dim_obs_f, obs_f, obserr_f)
   END SUBROUTINE c__PDAFomi_init_obserr_f_cb

   SUBROUTINE c__PDAFomi_prodRinvA_hyb_l_cb(domain_p, step, dim_obs_l, rank, obs_l, alpha, A_l, C_l) bind(c)
      ! Index of current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! Dimension of local observation vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! Rank of initial covariance matrix
      integer(c_int), INTENT(in) :: rank
      ! Local vector of observations
      real(c_double),  INTENT(in)    :: obs_l(dim_obs_l)
      ! Hybrid weight
      real(c_double),  INTENT(in)    :: alpha
      ! Input matrix
      real(c_double),  INTENT(inout) :: A_l(dim_obs_l, rank)
      ! Output matrix
      real(c_double),  INTENT(out)   :: C_l(dim_obs_l, rank)
      call PDAFomi_prodRinvA_hyb_l_cb(domain_p, step, dim_obs_l, rank, obs_l, alpha, A_l, C_l)
   END SUBROUTINE c__PDAFomi_prodRinvA_hyb_l_cb

   SUBROUTINE c__PDAFomi_likelihood_hyb_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, alpha, lhood_l) bind(c)
      ! Current local analysis domain
      integer(c_int), INTENT(in) :: domain_p
      ! Current time step
      integer(c_int), INTENT(in) :: step
      ! PE-local dimension of obs. vector
      integer(c_int), INTENT(in) :: dim_obs_l
      ! PE-local vector of observations
      real(c_double),  INTENT(in)    :: obs_l(dim_obs_l)
      ! Input vector of residuum
      real(c_double),  INTENT(inout) :: resid_l(dim_obs_l)
      ! Hybrid weight
      real(c_double),  INTENT(in)    :: alpha
      ! Output vector - log likelihood
      real(c_double),  INTENT(out)   :: lhood_l
      call PDAFomi_likelihood_hyb_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, alpha, lhood_l)
   END SUBROUTINE c__PDAFomi_likelihood_hyb_l_cb

end module PDAFomi_obs_c_binding
