import numpy as np
import pyPDAF.PDAF.PDAFomi as PDAFomi

class OBS:
    n_obs = 0
    
    def __init__(self, typename, mype_filter, 
                        nx, doassim, delt_obs, rms_obs):
        OBS.n_obs += 1

        self.i_obs = OBS.n_obs

        assert OBS.n_obs >= 1, 'observation count must start from 1'

        if (mype_filter==0):
            print(('Assimilate observations:', typename))

        self.doassim = doassim
        self.delt_obs = delt_obs
        self.rms_obs = rms_obs

        # Specify type of distance computation
        # 0=Cartesian 1=Cartesian periodic
        self.disttype = 0

        # Number of coordinates used for distance computation
        # The distance compution starts from the first row
        self.ncoord = len(nx)

        # Allocate process-local index array
        # This array has as many rows as required 
        # for the observation operator
        # 1 if observations are at grid points; 
        # >1 if interpolation is required
        self.nrows = 1

        # Size of domain for periodicity for disttype=1 
        # (<0 for no periodicity)
        self.domainsize = np.zeros(self.ncoord)
        self.domainsize[0] = nx[1]
        self.domainsize[1] = nx[0]

        # Type of observation error: (0) Gauss, (1) Laplace
        self.obs_err_type = None

        # Whether to use (1) global full obs. 
        # (0) obs. restricted to those relevant for a process domain
        self.use_global_obs = 1

        self.icoeff_p = None

    def init_dim_obs(self, step, dim_obs, local_range, 
                                mype_filter, nx, nx_p):

        obs_field = self.get_obs_field(step, nx)

        # Count valid observations that 
        # lie within the process sub-domain
        pe_start = nx_p[-1]*mype_filter
        pe_end = nx_p[-1]*(mype_filter+1)
        obs_field_p = obs_field[:, pe_start:pe_end]
        assert tuple(nx_p) == obs_field_p.shape, \
                'observation decomposition should be the same as the model decomposition'
        cnt_p = np.count_nonzero(obs_field_p > -999.0)
        self.dim_obs_p = cnt_p

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations 
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(nx_p, obs_field_p)
            self.set_id_obs_p(nx_p, obs_field_p)
            self.set_ocoord_p(obs_field_p, pe_start)
            self.set_ivar_obs_p()
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi(local_range)

    def set_obs_p(self, nx_p, obs_field_p):
        obs_field_tmp = obs_field_p.reshape(np.prod(nx_p), order='F')
        self.obs_p = np.zeros(self.dim_obs_p)
        self.obs_p[:self.dim_obs_p] = obs_field_tmp[obs_field_tmp > -999]

    def set_id_obs_p(self, nx_p, obs_field_p):
        self.id_obs_p = np.zeros((self.nrows, self.dim_obs_p))
        obs_field_tmp = obs_field_p.reshape(np.prod(nx_p), order='F')
        cnt0_p = np.where(obs_field_tmp > -999)[0] + 1
        assert len(cnt0_p) == self.dim_obs_p, 'dim_obs_p should equal cnt0_p'
        self.id_obs_p[0, :self.dim_obs_p] = cnt0_p

    def set_ocoord_p(self, obs_field_p, offset):
        self.ocoord_p = np.zeros((self.ncoord, self.dim_obs_p))
        self.ocoord_p[0, :self.dim_obs_p] = np.where(obs_field_p.T > -999)[0] + 1 + offset
        self.ocoord_p[1, :self.dim_obs_p] = np.where(obs_field_p.T > -999)[1] + 1 + offset

    def set_ivar_obs_p(self):
        self.ivar_obs_p = np.ones(
                                self.dim_obs_p
                                )/(self.rms_obs*self.rms_obs)

    def get_obs_field(self, step, nx):
        obs_field = np.zeros(nx)
        obs_field = np.loadtxt(f'inputs_online/obs_step{step}.txt')
        if self.i_obs == 1:
            # Make the observations at (8,5), (12,15) and (4,30) invalid
            # They will be used in observation type B
            obs_field[7, 4] = -1000.0 
            obs_field[11, 14] = -1000.0 
            obs_field[3, 29] = -1000.0
        else:
            # Make the observations at (8,5), (12,15) and (4,30) invalid
            # They will be used in observation type B
            obs_tmp = [obs_field[7, 4], obs_field[11, 14], obs_field[3, 29]]
            obs_field[:] = -1000.
            obs_field[7, 4] = obs_tmp[0]
            obs_field[11, 14] = obs_tmp[1]
            obs_field[3, 29] = obs_tmp[2]

        return obs_field

    def set_PDAFomi(self, local_range):
        PDAFomi.setOMIessential(self.i_obs, self.doassim, self.disttype, 
                                self.ncoord, self.id_obs_p)
        if self.domainsize is not None:
            PDAFomi.set_domainsize(self.i_obs, self.domainsize)
        if self.obs_err_type is not None:
            PDAFomi.set_obs_err_type(self.i_obs, self.obs_err_type)
        if self.use_global_obs is not None:
            PDAFomi.set_use_global_obs(self.i_obs, self.use_global_obs)
        if self.icoeff_p is not None:
            PDAFomi.set_icoeff_p(self.i_obs, self.icoeff_p)
        self.dim_obs = PDAFomi.gather_obs(self.i_obs, 
                                          self.dim_obs_p,
                                          self.nrows,
                                          self.obs_p, 
                                          self.ivar_obs_p, 
                                          self.ocoord_p, 
                                          local_range)

    def obs_op(self, step, state_p, ostate):
        if (self.doassim == 1):
            PDAFomi.obs_op_gridpoint(self.i_obs, state_p, ostate)

    def init_dim_obs_l(self, localization, domain_p, step, dim_obs, dim_obs_l):
        return PDAFomi.init_dim_obs_l(self.i_obs, localization.coords_l, 
                                localization.loc_weight, 
                                localization.local_range, 
                                localization.srange)

    def localize_covar(self, localization, HP_p, HPH, coords_p):
        PDAFomi.localize_covar(self.i_obs, localization.loc_weight, 
                                localization.local_range, 
                                localization.srange, 
                                coords_p, hp_p, hph)

    def deallocate_obs(self, step):
        PDAFomi.deallocate_obs(self.i_obs, step)