"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and
National Centre for Earth Observation

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import pyPDAF.PDAF as PDAF


class OBS_A:

    """observation information and user-supplied routines

    Attributes
    ----------
    delt_obs : int
        time step interval for observations
    dim_obs : int
        dimension size of the observation vector
    dim_obs_p : int
        dimension size of the PE-local observation vector
    disttype : int
        type of distance computation to use for localization
    doassim : int
        whether to assimilate this observation type
    domainsize : ndarray
        size of domain for periodicity (<=0 for no periodicity)
    i_obs : int
        index of the observation type
    icoeff_p : ndarray
        2d array for interpolation coefficients for obs. operator
    id_obs_p : ndarray
        indices of process-local observed field in state vector
    ivar_obs_p : ndarray
        vector of process-local inverse observation error variance
    n_obs : int
        number of observation types
    ncoord : int
        number of coordinate dimension
    nrows : int
        number of rows in ocoord_p
    obs_err_type : int
        type of observation error
    obs_p : ndarray
        vector of process-local observations
    ocoord_p : ndarray
        2d array of process-local observation coordinates
    rms_obs : float
        observation error standard deviation (for constant errors)
    use_global_obs : int
       Whether to use (1) global full obs. or
       (0) obs. restricted to those relevant for a process domain
    """

    n_obs = 0

    def __init__(self, typename, cnt, mype_filter, task_id,
                 dims, doassim, delt_obs, rms_obs):
        """constructor

        Parameters
        ----------
        typename : string
            name of the observation type
        mype_filter : int
            rank of the PE in filter communicator
        dims : ndarray
            grid size of the model domain
        doassim : int
            whether to assimilate this observation type
        delt_obs : int
            time step interval for observations
        rms_obs : float
            observation error standard deviation (for constant errors)
        """
        OBS_A.n_obs = cnt

        self.i_obs = OBS_A.n_obs

        assert OBS_A.n_obs >= 1, 'observation count must start from 1'

        if (mype_filter == 0 and task_id==1):
            print('Assimilate observations:', typename)

        self.doassim = doassim
        self.delt_obs = delt_obs
        self.rms_obs = rms_obs

        # Specify type of distance computation
        # 0=Cartesian 1=Cartesian periodic
        self.disttype = 0

        # Number of coordinates used for distance computation
        # The distance compution starts from the first row
        self.ncoord = len(dims)

        # Allocate process-local index array
        # This array has as many rows as required
        # for the observation operator
        # 1 if observations are at grid points;
        # >1 if interpolation is required
        self.nrows = 1

        # Size of domain for periodicity for disttype=1
        # (<0 for no periodicity)
        if self.i_obs == 1:
            self.domainsize = np.zeros(self.ncoord)
            self.domainsize[0] = dims[1]
            self.domainsize[1] = dims[0]
        else:
            self.domainsize = None

        # Type of observation error: (0) Gauss, (1) Laplace
        self.obs_err_type = None

        # Whether to use (1) global full obs.
        # (0) obs. restricted to those relevant for a process domain
        self.use_global_obs = 1

        self.icoeff_p = None

    def init_dim_obs(self, step, dim_obs, cradius,
                     mype_filter, dims, dims_p):
        """intialise PDAFomi and getting dimension of observation vector

        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension size of the observation vector
        cradius : float
            cut-off radius for local observation domain
        mype_filter : int
            rank of the PE in filter communicator
        dims : ndarray
            integer array for grid size
        dims_p : ndarray
            integer array for PE-local grid size
        """
        obs_field = self.get_obs_field(step, dims)

        # Count valid observations that
        # lie within the process sub-domain
        pe_start = dims_p[-1]*mype_filter
        pe_end = dims_p[-1]*(mype_filter+1)
        obs_field_p = obs_field[:, pe_start:pe_end]
        assert tuple(dims_p) == obs_field_p.shape, \
               'observation decomposition should be the same as' \
               ' the model decomposition'
        cnt_p = np.count_nonzero(obs_field_p > -999.0)
        self.dim_obs_p = cnt_p

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(dims_p, obs_field_p)
            self.set_id_obs_p(dims_p, obs_field_p)
            self.set_ocoord_p(obs_field_p, pe_start)
            self.set_ivar_obs_p()
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi(cradius)

    def set_obs_p(self, dims_p, obs_field_p):
        """set up PE-local observation vector

        Parameters
        ----------
        dims_p : ndarray
            PE-local model domain
        obs_field_p : ndarray
            PE-local observation field
        """
        obs_field_tmp = obs_field_p.reshape(np.prod(dims_p), order='F')
        self.obs_p = np.zeros(self.dim_obs_p)
        self.obs_p[:self.dim_obs_p] = obs_field_tmp[obs_field_tmp > -999]

    def set_id_obs_p(self, dims_p, obs_field_p):
        """set id_obs_p

        Parameters
        ----------
        dims_p : ndarray
            PE-local model domain
        obs_field_p : ndarray
            PE-local observation field
        """
        self.id_obs_p = np.zeros((self.nrows, self.dim_obs_p))
        obs_field_tmp = obs_field_p.reshape(np.prod(dims_p), order='F')
        cnt0_p = np.where(obs_field_tmp > -999)[0] + 1
        assert len(cnt0_p) == self.dim_obs_p, 'dim_obs_p should equal cnt0_p'
        self.id_obs_p[0, :self.dim_obs_p] = cnt0_p

    def set_ocoord_p(self, obs_field_p, offset):
        """set ocoord_p

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        offset : int
            PE-local offset starting from rank 0
        """
        self.ocoord_p = np.zeros((self.ncoord, self.dim_obs_p))
        ix, iy = np.where(obs_field_p.T > -999)
        self.ocoord_p[0, :self.dim_obs_p] = ix + 1 + offset
        self.ocoord_p[1, :self.dim_obs_p] = iy + 1

    def set_ivar_obs_p(self):
        """set ivar_obs_p
        """
        self.ivar_obs_p = np.ones(
                                self.dim_obs_p
                                )/(self.rms_obs*self.rms_obs)

    def get_obs_field(self, step, dims):
        """retrieve observation field

        Parameters
        ----------
        step : int
            current time step
        dims : ndarray
            grid size of the model domain

        Returns
        -------
        obs_field : ndarray
            observation field
        """
        obs_field = np.zeros(dims)
        obs_field = np.loadtxt(f'inputs_online/obs_step{step}.txt')
        return obs_field

    def set_PDAFomi(self, cradius):
        """set PDAFomi obs_f object

        Parameters
        ----------
        cradius : double
            localization cut-off radius (the maximum radius used in this process domain)
        """
        #print (self.obs_p)
        PDAF.omi_set_doassim(self.i_obs, self.doassim)
        PDAF.omi_set_disttype(self.i_obs, self.disttype)
        PDAF.omi_set_ncoord(self.i_obs, self.ncoord)
        PDAF.omi_set_id_obs_p(self.i_obs, self.id_obs_p)
        if self.domainsize is not None:
            PDAF.omi_set_domainsize(self.i_obs, self.domainsize)
        if self.obs_err_type is not None:
            PDAF.omi_set_obs_err_type(self.i_obs, self.obs_err_type)
        if self.use_global_obs is not None:
            PDAF.omi_set_use_global_obs(self.i_obs, self.use_global_obs)
        if self.icoeff_p is not None:
            PDAF.omi_set_icoeff_p(self.i_obs, self.icoeff_p)

        self.dim_obs = PDAF.omi_gather_obs(self.i_obs,
                                          self.obs_p,
                                          self.ivar_obs_p,
                                          self.ocoord_p,
                                          cradius)

    def obs_op(self, step, state_p, ostate):
        """convert state vector by observation operator

        Parameters
        ----------
        step : int
            current time step
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix
        """
        ostate = PDAF.omi_obs_op_gridpoint(self.i_obs, state_p, ostate)
        return ostate

    def init_dim_obs_l(self, localization, domain_p, step, dim_obs, dim_obs_l):
        """intialise local observation vector

        Parameters
        ----------
        localization : TYPE
            Description
        domain_p : int
            index of current local analysis domain
        step : int
            current time step
        dim_obs : int
            dimension of observation vector
        dim_obs_l : int
            dimension of local observation vector

        Returns
        -------
        dim_obs_l : int
            dimension of local observations
        """
        dim_obs_l = 0
        return PDAF.omi_init_dim_obs_l(self.i_obs, localization.coords_l,
                                      localization.loc_weight,
                                      localization.cradius,
                                      localization.sradius, dim_obs_l)

    def localize_covar(self, localization, HP_p, HPH, coords_p):
        """localze covariance matrix

        Parameters
        ----------
        localization : `Localization.Localization`
            the localization object
        HP_p : ndarray
            matrix HPH
        HPH : ndarray
            PE local part of matrix HP
        coords_p : ndarray
            coordinates of state vector elements
        """
        PDAF.omi_localize_covar(self.i_obs, localization.loc_weight,
                               localization.cradius,
                               localization.sradius,
                               coords_p, HP_p, HPH)

    def deallocate_obs(self):
        """deallocate PDAFomi object

        Parameters
        ----------
        step : int
            current time step
        """
        PDAF.omi_deallocate_obs(self.i_obs)
