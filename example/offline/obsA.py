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
import log

import numpy as np
import pyPDAF.PDAF as PDAF

import config_obsA as config
import localisation
import model
import parallelisation

class obsA:

    """observation functions for type-A observation

    Attributes
    ----------
    i_obs : int
        index of the observation type
    obs_name : str
        name of the current observation type
    filename : str
        observation filename to be read
    dim_obs_p : int
        dimension size of the PE-local observation vector
    disttype : int
        type of distance computation to use for localization
    doassim : int
        whether to assimilate this observation type
    ncoord : int
        number of coordinate dimension
    nrows : int
        number of rows in ocoord_p
    dtobs : int
        time interval between observations
    missing_value : float
        missing value in observations (filled by -9999.0 for example)
    model_t : model.model
        model object
    pe : parallelisation.parallelisation
        parallelization object
    local : localisation.localisation
        localisation object
    """

    def __init__(self, i_obs:int,
                 pe:parallelisation.parallelisation,
                 model_grid:model.model_grid,
                 local:localisation.localisation) -> None:
        # i_obs-th observations in the system starting from 1
        self.i_obs:int = i_obs
        self.obs_name:str = config.obs_name
        assert self.i_obs >= 1, 'observation count must start from 1'
        # observation filename
        self.filename:str = config.obs_path
        # whether this observation is assimilated
        self.doassim:int = 1
        # Type of distance computation to use for localization
        # It is mandatory for OMI even if we don't use localisation
        # 0=Cartesian 1=Cartesian periodic
        # details of this property can be seen at
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype
        # currently, only PDAF V2.1 is supported
        self.disttype:int = config.disttype
        # Number of coordinates use for distance computation
        # Here, it is a 2-dimensional domain
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsncoord
        self.ncoord:int = config.ncoord
        # nrows depends on the necessity of interpolating
        # observations onto model grid
        # if nrows = 1, observations are on the model grid points
        # when interpolation is required,
        # this is the number of grid points required for interpolation.
        # For example, nrows = 4 for bi-linear interpolation in 2D,
        # and 8 for 3D linear interpolation.
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsid_obs_p
        # More information about interpolation is available at
        # https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin
        self.nrows:int = config.nrows
        self.dtobs = config.dtobs
        # specify the missing/filled value of the observations
        self.missing_value:float = config.missing_value
        self.model_grid:model.model_grid = model_grid
        self.pe:parallelisation.parallelisation = pe
        self.local:localisation.localisation = local


    def init_dim(self, step:int, dim_obs:int) -> int:
        """In PDAFomi, init_dim takes the responsibility to
        set the obs_f object, gather information on the observations
        including the observation values, observation error, index of
        the state vector of the observation, the model domain and
        observation coordinate for interpolation, and the distance metric
        for localisation.


        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension of observation vector

        Returns
        -------
        int
            dimension of observation vector
        """
        if self.pe.mype_filter == 0:
            log.logger.info(f'Assimilate observations: {self.obs_name}')

        # switch for assimilation of the observation
        PDAF.omi_set_doassim(self.i_obs, self.doassim)
        # Type of distance computation to use for localization
        # It is mandatory for OMI even if we don't use localisation
        PDAF.omi_set_disttype(self.i_obs, self.disttype)
        # Number of coordinates used for distance computation
        PDAF.omi_set_ncoord(self.i_obs, self.ncoord)

        # read observations
        filename = self.filename
        obs_field:np.ndarray = np.loadtxt(filename)
        # when domain decomposition is used, we only need observations
        # within the local model domain
        pe_start:int = self.model_grid.nx_p*self.pe.mype_filter
        pe_end:int = pe_start + self.model_grid.nx_p
        obs_field_p:np.ndarray = obs_field[:, pe_start:pe_end].ravel()
        # get total number of valid observations
        # a mask for invalid observations
        condition:np.ndarray = np.logical_not(
            np.isclose(obs_field_p, self.missing_value)
            )
        dim_obs_p:int = np.count_nonzero(condition)
        # obtain the observation vector
        obs_p:np.ndarray = obs_field_p[condition]
        assert len(obs_p) == dim_obs_p, 'dimension of the ' \
            f'observation vector ({len(obs_p)})' \
                f'should be the same as the dim_obs_p ({dim_obs_p})'

        # inverse of observation variance
        # here we specify/hard-code the standard deviation of observation is 0.5
        ivar_obs_p:np.ndarray = (1./config.rms_obs/config.rms_obs)*np.ones_like(obs_p)

        # coordinate of each observations
        ocoord_p:np.ndarray = np.zeros((self.ncoord, len(obs_p)))
        ocoord_p[0] = np.tile(np.arange(self.model_grid.nx_p) +
                              pe_start,
                              self.model_grid.ny_p)[condition]
        ocoord_p[1] = np.repeat(np.arange(self.model_grid.ny_p),
                                self.model_grid.nx_p)[condition]
        ocoord_p = ocoord_p + 1.


        id_obs_p:np.ndarray = np.zeros((self.nrows, len(obs_p)),
                                       dtype=np.intc)
        if self.nrows == 1:
            # The relationship between observation and state vector
            # id_obs_p gives the indices of observed field in state vector
            # the index starts from 1
            # The index is based on the full state vector
            # instead of the index in the local domain
            # The following code is used because in our example,
            # observations are masked and
            # have the same shape as the model grid
            state_vector_index_p:np.ndarray = \
                np.arange(1,
                          self.model_grid.nx_p*self.model_grid.ny_p +
                            1, dtype=np.intc)
            id_obs_p[0] = state_vector_index_p[condition]
        else:
            # If interpolation is required for a 2D domain
            # id_obs_p has a dimension of (4, dim_obs_p)
            # id_obs_p[0] is the index of the grid point in the state vector at lower left of the observation
            # for more details of the interpolation see:
            # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsid_obs_p
            # and
            # https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin
            # In this case, we also need to specify the coefficients for linear interpolation
            # using PDAF.omi_get_interp_coeff_lin() function and set icoeff_p to PDAF
            icoeff_p:np.ndarray = np.zeros((self.nrows, len(obs_p)))
            for i in range(dim_obs_p):
                # here gcoords are set to 0 in this example
                # in real applications, it must be the actual coordinates
                gcoords:np.ndarray = np.zeros((self.nrows, self.ncoord))
                icoeff_p[:, i] = PDAF.omi_get_interp_coeff_lin(gcoords, ocoord_p[:, i], icoeff_p[:, i])
            PDAF.omi_set_icoeff_p(self.i_obs, icoeff_p)
        PDAF.omi_set_id_obs_p(self.i_obs, id_obs_p)

        # Type of observation error: (0) Gauss, (1) Laplace
        # This is optional
        # Without explicit setting, this is 0
        PDAF.omi_set_obs_err_type(self.i_obs, 0)

        # Whether to use (1) global full obs.
        # (0) obs. restricted to those relevant for a process domain
        # Without explicit setting, this is 1 (using all obs.)
        PDAF.omi_set_use_global_obs(self.i_obs, 1)

        # when localisation is used we need to
        if self.local.local_filter:
            # Size of domain for periodicity for disttype=1
            # (<0 for no periodicity)
            domainsize:np.ndarray = np.array([self.model_grid.nx,
                                              self.model_grid.ny],
                                              dtype=float)
            PDAF.omi_set_domainsize(self.i_obs, domainsize)

        # PDAF need to gather observation information
        dim_obs = PDAF.omi_gather_obs(self.i_obs,
                                      obs_p,
                                      ivar_obs_p,
                                      ocoord_p,
                                      self.local.cradius)
        return dim_obs

    def init_dim_obs_l(self, domain_p:int, step:int, dim_obs:int, dim_obs_l:int) -> int:
        """intialise local observation vector for domain localisation

        Parameters
        ----------
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
        # initialize coordinates of local domain
        # we use grid point indices as coordinates,
        #  but could e.g. use meters
        coords_l:np.ndarray = np.zeros(self.ncoord)
        # we comment out the smart way to calculate the index
        # offset = self.pe.mype_filter*self.model.nx_p*self.model.ny_p
        # coords_l[0] = domain_p + offset - 1
        # coords_l[0] = coords_l[0]//self.model.ny_p
        # coords_l[1] = domain_p + offset - 1
        # coords_l[1] = coords_l[1] - coords_l[0]*self.model.ny_p

        # here is a brutal force way where we simply list all coordinates on each local processor
        # and index them based on domain_p
        offset:int = self.pe.mype_filter*self.model_grid.nx_p
        coords_l[0] = np.tile(np.arange(self.model_grid.nx_p) +
                              offset,
                              self.model_grid.ny_p)[domain_p - 1]
        coords_l[1] = np.repeat(np.arange(self.model_grid.ny_p),
                                self.model_grid.nx_p)[domain_p - 1]
        coords_l = coords_l + 1.

        return PDAF.omi_init_dim_obs_l_iso(self.i_obs, coords_l,
                                      self.local.loc_weight,
                                      self.local.cradius,
                                      self.local.sradius, dim_obs_l)

    def localize_covar(self, dim_p:int, dim_obs:int, HP_p:np.ndarray, HPH:np.ndarray, coords_p:np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """get localised covariance matrix for covariance localisation

        Parameters
        ----------
        dim_p: int
            dimension of state vector on local processor
        dim_obs_p: int
            dimension of observation vector on local processor
        HP_p : ndarray
            matrix HPH
        HPH : ndarray
            PE local part of matrix HP
        coords_p : ndarray
            coordinates of state vector elements

        Returns
        -------
        HP_p : ndarray
            matrix HP
        HPH : ndarray
            matrix HPH
        """
        return PDAF.omi_localize_covar_iso(self.i_obs, self.local.loc_weight,
                                       self.local.cradius,
                                       self.local.sradius,
                                       coords_p, HP_p, HPH)


    def obs_op(self, step:int, state_p:np.ndarray, ostate:np.ndarray) -> np.ndarray:
        """convert state vector by observation operator

        Parameters
        ----------
        step : int
            current time step
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix

        Returns
        -------
        ostate : ndarray
            state vector transformed by identity matrix
        """
        if self.nrows == 1:
            return PDAF.omi_obs_op_gridpoint(self.i_obs, state_p, ostate)
        else:
            # if interpolation is required
            return PDAF.omi_obs_op_interp_lin(self.i_obs, self.nrows, state_p, ostate)