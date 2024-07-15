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
import logging

import numpy as np
import pyPDAF.PDAF as PDAF

import config_obsA
import config_obsB
import localisation
import model
import obsA
import obsB
import parallelisation

class obs_factory:
    """This class implements all user-supplied functions
    used by PDAFomi. These functions are called at every time steps

    Attributes
    ----------
    pe : `parallelisation.parallelisation`
        parallelization object
    model : `model.model`
        model object
    local : `localisation.localisation`
        localisation object
    obs_list : list
        list of observation types
    nobs : int
        total number of observation types
    """
    def __init__(self, pe:parallelisation.parallelisation, model_t:model.model, local:localisation.localisation) -> None:
        # Initialise observations
        self.pe: parallelisation.parallelisation =  pe
        self.model:model.model = model_t
        self.local:localisation.localisation = local
        self.obs_list:list = []
        self.nobs:int = 0
        if config_obsA.assim:
            self.nobs += 1
            self.obs_list.append(obsA.obsA(self.nobs, self.pe, self.model, self.local)
            )
        if config_obsB.assim:
            self.nobs += 1
            self.obs_list.append(obsB.obsB(self.nobs, self.pe, self.model, self.local)
            )
        logging.info (f'total number of observation types: {self.nobs}')
        PDAF.omi_init(self.nobs)

    def init_dim_obs_pdafomi(self, step:int, dim_obs:int) -> int:
        """initialise observation dimensions

        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension of observation vector

        Returns
        -------
        dim_obs : int
            dimension of observation vector
        """
        # it is possibly useful to add some checks on obs.doassim here
        # For example, set obs.doassim = 0 if one type of observation
        # is not used for this particular step.
        for obs in self.obs_list:
            if step % obs.dtobs == 0:
                obs.doassim = 1
        # calculate the dimension of full observation vector
        dim_obs = 0
        for obs in self.obs_list:
            if obs.doassim == 1:
                dim_obs_o:int = obs.init_dim(step, dim_obs)
                dim_obs = dim_obs + dim_obs_o

        return dim_obs

    def obs_op_pdafomi(self, step:int, dim_p:int, dim_obs_p:int, state_p:np.ndarray, ostate:np.ndarray) -> np.ndarray:
        """turn state vector to observation space

        Parameters
        ----------
        step : int
            current time step
        dim_p: int
            dimension of state vector on local processor
        dim_obs_p: int
            dimension of observation vector on local processor
        state_p : ndarray
            local PE state vector
        ostate : ndarray
            state vector in obs space

        Returns
        -------
        ostate : ndarray
            state vector in obs space
        """

        for obs in self.obs_list:
            if obs.doassim == 1:
                ostate = obs.obs_op(step, state_p, ostate)
        return ostate

    def init_dim_obs_l_pdafomi(self, domain_p:int, step:int, dim_obs:int, dim_obs_l:int) -> int:
        """initialise number of observations in each local domain

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
            dimension of local observation vector
        """
        dim_obs_l = 0
        for obs in self.obs_list:
            if obs.doassim == 1:
                dim_obs_l_o:int = obs.init_dim_obs_l(domain_p, step, dim_obs, dim_obs_l)
                dim_obs_l = dim_obs_l + dim_obs_l_o

        return dim_obs_l

    def localize_covar_pdafomi(self, dim_p:int, dim_obs:int, HP_p:np.ndarray, HPH:np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """getting localised covariance matrix

        Parameters
        ----------
        dim_p: int
            dimension of state vector on local processor
        dim_obs_p: int
            dimension of observation vector on local processor
        HP_p : ndarray
            matrix HP
        HPH : ndarray
            matrix HPH

        Returns
        -------
        HP_p : ndarray
            matrix HP
        HPH : ndarray
            matrix HPH
        """
        # coords_p is the coordinate of the state vector which should have the same
        # unit as the obervation coordinates
        coords_p = np.zeros((2, dim_p))
        offset = self.pe.mype_filter*self.model.nx_p

        coords_p[0] = np.tile(np.arange(self.model.nx_p) + offset, self.model.ny_p)
        coords_p[1] = np.repeat(np.arange(self.model.ny_p), self.model.nx_p)
        coords_p = coords_p + 1

        for obs in self.obs_list:
            obs.localize_covar(dim_p, dim_obs, HP_p, HPH, coords_p)

        return HP_p, HPH
