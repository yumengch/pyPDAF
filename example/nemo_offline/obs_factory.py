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
import configparser
import log

import numpy as np

import model
import obs_sst
import sfields

class ObsFactory:
    """Create observation modules and expose PDAF-OMI callbacks.

    ``[observations].obsnames`` selects which observation sections are active.
    The current factory implements ``sst`` through :class:`obs_sst.ObsSST`.
    To add another observation type, implement a class with the same callback
    methods and register it in this factory.

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
    def __init__(self, config: configparser.ConfigParser, model_grid:model.NemoDomain,
                 fields:dict[str, sfields.StateField]) -> None:
        # Initialise observations
        self.model = model_grid
        self.fields = fields
        self.obs_list = []
        obs_names = [v.strip() for v in config['observations']['obsnames'].split(",")]
        self.nobs = 0
        for obs_name in obs_names:
            # Register observation modules here. The PDAF-OMI observation index
            # starts at 1 and must increase for every assimilated type.
            if obs_name == 'sst':
                self.nobs += 1
                self.obs_list.append(obs_sst.ObsSST(self.nobs, config[obs_name]))
            # Add more observation types here as needed
        output_str = f'total number of observation types: {self.nobs}'
        log.info (output_str)

    def init_dim_obs_pdafomi(self, _step:int, dim_obs:int) -> int:
        """Initialise global observation dimensions for PDAF-OMI.

        Each observation module creates or reads its local observations,
        registers metadata with PDAF-OMI, and contributes to the total
        observation dimension.

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
        # calculate the dimension of full observation vector
        dim_obs = 0
        for obs in self.obs_list:
            if obs.doassim == 1:
                dim_obs_o = obs.init_dim(dim_obs, self.fields[obs.field_name], self.model)
                dim_obs = dim_obs + dim_obs_o
        return dim_obs

    def obs_op_pdafomi(self, step:int, _dim_p:int, _dim_obs_p:int,
                       state_p:np.ndarray, ostate:np.ndarray) -> np.ndarray:
        """Map a state vector to observation space.

        For the synthetic SST example, PDAF-OMI's grid-point observation
        operator simply selects state entries. More complex observations can
        interpolate, integrate, or otherwise transform the model state.

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

    def init_dim_obs_l_pdafomi(self, domain_p:int, _step:int, _dim_obs:int, dim_obs_l:int) -> int:
        """Initialise the number of observations used by one local domain.

        Local filters call this for every wet surface point. The current
        coordinates are longitude/latitude in radians, matching PDAF-OMI's
        spherical distance option.

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
        if self.model.wet_grid.nwet2d == 0:
            return 0
        j = self.model.wet_grid.wet_pts[5, domain_p - 1]
        i = self.model.wet_grid.wet_pts[6, domain_p - 1]
        coords_l = np.array([np.deg2rad(self.model.coords.nav_lon[j, i]),
                             np.deg2rad(self.model.coords.nav_lat[j, i])], dtype=float)

        dim_obs_l = 0
        for obs in self.obs_list:
            if obs.doassim == 1:
                dim_obs_l_o = obs.init_dim_obs_l(coords_l, dim_obs_l)
                dim_obs_l = dim_obs_l + dim_obs_l_o
        return dim_obs_l
