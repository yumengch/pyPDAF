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

import config_obsA
import config_obsB
import model

class Distributor:
    """This class implements the function where
    PDAF distributes ensemble to the model field

    Attributes
    ----------
    model: model.model
        model instance
    """
    def __init__(self, model_t:model.Model) -> None:
        # get the model insta
        self.model = model_t

    def distribute_state_ini(self, _dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """PDAF will distribute state vector (state_p) to model field

        Parameters
        ----------
        dim_p: int
            Dimension of the state vector on local processor
        state_p: np.ndarray
            state vector on local processor

        Returns
        -------
        state_p: np.ndarray
            state vector
        """
        self.model.field_p[:] = state_p.reshape(self.model.ny_p, self.model.nx_p)
        return state_p


    def distribute_state(self, _dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """PDAF will distribute state vector (state_p) to model field

        Parameters
        ----------
        dim_p: int
            Dimension of the state vector on local processor
        state_p: np.ndarray
            state vector on local processor

        Returns
        -------
        state_p: np.ndarray
            state vector
        """
        self.model.field_p[:] = state_p.reshape(self.model.ny_p, self.model.nx_p)
        print ("distribute state to model", state_p[:6])
        return state_p

    def next_observation(self, _stepnow:int, nsteps:int,
                         doexit:int, time:float) -> tuple[int, int, float]:
        """Providing PDAF the information on the number of model integration steps
        to next analysis
        """
        # next observation will arrive at `nsteps' steps
        nsteps = min(config_obsA.dtobs, config_obsB.dtobs)
        # doexit = 0 means that PDAF will continue to distribute state
        # to model for further integrations
        doexit = 0
        # model time is not used here as we only use steps to define the time
        return nsteps, doexit, time
