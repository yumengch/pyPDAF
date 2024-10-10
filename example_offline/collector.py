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

import config
import model
import parallelisation


class collector:
    """This class implements functions where PDAF collects the state vector from model ensemble

    Attributes
    ----------
    model: model.model
        model instance
    """
    def __init__(self, model_t:model.model, pe: parallelisation.parallelisation) -> None:
        # initialise the model instance
        self.model: model.model = model_t
        self.pe: parallelisation = pe

    def init_ens_pdaf(self, filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray, ens_p:np.ndarray,
                      status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """Here, only ens_p variable matters while dim_p and dim_ens defines the
        size of the variables. uinv, state_p are not used in this example.

        status_pdaf is used to handle errors which we will not do it in this example.
        """
        # The initial ensemble is read here and will be distributed to 
        # the model in the PDAF.get_state functtion by a distributor.

        # If your ensemble is read from a restart file, you can simply set this
        # function as a dummy function without doing anything
        # However, you still need to set a distributor to call PDAF.get_state, which
        # does nothing as well. 
        nx_p:int = self.model.nx_p
        offset:int = self.pe.mype_filter*nx_p
        for i in range(dim_ens):
            ens_p[:, i] = np.loadtxt(config.init_ens_path.format(i=i+1))[:, offset:offset+nx_p].ravel()
        return state_p, uinv, ens_p, status_pdaf

    def collect_state(self, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """PDAF will collect state vector (state_p) from model field.


        Parameters
        ----------
        dim_p: int
            Dimension of the state vector on local processor
        state_p: np.ndarray
            state vector on local processor.
            This argument is used by PDAF to form the state vector and ensemble matrix

        Returns
        -------
        state_p: np.ndarray
            state vector filled with model field
        """
        # The [:] treatment ensures that we only change values of state_p not the memory address
        if not config.offline_mode:
            state_p[:] = self.model.field_p.ravel()
        return state_p
