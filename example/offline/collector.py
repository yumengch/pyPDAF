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
    """Here, the background/forecast ensemble is read from files.

    Attributes
    ----------
    model: model.model_grid
        model instance
    """
    def __init__(self, model_grid:model.model_grid,
                 pe: parallelisation.parallelisation) -> None:
        # initialise the model instance
        self.model_grid: model.model_grid = model_grid
        self.pe: parallelisation.parallelisation = pe

    def init_ens_pdaf(self, filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray, ens_p:np.ndarray,
                      status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """read forecast ensemble for DA

        In offline DA, the ensemble is read from disk.
        One does not need collect/distribute ensemble to model.

        Parameters
        ----------
        filtertype: int
            filter type
        dim_p: int
            dimension of state vector on local processor
        dim_ens: int
            ensemble size
        state_p: np.ndarray
            state vector on local processor
        uinv: np.ndarray
            inverse of of
        ens_p: np.ndarray
            ensemble matrix on local processor
        """
        # The initial ensemble is read here and will be distributed to
        # the model in the PDAF.get_state functtion by a distributor.

        # If your ensemble is read from a restart file, you can simply set this
        # function as a dummy function without doing anything
        # However, you still need to set a distributor to call PDAF.get_state, which
        # does nothing as well.
        nx_p:int = self.model_grid.nx_p
        offset:int = self.pe.mype_filter*nx_p
        for i in range(dim_ens):
            ens_p[:, i] = np.loadtxt(
                config.init_ens_path.format(i=i+1)
                )[:, offset:offset+nx_p].ravel()
        return state_p, uinv, ens_p, status_pdaf

    def collect_state(self, dim_p:int, state_p:np.ndarray
                      ) -> np.ndarray:
        """PDAF will collect state vector (state_p) from model field.


        Parameters
        ----------
        dim_p: int
            Dimension of the state vector on local processor
        state_p: np.ndarray
            state vector on local processor.
            This argument is used by PDAF to form the state vector
            of a single ensemble member.

        Returns
        -------
        state_p: np.ndarray
            state vector filled with model field
        """
        return state_p