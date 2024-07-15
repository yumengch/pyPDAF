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

import model


class collector:
    """This class implements functions where PDAF collects the state vector from model ensemble

    Attributes
    ----------
    model: model.model
        model instance
    """
    def __init__(self, model_t:model.model) -> None:
        # initialise the model instance
        self.model: model.model = model_t

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
        state_p[:] = self.model.field_p.ravel()
        return state_p
