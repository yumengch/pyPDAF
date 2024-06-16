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


class AssimilationDimensions:

    """Dimension of state vector and ensemble size

    Attributes
    ----------
    dim_ens : int
        ensemble size
    dim_state : int
        dimension of global state vector
    dim_state_p : int
        dimension of PE-local state vector
    """

    def __init__(self, model, dim_ens):
        """AssimilationDimensions constructor

        Parameters
        ----------
        model : `Model.Model`
            model object
        dim_ens : int
            ensemble size
        """
        self.dim_state_p = np.prod(model.nx_p)
        self.dim_state = np.prod(model.nx)
        self.dim_ens = dim_ens
