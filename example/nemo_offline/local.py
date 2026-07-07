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

import pyPDAF

import model
import sfields


class Local:
    """Local-analysis callback functions for local PDAF filters.

    Local filters such as LETKF/LESTKF analyse one local domain at a time. In
    this NEMO adapter each wet surface point is one local domain; the local
    state contains all configured 2D variables at that point plus all wet
    vertical levels for configured 3D variables in the same water column.

    Attributes
    ----------
    loc_weight : int
        - (0) constant weight of 1
        - (1) exponentially decreasing with sradius
        - (2) use 5th-order polynomial
        - (3) regulated localization of R with mean error variance
        - (4) regulated localization of R with single-point error variance

    local_filter : bool
        a boolean variable determining whether a local filter is used
        by default, it is false and will be later determined by PDAF function
        in init_pdaf
    sv : state_vector.state_vector
        a reference to the state vector object

    Methods
    -------
    init_n_domains(self, n_domains:int)
        initialise the number of local domains
    """

    def __init__(self, nemo_grid:model.NemoDomain, fields:dict[str, sfields.StateField]) -> None:
        """class constructor

        Parameters
        ----------
        wet_grid : model.WetGrid
            a reference to the wet grid object
        """
        self.nemo_grid = nemo_grid
        self.fields = fields
        self.local_filter = False

    def set_domain_limits(self):
        """Set the local domain limits for the PDAFomi used for local filters
        with use_global=0."""
        lim_coords = np.array(
            [
                [np.deg2rad(np.nanmin(self.nemo_grid.coords.glamt)),
                    np.deg2rad(np.nanmax(self.nemo_grid.coords.glamt))],
                [np.deg2rad(np.nanmax(self.nemo_grid.coords.gphit)),
                    np.deg2rad(np.nanmin(self.nemo_grid.coords.gphit))],
            ], order='F'
        )
        pyPDAF.PDAFomi.set_domain_limits(lim_coords)

    def init_n_domains_pdaf(self, _step:int, _ndomains:int) -> int:
        """Initialise the number of local analysis domains.

        Parameters
        ----------
        step : int
            current time step
        ndomains: int
            number of local domains on current processor

        Returns
        -------
        n_domains_p : int
            PE-local number of analysis domains
        """
        return self.nemo_grid.wet_grid.nwet2d

    def init_dim_l_pdaf(self, _step:int, domain_p:int, dim_l:int) -> int:
        """Initialise the local state dimension and local-to-global indices.

        The function returns
        the dimension of local state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            index of current local analysis domain
        dim_l : int
            dimension of local state vector

        Returns
        -------
        dim_l : int
            dimension of local state vector
        """
        # initialize local state dimension
        wet_pts = self.nemo_grid.wet_grid.wet_pts
        point = domain_p - 1
        n_layers = wet_pts[2, point]
        col_start = wet_pts[4, point]

        dim_l = sum(1 if field.ndims == 2 else n_layers for field in self.fields.values())
        id_lstate_in_pstate = np.empty(dim_l, dtype=np.intc)
        cnt = 0
        for field in self.fields.values():
            if field.ndims == 2:
                id_lstate_in_pstate[cnt] = field.off + domain_p
                cnt += 1
            else:
                ids = field.off + col_start + np.arange(1, n_layers + 1, dtype=np.intc)
                id_lstate_in_pstate[cnt:cnt + n_layers] = ids
                cnt += n_layers
        pyPDAF.PDAFlocal.set_indices(dim_l, id_lstate_in_pstate)
        return dim_l
