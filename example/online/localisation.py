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

import config
import log
import parallelisation
import state_vector


class Localisation:
    """class for localization information and user-supplied functions

    Attributes
    ----------
    loc_weight : int
        - (0) constant weight of 1
        - (1) exponentially decreasing with sradius
        - (2) use 5th-order polynomial
        - (3) regulated localization of R with mean error variance
        - (4) regulated localization of R with single-point error variance
    cradius : float
        range for local observation domain
    sradius : float
        support range for 5th order polynomial
        or radius for 1/e for exponential weighting
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

    def __init__(self, sv:state_vector.StateVector) -> None:
        """class constructor

        Parameters
        ----------
        sv : 'state_vector.state_vector'
            state vector object
        """
        self.loc_weight:int = config.loc_weight
        self.cradius:float = config.cradius
        self.sradius:float = config.sradius
        # a boolean variable determining whether a local filter is used
        # by default, it is false and will be later determined by PDAF function
        # in init_pdaf
        self.local_filter:bool = False
        # a reference to the state vector object
        self.sv:state_vector.StateVector = sv

    def init_n_domains_pdaf(self, _step:int, _ndomains:int) -> int:
        """initialize the number of local domains

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
        output_str = f'ndomains: ndomains {self.sv.dim_state_p}'
        log.logger.info(output_str)
        return self.sv.dim_state_p

    def set_lim_coords(self, nx_p:int, ny_p:int, pe:parallelisation.Parallelisation) -> None:
        """set local domain

        Parameters
        ----------
        nx_p : int
            size of PE-local state vector in x-direction
        ny_p : int
            size of PE-local state vector in y-direction
        pe : 'parallelisation.parallelisation'
            parallelization object
        """
        # Get offset of local domain in global domain in x-direction
        # The the following link for more information
        # https://pdaf.awi.de/trac/wiki/PDAFomi_additional_functionality#PDAFomi_set_domain_limits
        # note that this is Fortran-based documentation where array index
        # starts from 1
        off_nx = nx_p*pe.mype_filter

        lim_coords = np.zeros((2, 2), order='F')
        lim_coords[0, 0] = float(off_nx + 1)
        lim_coords[0, 1] = float(off_nx + nx_p)
        lim_coords[1, 0] = ny_p
        lim_coords[1, 1] = 1

        pyPDAF.PDAFomi.set_domain_limits(lim_coords)

    def init_dim_l_pdaf(self, _step:int, domain_p:int, dim_l:int) -> int:
        """initialise the local dimension of PDAF.

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
        dim_l = 1
        id_lstate_in_pstate:np.ndarray = domain_p*np.ones((dim_l), dtype=np.intc)
        pyPDAF.PDAFlocal.set_indices(dim_l, id_lstate_in_pstate)
        return dim_l
