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

import config
import parallelisation
import state_vector


class localisation:

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

    def __init__(self, sv:state_vector.state_vector) -> None:
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
        self.sv:state_vector.state_vector = sv

    def init_n_domains_pdaf(self, step:int, ndomains:int) -> int:
        """initialize the number of local domains

        Parameters
        ----------`
        step : int
            current time step
        ndomains: int
            number of local domains on current processor

        Returns
        -------
        n_domains_p : int
            PE-local number of analysis domains
        """
        log.logger.info(f'ndomains: ndomains {self.sv.dim_state_p}')
        return self.sv.dim_state_p

    def set_lim_coords(self, nx_p:int, ny_p:int, pe:parallelisation.parallelisation) -> None:
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

        lim_coords = np.zeros((2, 2))
        lim_coords[0, 0] = float(off_nx + 1)
        lim_coords[0, 1] = float(off_nx + nx_p)
        lim_coords[1, 0] = ny_p
        lim_coords[1, 1] = 1

        PDAF.omi_set_domain_limits(lim_coords)

    def init_dim_l_pdaf(self, step:int, domain_p:int, dim_l:int) -> int:
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
        return dim_l

    def g2l_state_pdaf(self, step:int, domain_p:int, dim_p:int, state_p:np.ndarray, dim_l:int, state_l:np.ndarray) -> np.ndarray:
        """convert state vector (on each processor) to domain local state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        dim_p : int
            dimension of state vector
        state_p : ndarray
            state vector
        dim_l : int
            dimension of domain local state vector
        state_l : ndarray
            domain local state vector for local analysis
        """
        # generic initialization
        # domain_p is the *domain_p*-th local domain on the local processor
        # The dimension of the state vector on the local domain is dim_l defined in init_dim_l_pdaf
        # Here, because we set dim_l = 1, each local domain is one element of the state vector
        # In more complex cases, it is possible to define a relationship matrix between local domain
        # and the element of the state vector.
        state_l[:] = state_p[domain_p - 1]
        return state_l

    def l2g_state_pdaf(self, step:int, domain_p:int, dim_l:int, state_l:np.ndarray, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """convert local state vector to PE-local global state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        dim_l : int
            dimension of local state vector
        state_l : ndarray
            local state vector for local analysis
        dim_p : int
            dimension of state vector
        state_p : ndarray
            state vector
        """
        state_p[domain_p -1] = state_l[:]
        return state_p
