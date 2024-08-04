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
import pyPDAF.PDAF as PDAF


class Localization:

    """class for localization information and user-supplied functions

    Attributes
    ----------
    coords_l : ndarray
        coordinates of local analysis domain
    id_lstate_in_pstate : ndarray
        indices of local state vector in PE-local global state vector
    loc_weight : int
        - (0) constant weight of 1
        - (1) exponentially decreasing with SRADIUS
        - (2) use 5th-order polynomial
        - (3) regulated localization of R with mean error variance
        - (4) regulated localization of R with single-point error variance
    cradius : float
        cut-off radius for local observation domain
    sradius : float
        support radius for 5th order polynomial
        or radius for 1/e for exponential weighting
    """

    def __init__(self, loc_weight, cradius, sradius):
        """class constructor

        Parameters
        ----------
        loc_weight : int
            - (0) constant weight of 1
            - (1) exponentially decreasing with SRADIUS
            - (2) use 5th-order polynomial
            - (3) regulated localization of R with mean error variance
            - (4) regulated localization of R with single-point error variance
        cradius : float
            cut-off radius for local observation domain
        sradius : float
            support radius for 5th order polynomial
            or radius for 1/e for exponential weighting
        """
        self.loc_weight = loc_weight
        self.cradius = cradius
        self.sradius = sradius

    def set_lim_coords(self, dims_p, pe):
        """set local domain

        Parameters
        ----------
        dims_p : int
            size of PE-local state vector
        pe : 'parallelization.parallelization'
            parallelization object
        """
        # Get offset of local domain in global domain in x-direction
        off_nx = dims_p[-1]*pe.mype_filter

        lim_coords = np.zeros((2, 2))
        lim_coords[0, 0] = float(off_nx + 1)
        lim_coords[0, 1] = float(off_nx + dims_p[-1])
        lim_coords[1] = 0

        PDAF.omi_set_domain_limits(lim_coords)

    def init_n_domains_pdaf(self, assim_opt, step, ndomains):
        """get the number of analysis domains

        Parameters
        ----------
        assim_opt : `AssimilationOptions.AssimilationOptions`
            assim_opt object
        step : int
            current time step

        Returns
        -------
        n_domains_p : int
            PE-local number of analysis domains
        """
        return assim_opt.dim_state_p

    def init_dim_l_pdaf(self, dims_p, mype_filter, step, domain_p, dim_l):
        """initialise the local dimension of PDAF.

        The function set coordinates of local analysis domain
        and its conversion to global state vector. It also returns
        the dimension of local state vector

        Parameters
        ----------
        dims_p : int
            size of PE-local state vector
        mype_filter : int
            rank of the process in filter communicator
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
        # initialize coordinates of local domain
        # we use grid point indices as coordinates,
        #  but could e.g. use meters
        self.coords_l = np.zeros(2)
        offset = mype_filter*np.prod(dims_p)
        self.coords_l[0] = domain_p + offset
        self.coords_l[0] = np.ceil(self.coords_l[0]/dims_p[0])
        self.coords_l[1] = domain_p + offset
        self.coords_l[1] = self.coords_l[1] \
            - (self.coords_l[0] - 1)*dims_p[0]
        # initialize array of indices of the local state
        #  vector elements in the global state vector.

        # allocate array
        self.id_lstate_in_pstate = np.zeros(dim_l, dtype=int)

        # here the local domain is a single grid point
        # and variable given by domain_p
        #print (domain_p)
        self.id_lstate_in_pstate[0] = domain_p

        return dim_l

    def g2l_state_pdaf(self, step, domain_p, dim_p, state_p, dim_l, state_l):
        """convert PE-local global state vector to local state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        state_p : ndarray
            PE-local global state vector
        state_l : ndarray
            local state vector for local analysis
        """
        # generic initialization
        # using id_lstate_in_pstate set in init_dim_l_pdaf
        state_l = state_p[self.id_lstate_in_pstate -1]
        return state_l

    def l2g_state_pdaf(self, step, domain_p, dim_l, state_l, dim_p, state_p):
        """convert local state vector to PE-local global state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        state_l : ndarray
            local state vector for local analysis
        state_p : ndarray
            PE-local global state vector
        """
        state_p[self.id_lstate_in_pstate -1] = state_l
        return state_p
