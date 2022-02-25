"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and National Centre for Earth Observation

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

This module calls PDAFomi iso_C_binding wrappers.

Function calls convert Python objects for Fortran routines.
Detailed information on the PDAFomi properties or subroutines can be
found at https://pdaf.awi.de/trac/wiki/OverviewOfOMIRoutines and 
http://pdaf.awi.de/trac/wiki/ImplementFilterAnalysisOverview 
"""
import numpy as np


def init(int n_obs):
    """Python wrapper for PDAFomi init.

        Parameters
        ----------
        n_obs : int
            number of observation types
    """
    c__init_pdafomi(&n_obs)


def setOMIessential(int i_obs, int doassim, 
                    int disttype, int ncoord, id_obs_p):
    """Python wrapper for setting PDAFomi objects.

        Parameters
        ----------
        i_obs : int
            index of observation types
        doassim : int
            whether the observation is assimilated
        disttype : int
            type of distance measure for localisation
        ncoord : int
            number of coordinate dimension
        id_obs_p : ndarray
            indices of process-local observed field in state vector
    """
    c__set_pdafomi_doassim(&i_obs, &doassim)
    c__set_pdafomi_disttype(&i_obs, &disttype)
    c__set_pdafomi_ncoord(&i_obs, &ncoord)
    shape = id_obs_p.shape
    cdef int nrows = shape[0]
    cdef int dim_obs_p = shape[1]
    cdef int[::1, ::] id_obs_p_view = np.array(id_obs_p, order='F', dtype=np.intc)
    c__set_pdafomi_id_obs_p(&i_obs, &nrows, &dim_obs_p, &id_obs_p_view[0][0])


def set_icoeff_p(int i_obs, icoeff_p):
    """Python wrapper for setting PDAFomi icoeff_p.

        Parameters
        ----------
        i_obs : int
            index of observation types
        icoeff_p : ndarray
            2d array for interpolation coefficients for obs. operator
    """
    cdef int nrows = icoeff_p.shape[0]
    cdef int dim_obs_p = icoeff_p.shape[1]
    cdef double[::1, ::] icoeff_p_view = np.array(icoeff_p, order='F')
    c__set_pdafomi_icoeff_p(&i_obs, &nrows, &dim_obs_p, &icoeff_p_view[0][0])


def set_domainsize(int i_obs, domainsize):
    """Python wrapper for setting PDAFomi domainsize.

        Parameters
        ----------
        i_obs : int
            index of observation types
        domainsize : ndarray
            1d array for size of domain for periodicity (<=0 for no periodicity)
    """
    cdef double[::1] domainsize_view = domainsize
    cdef int ncoord = len(domainsize)
    c__set_pdafomi_domainsize(&i_obs, &ncoord, &domainsize_view[0])


def set_obs_err_type(int i_obs, int obs_err_type):
    """Python wrapper for setting PDAFomi obs_err_type.

        Parameters
        ----------
        i_obs : int
            index of observation types
        obs_err_type : int
            type of observation error
    """
    c__set_pdafomi_obs_err_type(&i_obs, &obs_err_type)


def set_use_global_obs(int i_obs, int use_global_obs):
    """Python wrapper for setting PDAFomi use_global_obs.

        Parameters
        ----------
        i_obs : int
            index of observation types
        use_global_obs : int
           Whether to use (1) global full obs. or 
           (0) obs. restricted to those relevant for a process domain
    """
    c__set_pdafomi_use_global_obs(&i_obs, &use_global_obs)


def gather_obs(int i_obs, int dim_obs_p, int nrows, obs_p, ivar_obs_p, 
                ocoord_p, double local_range):
    """Python wrapper for calling PDAFomi_gather_obs.

        Parameters
        ----------
        i_obs : int
            index of observation types
        dim_obs_p : int
            PE_local observation dimension
        nrows : int
            number of rows in ocoord_p
        obs_p : ndarray
            vector of process-local observations
        ivar_obs_p : ndarray
            vector of process-local inverse observation error variance
        ocoord_p : ndarray
            2d array of process-local observation coordinates
        local_range : double
            lcalization radius (the maximum radius used in this process domain) 

        Returns
        -------
        dim_obs : int
            dimension of the entire observation vector
    """
    cdef double[::1, ::] ocoord_p_view = np.array(ocoord_p, order='F')
    cdef double[::1] obs_p_view = obs_p
    cdef double[::1] ivar_obs_p_view = ivar_obs_p
    cdef int dim_obs
    c__pdafomi_gather_obs(&i_obs, &nrows, &dim_obs_p,
                           &obs_p_view[0], &ivar_obs_p_view[0], 
                           &ocoord_p_view[0][0], &local_range, &dim_obs)
    return dim_obs


def set_domain_limits(lim_coords):
    """Python wrapper for calling PDAFomi_set_domain_limits.

        Parameters
        ----------
        lim_coords : ndarray
            sets the limiting coordinates of a process domain
            
    """
    cdef double[::1, ::] lim_coords_view = np.array(
                                            lim_coords, order='F')
    c__pdafomi_set_domain_limits(&lim_coords_view[0][0])


def obs_op_gridpoint(int i_obs, state_p, ostate):
    """Python wrapper for calling PDAFomi_obs_op_gridpoint.

        Parameters
        ----------
        i_obs : int
            index of observation types
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix
    """
    cdef int dim_p, dim_obs
    dim_p = len(state_p) 
    dim_obs = len(ostate)
    cdef double[::1] state_p_view = state_p
    cdef double[::1] ostate_view = ostate
    c__pdafomi_obs_op_gridpoint(&i_obs, &dim_p, &dim_obs, 
                                &state_p_view[0], &ostate_view[0])


def init_dim_obs_l(int i_obs, coords_l, int loc_weight, 
                    double local_range, double srange):
    """Python wrapper for calling PDAFomi_init_dim_obs_l.

        Parameters
        ----------
        i_obs : int
            index of observation types
        coords_l : ndarray
            coordinates of local domain
        loc_weight : int
            type of localizing weighting of observations
        local_range : float
            range for local observation domain
        srange : float
            support range for 5th order polynomial
            or radius for 1/e for exponential weighting

        Returns
        -------
        dim_obs_l : int
            dimension of local observations
            
    """
    cdef double[::1, ::] coords_l_view = np.array(coords_l, order='F')
    cdef int dim_obs_l
    c__pdafomi_init_dim_obs_l(&i_obs, &coords_l_view[0][0], &loc_weight, 
                              &local_range, &srange, &dim_obs_l)
    return dim_obs_l


def localize_covar(int i_obs, int loc_weight, 
                   double local_range, double srange, 
                   coords_p, hp_p, hph):
    """Python wrapper for calling PDAFomi_localize_covar.

        Parameters
        ----------
        i_obs : int
            index of observation types
        loc_weight : int
            type of localizing weighting of observations
        local_range : float
            range for local observation domain
        srange : float
            support range for 5th order polynomial
            or radius for 1/e for exponential weighting
        coords_p : ndarray
            coordinates of state vector elements
        hp_p : ndarray
            matrix HPH
        hph : ndarray
            PE local part of matrix HP
    """
    cdef double[::1, ::] coords_p_view = np.array(coords_p, order='F')
    cdef int dim_coords, dim_p, dim_obs
    dim_coords, dim_p = coords_p.shape
    dim_obs = hp_p.shape[0]
    cdef double[::1, ::] hp_p_view = np.array(hp_p, order='F')
    cdef double[::1, ::] hph_view = np.array(hph, order='F')
    c__pdafomi_localize_covar(&i_obs, &dim_p, &dim_obs, &dim_coords,
                              &loc_weight, &local_range, &srange, 
                              &coords_p_view[0][0], 
                              &hp_p_view[0][0], 
                              &hph_view[0][0]);


def deallocate_obs(int i_obs, int step):
    """Python wrapper for calling PDAFomi_deallocate_obs.

        Parameters
        ----------
        i_obs : int
            index of observation types
        step : int
            current time step
    """
    c__pdafomi_deallocate_obs(&i_obs, &step);
