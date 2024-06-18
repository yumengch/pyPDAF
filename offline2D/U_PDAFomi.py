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


def init_dim_obs_pdafomi(list_of_obs, cradius,
                         mype_filter, nx, nx_p, step, dim_obs):
    """initialise observation dimensions

    Parameters
    ----------
    list_of_obs : list
        a list of all types of observations
    cradius : float
        cut-off radius for local observation domain
    mype_filter : int
        rank of the PE in filter communicator
    nx : ndarray
        integer array for grid size
    nx_p : ndarray
        integer array for PE-local grid size
    step : int
        current time step
    dim_obs : int
        dimension of observation vector

    Returns
    -------
    TYPE
        Description
    """
    d_obs = 0
    for obs in list_of_obs:
        if(obs.doassim):
            obs.init_dim_obs(step, dim_obs, cradius,
                             mype_filter, nx, nx_p)
            d_obs += obs.dim_obs
    dim_obs = d_obs
    return dim_obs


def obs_op_pdafomi(list_of_obs, step, dim_p, dim_obs_p, state_p, ostate):
    """turn state vector to observation space

    Parameters
    ----------
    list_of_obs : list
        a list of all types of observations
    step : int
        current time step
    state_p : ndarray
        local PE state vector
    ostate : ndarray
        state vector in obs space
    """

    for obs in list_of_obs:
        ostate = obs.obs_op(step, state_p, ostate)
    return ostate


def init_dim_obs_l_pdafomi(list_of_obs, localization,
                           domain_p, step, dim_obs, dim_obs_l):
    """initialise local observation dimension

    Parameters
    ----------
    list_of_obs : list
        a list of all types of observations
    localization : `Localization.Localization`
        a localization info obejct
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    dim_obs_l : int
        dimension of local observation vector
    """

    dim_obs_lA = 0
    for obs in list_of_obs:
        dim_obs_lB=obs.init_dim_obs_l(localization,
                           domain_p, step, dim_obs, dim_obs_lA)
        dim_obs_l = dim_obs_l + dim_obs_lB
    return dim_obs_l


def localize_covar_pdafomi(list_of_obs, localization,
                           mype_filter, nx_p, HP_p, HPH):
    """localize covariance matrix

    Parameters
    ----------
    list_of_obs : list
        a list of all types of observations
    localization : `Localization.Localization`
        a localization info obejct
    mype_filter : int
        rank of the PE in filter communicator
    nx_p : ndarray
        integer array for PE-local grid size
    HP_p : ndarray
        matrix HP
    HPH : ndarray
        matrix HPH
    """
    dim_p = HPH.shape[0]
    coords_p = np.zeros((2, dim_p))
    offset = mype_filter*nx_p

    coords_p[0] = np.where(np.ones(nx_p))[1] + offset
    coords_p[1] = np.where(np.ones(nx_p))[0]

    for i_obs_f in list_of_obs:
        i_obs_f.localize_covar(localization, HP_p, HPH, coords_p)
