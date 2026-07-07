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
import configparser
import log

import numpy as np
import pyPDAF.PDAFomi

import model
import sfields

class ObsSST:
    """PDAF-OMI observation functions for synthetic SST-like observations.

    The current implementation does not read an external SST product. It
    creates synthetic observations by taking the background ensemble mean at
    wet surface grid points and adding Gaussian noise with standard deviation
    ``noise_amp``. This makes the class useful as a template for real
    satellite/in-situ observations: replace observation generation in
    :meth:`init_dim` while keeping the PDAF-OMI registration calls.

    Attributes
    ----------
    i_obs : int
        index of the observation type
    obs_name : str
        name of the current observation type
    filename : str
        observation filename to be read
    dim_obs_p : int
        dimension size of the PE-local observation vector
    disttype : int
        type of distance computation to use for localization
    doassim : int
        whether to assimilate this observation type
    ncoord : int
        number of coordinate dimension
    nrows : int
        number of rows in ocoord_p
    cradius : float
        range for local observation domain
    sradius : float
        support range for 5th order polynomial
        or radius for 1/e for exponential weighting
    pe : parallelisation.parallelisation
        parallelization object
    local : localisation.localisation
        localisation object
    """

    def __init__(self, i_obs:int, config:configparser.SectionProxy) -> None:
        # i_obs-th observations in the system starting from 1
        self.i_obs = i_obs
        assert self.i_obs >= 1, 'observation count must start from 1'
        # Identifier used in logs and by the observation factory.
        self.obs_name = 'sst'
        # Name of the state-vector section observed by this module.
        self.field_name = config.get('field_name', 'sst')
        # Placeholder for future real-observation readers.
        self.filename = config.get('obs_path', 'filename')
        # Assimilation switch: 1=active, 0=registered but skipped.
        self.doassim = config.getint('doassim', 1)
        # Type of distance computation to use for localization
        # 3 is for the great-circle distance on a sphere
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype
        self.disttype = 3
        # Number of coordinates use for distance computation
        # Here, this observation only uses a 2-dimensional surface domain
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsncoord
        self.ncoord = 2
        # nrows depends on the necessity of interpolating
        # observations onto model grid
        # if nrows = 1, observations are on the model grid points
        # Here, we use 1 for synthetic observations.
        # when interpolation is required,
        # this is the number of grid points required for interpolation.
        # For example, nrows = 4 for bi-linear interpolation in 2D,
        # and 8 for 3D linear interpolation.
        # https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsid_obs_p
        # More information about interpolation is available at
        # https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin
        self.nrows = 1
        # Synthetic observation standard deviation.
        self.rms = config.getfloat('noise_amp', 0.5)
        # Localisation weighting: 0=constant, 1=exponential,
        # 2=fifth-order polynomial, 3/4=regulated PDAF variants.
        self.loc_weight = config.getint('loc_weight', 2)
        # support range for 5th order polynomial or radius for 1/e for exponential weighting
        self.sradius = config.getfloat('sradius', 1.)
        # radius for local observation domain
        self.cradius = config.getfloat('cradius', 1.)


    def init_dim(self, dim_obs:int, field:sfields.StateField, model_grid:model.NemoDomain) -> int:
        """Create/register SST observations for PDAF-OMI.

        This method supplies observation values, inverse variances, observed
        state indices, coordinates, distance metric, and localisation settings.
        For real observations, replace the synthetic ``obs_p`` construction
        with file/database input and interpolation metadata.

        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension of observation vector
        field : sfields.StateField
            state field object
        model_grid : model.NemoDomain
            model grid object

        Returns
        -------
        int
            dimension of observation vector
        """
        output_str = f'Assimilate observations: {self.obs_name}'
        log.info(output_str)

        # switch for assimilation of the observation
        pyPDAF.PDAFomi.set_doassim(self.i_obs, self.doassim)
        # Type of distance computation to use for localization
        # It is mandatory for OMI even if we don't use localisation
        pyPDAF.PDAFomi.set_disttype(self.i_obs, self.disttype)
        # Number of coordinates used for distance computation
        pyPDAF.PDAFomi.set_ncoord(self.i_obs, self.ncoord)

        dim_obs_p = model_grid.wet_grid.nwet2d

        if dim_obs_p == 0:
            obs_p = np.zeros((1))
            ocoord_p = np.zeros((self.ncoord, 1))
            ivar_obs_p = np.zeros((1))
            id_obs_p = np.zeros((self.nrows, 1), dtype=np.intc)
            pyPDAF.PDAFomi.set_id_obs_p(self.i_obs, self.nrows, 1, id_obs_p)
        else:
            # obtain the observation vector
            rng = np.random.default_rng(seed=12345)
            obs_p = model_grid.bkg_ens.mean(1)[field.off : field.off + field.dim]
            obs_p += rng.normal(scale=self.rms, size=dim_obs_p)
            # inverse of observation variance
            # here we specify/hard-code the standard deviation of observation is 0.5
            ivar_obs_p = (1./self.rms/self.rms)*np.ones_like(obs_p)
            # coordinate of each observations
            j = model_grid.wet_grid.wet_pts[5]
            i = model_grid.wet_grid.wet_pts[6]
            ocoord_p = np.zeros((self.ncoord, dim_obs_p), order='F')
            ocoord_p[0] = np.deg2rad(model_grid.coords.nav_lon[j, i])
            ocoord_p[1] = np.deg2rad(model_grid.coords.nav_lat[j, i])

            id_obs_p = np.zeros((self.nrows, len(obs_p)), dtype=np.intc, order='F')
            # pyPDAF.PDAFomi.set_debug_flag(10)
            # pyPDAF.PDAF.set_debug_flag(10)
            # The relationship between observation and state vector
            # id_obs_p gives the indices of observed field in state vector
            # the index starts from 1
            # The index is based on the full state vector instead of the index in the local domain
            # The following code is used because in our example, observations are masked and
            # have the same shape as the model grid
            id_obs_p[0] = np.arange(field.off + 1, field.off + dim_obs_p + 1, dtype=np.intc)
            pyPDAF.PDAFomi.set_id_obs_p(self.i_obs, self.nrows, dim_obs_p, id_obs_p)

        # PDAF need to gather observation information
        dim_obs = pyPDAF.PDAFomi.gather_obs(self.i_obs, dim_obs_p,
                                            obs_p,
                                            ivar_obs_p,
                                            ocoord_p, self.ncoord,
                                            self.cradius)
        # set covariance localisation for stochastic EnKF/EAKF/EnSRF
        if pyPDAF.PDAF.get_local_type() == 2:
            dim_p = model_grid.wet_grid.nwet3d if field.ndims == 3 else model_grid.wet_grid.nwet2d
            coords_p = np.zeros((2, dim_p), order='F')
            coords_p[0] = np.deg2rad(model_grid.coords.nav_lon[j, i])
            coords_p[1] = np.deg2rad(model_grid.coords.nav_lat[j, i])
            pyPDAF.PDAFomi.set_localize_covar_iso(self.i_obs, dim_p,
                                                  self.ncoord, coords_p,
                                                  self.loc_weight,
                                                  self.cradius,
                                                  self.sradius)

        return dim_obs

    def init_dim_obs_l(self, coords_l:np.ndarray, dim_obs_l:int) -> int:
        """Initialise local observation vector size for domain localisation.

        Parameters
        ----------
        domain_p : int
            index of current local analysis domain
        coords_l : ndarray
            coordinates of local domain
        dim_obs_l : int
            dimension of local observation vector

        Returns
        -------
        dim_obs_l : int
            dimension of local observation vector
        """
        # initialize coordinates of local domain

        return pyPDAF.PDAFomi.init_dim_obs_l_iso(self.i_obs, coords_l, self.loc_weight,
                                                 self.cradius, self.sradius, dim_obs_l)

    def obs_op(self, _step:int, state_p:np.ndarray, ostate:np.ndarray) -> np.ndarray:
        """Apply the SST observation operator.

        Synthetic SST observes grid-point state entries directly, so this uses
        PDAF-OMI's grid-point operator. Replace this method when an observation
        requires interpolation or a nonlinear operator.

        Parameters
        ----------
        step : int
            current time step
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix

        Returns
        -------
        ostate : ndarray
            state vector transformed by identity matrix
        """
        return pyPDAF.PDAFomi.obs_op_gridpoint(self.i_obs, state_p, ostate)
