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

This is a template file for PDAF configurations.

The module can be modified to be adjusted to read configparser/YAML/JSON etc. to initialise
options in PDAF and the DA system.
"""
import os


# name of the observation
obs_name = 'A'
# path to the observation files
obs_path:str = os.path.join('/home/runner/work/pyPDAF/pyPDAF/inputs_online', 'obs_step{i}.txt')
# time steps between observations / assimilation frequency
dtobs:int = 2
# Observation error standard deviation
rms_obs:float = 0.5
# Switch for assimilating observation type A
assim:bool = True
# Type of distance computation to use for localization
# It is mandatory for OMI even if we don't use localisation
# 0=Cartesian 1=Cartesian periodic
# details of this property can be seen at
# https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsdisttype
# currently, only PDAF V2.1 is supported
disttype:int = 0
# Number of coordinates use for distance computation
# Here, it is a 2-dimensional domain
# https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsncoord
ncoord:int = 2
# nrows depends on the necessity of interpolating
# observations onto model grid
# if nrows = 1, observations are on the model grid points
# when interpolation is required, 
# this is the number of grid points required for interpolation.
# For example, nrows = 4 for bi-linear interpolation in 2D,
# and nrows = 8 for 3D linear interpolation.
# https://pdaf.awi.de/trac/wiki/OMI_observation_modules#thisobsid_obs_p
# More information about interpolation is available at
# https://pdaf.awi.de/trac/wiki/OMI_observation_operators#PDAFomi_get_interp_coeff_lin
nrows:int = 1
# missing value in the observation
missing_value:float = -999.
