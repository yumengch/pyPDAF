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


# verbosity of the PDAF screen output
screen:int = 2

### Filepath ###
# path to initial ensemble, the filename is formatted for different time step
# here a relative path is given
init_ens_path:str = os.path.join('inputs_online', 'ens_{i}.txt')
# path to initial truth
init_truth_path:str = os.path.join('inputs_online', 'true_initial.txt')


#### Model configurations ####
# number of grid points in x-direction
nx = 36
# number of grid points in y-direction
ny = 18
# initial time step
init_step = 0
# total number of time steps
nsteps:int = 18

#### filter options ####
# number of ensemble members
dim_ens:int = 4
# number of parallel tasks run simultanously
n_modeltasks:int = 4
# the type of filter used
# 1=SEIK, 2=EnKF, 3=LSEIK, 4=ETKF, 5=LETKF, 6=ESTKF, 7=LESTKF
# 8=LEnKF, 9=NETF, 10=LNETF, 11=LKNETF, 12=PF, 100=GENOBS,
# 200=3DVar, 0=SEEK
# For a simplified documentation, see:https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF
# Different DA scheme requires different user-supplied functions
# More information can be found in the PDAF documentation
filtertype:int = 6
# Variants of each DA scheme check https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF
subtype:int = 0
# type of forgetting factor
# - (0) fixed
# - (1) global adaptive
# - (2) local adaptive for LSEIK/LETKF/LESTKF
type_forget:int = 0
# forgetting factor
forget:float = 1.0

#### transformation-related options for Kalman filter (KF) #####
# type of ensemble transformation
type_trans:int = 0
# Ways of doing square-root of thetransform matrix
# (0) symmetric square root, (1) Cholesky decomposition
type_sqrt:int = 0
# (1) to perform incremental updating (only in SEIK/LSEIK!)
incremental:int = 0
# Definition of factor in covar. matrix used in SEIK
# - (0) for dim_ens^-1 (old SEIK)
# - (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
# This parameter has also to be set internally in PDAF_init.
covartype:int = 1
#  rank to be considered for inversion of HPH in analysis of EnKF
#     (0) for analysis w/o eigendecomposition
#     if set to >=ensemble size, it is reset to ensemble size - 1
rank_analysis_enkf:int = 0


#### localisation options ####
# Type of localization function (0: constant, 1: exponential decay, 2: 5th order polynomial)
loc_weight:int = 3
# localization cut-off radius in grid points
cradius:float = 6.0
# Support radius for localization function
sradius:float = cradius
