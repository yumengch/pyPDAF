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
from parallelization import parallelization
import Model
import OBS_A
import OBS_B
import pyPDAF.PDAF as PDAF
import PDAF_caller


from Localization import Localization
from AssimilationDimensions import AssimilationDimensions
from FilterOptions import FilterOptions


def main():
    USE_PDAF = True

        ###### Model configuration ######################
    # grid dimensions
    nx = 36
    ny = 18

    # Number of time steps to compute
    nt = 18
    
    ###### PDAF Configuration #######################

    # number of ensemble members
    dim_ens = 4
    # using error space transform Kalman filter (4=ETKF, 5=LETKF)
    filtertype = 4
    # standard form
    subtype = 0
    # forgetting factor
    forget = 0.1

    # time interval between observations
    dtobs = 2
    # Observation error standard deviation
    rms_obs = 0.5
    # Assimilate observation type A
    assim_A = True
    # Assimilate observation type B
    assim_B = True
    
    # Type of localization function (0: constant, 1: exponential decay, 2: 5th order polynomial)
    loc_weight = 0
    # localization cut-off radius in grid points
    cradius = 0
    # Support radius for localization function
    sradius = cradius
    
    ###############################


    pe = parallelization(dim_ens=dim_ens, n_modeltasks=dim_ens, screen=2)

    # Initial Screen output
    if (pe.mype_world == 0):
        print('     +++++ PDAF online mode +++++')
        print('     2D model with parallelization')
        print('')
        if pe.npes_world==1:
            print('Running on ', pe.npes_world, 'PE')
        else:
            print('Running on ', pe.npes_world, 'PEs')

    # Initialize model
    model = Model.Model((ny, nx), nt=nt, pe=pe)
    if not USE_PDAF:
        model.init_field('inputs_online/true_initial.txt',
                          pe.mype_model)
    
    # Initialize observations
    obs = []
    obscnt = 0

    if assim_A:
        obscnt += 1
        obs.append(OBS_A.OBS_A('A', obscnt, pe.mype_filter, pe.task_id,
                           model.dims, 1, dtobs, rms_obs))
    if assim_B:
        obscnt += 1
        obs.append(OBS_B.OBS_B('B', obscnt, pe.mype_filter, pe.task_id,
                               model.dims, 1, dtobs, rms_obs))

    PDAF.omi_init(len(obs))


    # Set state dimension and ensemble size for PDAF
    assim_dim = AssimilationDimensions(model=model,
                                       dim_ens=dim_ens,
                                       experiment=f'out_N{dim_ens}_f{forget}')
    # Set options for PDAF
    filter_options = FilterOptions(filtertype=filtertype,
                                   subtype=subtype,
                                   forget=forget)
    # Set localization parameters
    localization = Localization(loc_weight=loc_weight,
                                cradius=cradius,
                                sradius=sradius)

    # init PDAF
    PDAF_caller.init_pdaf(assim_dim,
                          filter_options,
                          localization,
                          model, obs, pe,
                          screen=2)

    # Run model and assimilate
    for it in range(nt):
        model.step(pe, it, USE_PDAF)
        if USE_PDAF:
            PDAF_caller.assimilate_pdaf(assim_dim,
                                        filter_options,
                                        localization,
                                        model, obs, pe)

    pe.finalize_parallel()


if __name__ == '__main__':
    main()
