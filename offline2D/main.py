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
import OBS
import pyPDAF.PDAF as PDAF
import PDAF_caller


from Localization import Localization
from AssimilationDimensions import AssimilationDimensions
from FilterOptions import FilterOptions
from Inflation import Inflation


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
    subtype = 5
    # forgetting factor
    forget = 1.0

    # time interval between observations
    dtobs = 1
    # Observation error standard deviation
    rms_obs = 0.5
    
    # Type of localization function (0: constant, 1: exponential decay, 2: 5th order polynomial)
    loc_weight = 0
    # localization cut-off radius in grid points
    local_range = 0
    # Support radius for localization function
    srange = local_range
    
    ###############################

    # dimension of the state vector
    # if model is parallelised, this is the dimension of state vector on each process
    dim_state_p = nx*ny

    if USE_PDAF:
        pe = parallelization(dim_ens=dim_ens, n_modeltasks=1, screen=1)

    # Initial Screen output
    if (pe.mype_world == 0):
        print('     +++++ PDAF offline mode +++++')
        print('      Data assimilation with PDAF')
        print('')
        if pe.npes_world==1:
            print('Running on ', pe.npes_world, 'PE')
        else:
            print('Running on ', pe.npes_world, 'PEs')

    # Initialize model
    model = Model.Model((ny, nx), nt=nt, pe=pe)
    
    obs = []
    for typename in ['A',]:
        obs.append(OBS.OBS(typename, pe.mype_filter, pe.task_id,
                           model.nx, 1, dtobs, rms_obs))

    assim_dim = AssimilationDimensions(model, dim_ens=dim_ens)
    filter_options = FilterOptions(filtertype=filtertype, subtype=subtype)
    filter_options.setTransformTypes(type_trans=0, type_sqrt=0,
                                     incremental=0, covartype=1,
                                     rank_analysis_enkf=0)
    infl = Inflation(type_forget=0, forget=forget)
    localization = Localization(loc_weight=loc_weight, local_range=local_range,
                                srange=srange)


    # init observations
    PDAF.omi_init(len(obs))

    # init PDAF
    PDAF_caller.init_pdaf(assim_dim, infl,
                          filter_options,
                          localization,
                          model, pe,
                          obs, screen=2)

    # Assimilate
    PDAF_caller.assimilate_pdaf(model, obs, pe,
                                assim_dim, localization,
                                filtertype)

    pe.finalize_parallel()


if __name__ == '__main__':
    main()
