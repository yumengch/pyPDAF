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
import DAS


from Localization import Localization
from AssimilationDimensions import AssimilationDimensions
from FilterOptions import FilterOptions
from Inflation import Inflation


def main():
    USE_PDAF = True
    nt = 100

    if USE_PDAF:
        pe = parallelization(dim_ens=4, n_modeltasks=4, screen=2)

    # Initial Screen output
    if (pe.mype_world == 0):
        print('+++++ PDAF online mode +++++')
        print('2D model without parallelization')

    # Initialize model
    model = Model.Model((18, 36), nt=18, pe=pe)
    obs = []
    for typename in ['A', 'B']:
        obs.append(OBS.OBS(typename, pe.mype_filter,
                           model.nx, 1, 2, 0.5))

    das = DAS.DAS(pe, model, obs, screen=2)
    assim_dim = AssimilationDimensions(das.model, dim_ens=das.pe.n_modeltasks)
    filter_options = FilterOptions(filtertype=6, subtype=0)
    filter_options.setTransformTypes(type_trans=0, type_sqrt=0,
                                     incremental=0, covartype=1,
                                     rank_analysis_enkf=0)
    infl = Inflation(type_forget=0, forget=1.)
    localization = Localization(loc_weight=0, local_range=0,
                                srange=0)
    das.init(assim_dim, filter_options, infl, localization)

    for it in range(nt):
        das.forward(it, USE_PDAF)

    pe.finalize_parallel()


if __name__ == '__main__':
    main()
