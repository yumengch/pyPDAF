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
import pyPDAF.PDAF.PDAFomi as PDAFomi
import PDAF_caller

from Localization import Localization
from AssimilationDimensions import AssimilationDimensions
from FilterOptions import FilterOptions
from Inflation import Inflation


class DAS:

    """Data assimilation system

    Attributes
    ----------
    assim_dim : `AssimilationDimensions.AssimilationDimensions`
        an object of AssimilationDimensions
    filter_options : `FilterOptions.FilterOptions`
        filtering options
    infl : `Inflation.Inflation`
        inflation object
    localization : `Localization.Localization`
        localization object
    model : `Model.Model`
        model object
    obs : `OBS.OBS`
        observation object
    pe : `parallelization.parallelization`
        parallelization object
    screen : int
        verbosity of PDAF screen output
    """

    def __init__(self, pe, model, obs, screen):
        """constructor

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        model : `Model.Model`
            model object
        obs : `OBS.OBS`
            observation object
        screen : int
            verbosity of PDAF screen output
        """
        self.pe = pe
        self.model = model
        self.obs = obs
        self.screen = screen

    def init(self):
        """initialise DA system
        """
        # init model
        self.model.init_field('inputs_online/true_initial.txt',
                              self.pe.mype_model)

        # init observations
        PDAFomi.init(len(self.obs))

        # Initialize PDAF
        self.assim_dim = AssimilationDimensions(self.model,
                                                dim_ens=self.pe.n_modeltasks)
        self.filter_options = FilterOptions(filtertype=6, subtype=0)
        self.filter_options.setTransformTypes(type_trans=0,
                                              type_sqrt=0,
                                              incremental=0,
                                              covartype=1,
                                              rank_analysis_enkf=0)
        self.infl = Inflation(type_forget=0, forget=1.)
        self.localization = Localization(loc_weight=0, local_range=0,
                                         srange=0)
        PDAF_caller.init_pdaf(self.assim_dim, self.infl,
                              self.filter_options,
                              self.localization,
                              self.model, self.pe,
                              self.obs, self.screen)

    def forward(self, step, USE_PDAF):
        """time forward DA system

        Parameters
        ----------
        step : int
            current time step
        USE_PDAF : bool
            whether PDAF is used
        """
        self.model.step(self.pe, step, USE_PDAF)
        if USE_PDAF:
            PDAF_caller.assimilate_pdaf(self.model, self.obs, self.pe,
                                        self.assim_dim, self.localization,
                                        self.filter_options.filtertype)
