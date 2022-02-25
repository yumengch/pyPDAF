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

    def init(self, assim_dim, filter_options, infl, localization):
        """initialise DA system

        Parameters
        ----------
        assim_dim : `AssimilationDimensions.AssimilationDimensions`
            an object of AssimilationDimensions
        filter_options : `FilterOptions.FilterOptions`
            filtering options
        infl : `Inflation.Inflation`
            inflation object
        localization : `Localization.Localization`
            a localization info obejct
        """
        # init model
        self.model.init_field('inputs_online/true_initial.txt',
                              self.pe.mype_model)

        # init observations
        PDAFomi.init(len(self.obs))

        # Initialize PDAF
        self.assim_dim = assim_dim
        self.filter_options = filter_options
        self.infl = infl
        self.localization = localization
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
            _ = PDAF_caller.assimilate_pdaf(self.model, self.obs, self.pe,
                                            self.assim_dim, self.localization,
                                            self.filter_options.filtertype)
