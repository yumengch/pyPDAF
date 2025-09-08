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
import pyPDAF

import collector
import config
import filter_options
import localisation
import model
import obs_factory
import parallelisation
import prepost_processing
import state_vector


class PDAFsystem:

    """PDAF system

    Attributes
    ----------
    pe : parallelisation.parallelisation
        parallelisation instance
    model_grid : model.model_grid
        list of model instances
    sv : state_vector.state_vector
        state vector
    local : localisation.localisation
        localisation
    obs : obs_factory.obs_factory
        observation factory
    filter_options : filter_options.filter_options
        filter options
    """
    def __init__(self, pe:parallelisation.Parallelisation,
                 model_grid:model.ModelGrid) -> None:
        self.pe:parallelisation.Parallelisation = pe
        self.model_grid:model.ModelGrid = model_grid

        self.filter_options = filter_options.FilterOptions()
        self.sv = state_vector.StateVector(model_grid,
                                            dim_ens=pe.dim_ens)
        self.local = localisation.Localisation(sv=self.sv)
        # here, observation only uses the domain observation of the model ensemble.
        # Therefore, only one ensemble member (i.e., model_ens[0]) is passed to obs_factory.
        # In a more complicated system, it is possible to have domain class for model,
        # in which case only the domain object is required here.
        self.obs = obs_factory.ObsFactory(self.pe,
                                           self.model_grid, self.local)
        # initial time step
        self.steps_for = config.init_step

    def init_pdaf(self, screen:int) -> None:
        """constructor

        Parameters
        ----------
        screen : int
            verbosity of PDAF screen output
        """
        filter_param_i = np.array([self.sv.dim_state_p, self.sv.dim_ens], dtype=np.intc)
        filter_param_r = np.array([self.filter_options.forget, ])
        status:int = 0
        cltor:collector.Collector = collector.Collector(
            self.model_grid, self.pe)
        # initialise PDAF filters, communicators, ensemble
        # initialise PDAF filters, communicators, ensemble
        _, _, status = pyPDAF.init(self.filter_options.filtertype, self.filter_options.subtype,
                                   0, filter_param_i, len(filter_param_i), filter_param_r, 2,
                                   cltor.init_ens_pdaf, screen)

        assert status == 0, f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{self.pe.mype_ens})'

        pyPDAF.PDAFomi.init(self.obs.nobs)

        lfilter:int = pyPDAF.PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1
        # set local domain on each model process
        if self.local.local_filter:
            pyPDAF.PDAFomi.init_local()
            self.local.set_lim_coords(self.model_grid.nx_p,
                                      self.model_grid.ny_p, self.pe)


    def assimilate(self) -> None:
        """offline assimilation function.

        The put_state_XXX functions execute the DA scheme.
        Here, the state vector is read from data in `init_ens_pdaf`.
        """
        status:int = 0
        prepost = prepost_processing.Prepost(self.model_grid, self.pe)
        pyPDAF.assim_offline(self.obs.init_dim_obs_pdafomi,
                             self.obs.obs_op_pdafomi,
                             self.local.init_n_domains_pdaf,
                             self.local.init_dim_l_pdaf,
                             self.obs.init_dim_obs_l_pdafomi,
                             prepost.prepostprocess, status)


        assert status == 0, f'ERROR {status} in PDAF_put_state' \
            f' - stopping! (PE {self.pe.mype_ens})'

    def finalise(self) -> None:
        """finalise PDAF system"""
        pyPDAF.PDAF.print_info(11)
        if self.pe.mype_ens == 0:
            pyPDAF.PDAF.print_info(3)
        pyPDAF.PDAF.deallocate()

