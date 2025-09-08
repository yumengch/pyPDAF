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

import collector
import config
import distributor
import filter_options
import localisation
import model
import obs_factory
import parallelisation
import prepost_processing
import state_vector

import pyPDAF

class PDAFsystem:

    """PDAF system

    Attributes
    ----------
    pe : parallelisation.parallelisation
        parallelisation instance
    model_ens : list[model.model]
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
    def __init__(self, pe:parallelisation.Parallelisation, model_ens:model.Model) -> None:
        self.pe = pe
        self.model_ens = model_ens

        self.filter_options = filter_options.FilterOptions()
        self.sv = state_vector.StateVector(model_ens, dim_ens=pe.dim_ens)
        self.local = localisation.Localisation(sv=self.sv)
        # here, observation only uses the domain observation of the model ensemble.
        # In a more complicated system, it is possible to have domain class for model,
        # in which case only the domain object is required here.
        self.obs = obs_factory.ObsFactory(self.pe, self.model_ens, self.local)
        # initial time step
        self.steps_for = config.init_step
        self.cltor = collector.Collector(self.model_ens, self.pe)
        self.prepost = prepost_processing.Prepost(self.model_ens, self.pe)
        self.dist = distributor.Distributor(self.model_ens)

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
        # initialise PDAF filters, communicators, ensemble
        _, _, status = pyPDAF.init(self.filter_options.filtertype, self.filter_options.subtype,
                                   0, filter_param_i, 2, filter_param_r, 1,
                                   self.cltor.init_ens_pdaf, screen)

        assert status == 0, f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{self.pe.mype_ens})'

        pyPDAF.PDAFomi.init(self.obs.nobs)
        pyPDAF.PDAFomi.init_local()
        lfilter = pyPDAF.PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1

        # PDAF distribute the initial ensemble to model field
        status = pyPDAF.init_forecast(self.dist.next_observation,
                                      self.dist.distribute_state,
                                      self.prepost.initial_process,
                                      status)
        # set local domain on each model process
        if self.local.local_filter:
            self.local.set_lim_coords(self.model_ens.nx_p, self.model_ens.ny_p, self.pe)

    def assimilate(self) -> None:
        """Calling assimilation function of PDAF

        Parameters
        ----------
        i : int
            index of the ensemble in current model task.
        """
        status:int = 0

        status = \
                pyPDAF.assimilate(self.cltor.collect_state,
                                  self.dist.distribute_state,
                                  self.obs.init_dim_obs_pdafomi,
                                  self.obs.obs_op_pdafomi,
                                  self.local.init_n_domains_pdaf,
                                  self.local.init_dim_l_pdaf,
                                  self.obs.init_dim_obs_l_pdafomi,
                                  self.prepost.prepostprocess,
                                  self.dist.next_observation, status)

        assert status == 0, f'ERROR {status} in PDAF_put_state - stopping! (PE {self.pe.mype_ens})'

    def finalise(self) -> None:
        """finalise PDAF
        """
        pyPDAF.PDAF.print_info(11)
        if self.pe.mype_ens == 0:
            pyPDAF.PDAF.print_info(3)
        pyPDAF.deallocate()
