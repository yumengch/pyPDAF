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
    def __init__(self, pe:parallelisation.Parallelisation, model_ens:list[model.Model]) -> None:
        self.pe = pe
        self.model_ens = model_ens

        self.filter_options = filter_options.FilterOptions()
        self.sv = state_vector.StateVector(model_ens[0], dim_ens=pe.dim_ens)
        self.local = localisation.Localisation(sv=self.sv)
        # here, observation only uses the domain observation of the model ensemble.
        # Therefore, only one ensemble member (i.e., model_ens[0]) is passed to obs_factory.
        # In a more complicated system, it is possible to have domain class for model,
        # in which case only the domain object is required here.
        self.obs = obs_factory.ObsFactory(self.pe, self.model_ens[0], self.local)
        # initial time step
        self.steps_for = config.init_step

    def init_pdaf(self, screen:int) -> None:
        """constructor

        Parameters
        ----------
        screen : int
            verbosity of PDAF screen output
        """
        filter_param_i:np.ndarray
        filter_param_r:np.ndarray
        if self.filter_options.filtertype == 2:
            # EnKF with Monte Carlo init
            filter_param_i, filter_param_r = self.set_enkf_options(6, 2)
        else:
            # All other filters
            filter_param_i, filter_param_r = self.set_etkf_options(7, 2)

        status:int = 0
        cltor = collector.Collector(self.model_ens[0], self.pe)
        # initialise PDAF filters, communicators, ensemble
        _, _, status = pyPDAF.init(self.filter_options.filtertype, self.filter_options.subtype,
                                   0, filter_param_i, len(filter_param_i), filter_param_r, 2,
                                   cltor.init_ens_pdaf, screen)

        assert status == 0, f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{self.pe.mype_ens})'

        pyPDAF.PDAFomi.init(self.obs.nobs)
        pyPDAF.PDAFomi.init_local()
        lfilter = pyPDAF.PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1

        # PDAF distribute the initial ensemble to model field
        prepost = prepost_processing.Prepost(self.model_ens[0], self.pe)
        for i in range(self.pe.dim_ens_l):
            dist = distributor.Distributor(self.model_ens[i])
            status = pyPDAF.init_forecast(dist.next_observation,
                                          dist.distribute_state_ini,
                                          prepost.initial_process,
                                          status)
        # set local domain on each model process
        if self.local.local_filter:
            self.local.set_lim_coords(self.model_ens[0].nx_p, self.model_ens[0].ny_p, self.pe)

    def set_enkf_options(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray, np.ndarray]:
        """set ensemble kalman filter options

        Parameters
        ----------
        dim_pint : int
            size of integer filter options
        dim_preal : int
            size of float filter options
        """
        filter_param_i:np.ndarray = np.zeros(dim_pint, dtype=np.intc)
        filter_param_r:np.ndarray = np.zeros(dim_preal)

        filter_param_i[0] = self.sv.dim_state_p
        filter_param_i[1] = self.sv.dim_ens
        filter_param_i[2] = self.filter_options.rank_analysis_enkf
        filter_param_i[3] = self.filter_options.incremental
        filter_param_i[4] = 0

        filter_param_r[0] = self.filter_options.forget

        return filter_param_i, filter_param_r

    def set_etkf_options(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray, np.ndarray]:
        """Summary

        Parameters
        ----------
        dim_pint : int
            size of integer filter options
        dim_preal : int
            size of float filter options
        """
        filter_param_i:np.ndarray = np.zeros(dim_pint, dtype=np.intc)
        filter_param_r:np.ndarray = np.zeros(dim_preal)

        filter_param_i[0] = self.sv.dim_state_p
        filter_param_i[1] = self.sv.dim_ens
        filter_param_i[2] = 0
        filter_param_i[3] = self.filter_options.incremental
        filter_param_i[4] = self.filter_options.type_forget
        filter_param_i[5] = self.filter_options.type_trans
        filter_param_i[6] = self.filter_options.type_sqrt

        filter_param_r[0] = self.filter_options.forget

        return filter_param_i, filter_param_r

    def assimilate(self, i:int) -> None:
        """Calling assimilation function of PDAF

        Parameters
        ----------
        i : int
            index of the ensemble in current model task.
        """
        status:int = 0
        cltor = collector.Collector(self.model_ens[i], self.pe)
        prepost = prepost_processing.Prepost(self.model_ens[i], self.pe)
        dist = distributor.Distributor(self.model_ens[i])
        status = \
                pyPDAF.assimilate(cltor.collect_state,
                                  dist.distribute_state,
                                  self.obs.init_dim_obs_pdafomi,
                                  self.obs.obs_op_pdafomi,
                                  self.local.init_n_domains_pdaf,
                                  self.local.init_dim_l_pdaf,
                                  self.obs.init_dim_obs_l_pdafomi,
                                  prepost.prepostprocess,
                                  dist.next_observation, status)

        assert status == 0, f'ERROR {status} in PDAF_put_state - stopping! (PE {self.pe.mype_ens})'

    def finalise(self) -> None:
        """finalise PDAF
        """
        pyPDAF.PDAF.print_info(11)
        if self.pe.mype_ens == 0:
            pyPDAF.PDAF.print_info(3)
        pyPDAF.deallocate()
