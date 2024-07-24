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

import pyPDAF.PDAF as PDAF

class PDAF_system:

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
    def __init__(self, pe:parallelisation.parallelisation, model_ens:list[model.model]) -> None:
        self.pe:parallelisation.parallelisation = pe
        self.model_ens:list[model.model] = model_ens
        
        self.filter_options = filter_options.filter_options()
        self.sv = state_vector.state_vector(model_ens[0], dim_ens=pe.dim_ens)
        self.local = localisation.localisation(sv=self.sv)
        # here, observation only uses the domain observation of the model ensemble.
        # Therefore, only one ensemble member (i.e., model_ens[0]) is passed to obs_factory.
        # In a more complicated system, it is possible to have domain class for model,
        # in which case only the domain object is required here.
        self.obs = obs_factory.obs_factory(self.pe, self.model_ens[0], self.local)
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
            filter_param_i, filter_param_r = self.setEnKFOptions(6, 2)
        else:
            # All other filters
            filter_param_i, filter_param_r = self.setETKFOptions(7, 2)

        status:int = 0
        cltor:collector.collector = collector.collector(self.model_ens[0], self.pe)
        # initialise PDAF filters, communicators, ensemble
        _, _, status = PDAF.init(self.filter_options.filtertype,
                                 self.filter_options.subtype,
                                 0,
                                 filter_param_i,
                                 filter_param_r,
                                 self.pe.comm_model.py2f(),
                                 self.pe.comm_filter.py2f(),
                                 self.pe.comm_couple.py2f(), self.pe.task_id,
                                 self.pe.n_modeltasks, self.pe.filter_pe,
                                 cltor.init_ens_pdaf, screen)

        assert status == 0, f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{self.pe.mype_ens})'

        lfilter = PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1

        # PDAF distribute the initial ensemble to model field
        doexit:int = 0
        prepost:prepost_processing.prepost = prepost_processing.prepost(self.model_ens[0], self.pe)
        for i in range(self.pe.dim_ens_l):
            dist = distributor.distributor(self.model_ens[i])
            self.steps_for, time, doexit, status = PDAF.get_state(self.steps_for, doexit,
                                             dist.next_observation,
                                             dist.distribute_state,
                                             prepost.initial_process,
                                             status)
        # set local domain on each model process
        if self.local.local_filter:
            self.local.set_lim_coords(self.model_ens[i].nx_p, self.model_ens[i].ny_p, self.pe)

    def setEnKFOptions(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray, np.ndarray]:
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

    def setETKFOptions(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray, np.ndarray]:
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

    def assimilate_full_parallel(self) -> None:
        """Assimilation function for the full parallel implementation
        """
        doexit:int = 0
        status:int = 0
        cltor:collector.collector = collector.collector(self.model_ens[0], self.pe)
        prepost:prepost_processing.prepost = prepost_processing.prepost(self.model_ens[0], self.pe)
        dist = distributor.distributor(self.model_ens[0])
        if self.local.local_filter:
            status = \
                   PDAF.omi_assimilate_local(cltor.collect_state,
                                           dist.distribute_state,
                                           self.obs.init_dim_obs_pdafomi,
                                           self.obs.obs_op_pdafomi,
                                           prepost.prepostprocess,
                                           self.local.init_n_domains_pdaf,
                                           self.local.init_dim_l_pdaf,
                                           self.obs.init_dim_obs_l_pdafomi,
                                           self.local.g2l_state_pdaf,
                                           self.local.l2g_state_pdaf,
                                           dist.next_observation)
        else:
            if self.filter_options.filtertype == 8:
                status = \
                       PDAF.omi_assimilate_lenkf(cltor.collect_state,
                                               dist.distribute_state,
                                               self.obs.init_dim_obs_pdafomi,
                                               self.obs.obs_op_pdafomi,
                                               prepost.prepostprocess,
                                               self.obs.localize_covar_pdafomi,
                                               dist.next_observation)
            else:
                status = \
                       PDAF.omi_assimilate_global(cltor.collect_state,
                                               dist.distribute_state,
                                               self.obs.init_dim_obs_pdafomi,
                                               self.obs.obs_op_pdafomi,
                                               prepost.prepostprocess,
                                               dist.next_observation)

        assert status == 0, f'ERROR {status} in PDAF_put_state - stopping! (PE {self.pe.mype_ens})'

    def assimilate_flexible(self) -> None:
        """This function implement the assimilation in flexibble implementation.

        The put_state_XXX functions put model fields into PDAF state vectors using
        i.e. PDAF will collect state vectors from models (from a user-supplied functions p.o.v.).
        When all ensemble members are collected, the PDAF distribute each model fields 
        """
        doexit:int = 0
        status:int = 0
        cltor: collector.collector
        prepost:prepost_processing.prepost = prepost_processing.prepost(self.model_ens[0], self.pe)
        if self.local.local_filter:
            for i in range(self.pe.dim_ens_l):
                cltor = collector.collector(self.model_ens[i], self.pe)
                status = PDAF.omi_put_state_local(cltor.collect_state,
                    self.obs.init_dim_obs_pdafomi,
                    self.obs.obs_op_pdafomi,
                    prepost.prepostprocess,
                    self.local.init_n_domains_pdaf,
                    self.local.init_dim_l_pdaf,
                    self.obs.init_dim_obs_l_pdafomi,
                    self.local.g2l_state_pdaf,
                    self.local.l2g_state_pdaf)
        else:
            if self.filter_options.filtertype == 8:
                for i in range(self.pe.dim_ens_l):
                    cltor = collector.collector(self.model_ens[i], self.pe)
                    status = PDAF.omi_put_state_lenkf(cltor.collect_state,
                                  self.obs.init_dim_obs_pdafomi, self.obs.obs_op_pdafomi,
                                  prepost.prepostprocess,
                                  self.obs.localize_covar_pdafomi)
            else:
                for i in range(self.pe.dim_ens_l):
                    cltor = collector.collector(self.model_ens[i], self.pe)
                    status = PDAF.omi_put_state_global(cltor.collect_state,
                                  self.obs.init_dim_obs_pdafomi, self.obs.obs_op_pdafomi,
                                  prepost.prepostprocess)

        for i in range(self.pe.dim_ens_l):
            dist = distributor.distributor(self.model_ens[i])
            time:float
            self.steps_for, time, doexit, status = PDAF.get_state(self.steps_for, doexit,
                                              dist.next_observation,
                                              dist.distribute_state,
                                              prepost.prepostprocess,
                                              status)

        assert status == 0, f'ERROR {status} in PDAF_put_state - stopping! (PE {self.pe.mype_ens})'


    def finalise(self) -> None:
        PDAF.print_info(11)
        if (self.pe.mype_ens == 0): PDAF.print_info(3)
        PDAF.deallocate()

