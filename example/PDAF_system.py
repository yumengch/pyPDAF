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
import step.shift

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
        self.initial_ensemble_filename:str = config.init_ens_path
        
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
                                 self.init_ens_pdaf, screen)

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

    def setEnKFOptions(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray[int], np.ndarray[float]]:
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

    def setETKFOptions(self, dim_pint:int, dim_preal:int) -> tuple[np.ndarray[int], np.ndarray[float]]:
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

    def init_ens_pdaf(self, filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray, ens_p:np.ndarray,
                      status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """Here, only ens_p variable matters while dim_p and dim_ens defines the
        size of the variables. uinv, state_p are not used in this example.

        status_pdaf is used to handle errors which we will not do it in this example.
        """
        # The initial ensemble is read here and will be distributed to 
        # the model in the PDAF.get_state functtion by a distributor.

        # If your ensemble is read from a restart file, you can simply set this
        # function as a dummy function without doing anything
        # However, you still need to set a distributor to call PDAF.get_state, which
        # does nothing as well. 
        nx_p:int = self.model_ens[0].nx_p
        offset:int = self.pe.mype_filter*nx_p
        for i in range(dim_ens):
            ens_p[:, i] = np.loadtxt(self.initial_ensemble_filename.format(i=i+1))[:, offset:offset+nx_p].ravel()
        return state_p, uinv, ens_p, status_pdaf

    def forward(self, nsteps:int) -> None:
        # When each model task runs one ensemble member,
        # i.e. no need to run each ensemble member sequentially,
        # we call this full parallel implementation
        if self.pe.dim_ens_l == 1:
            self.forward_full_parallel(nsteps)
        else:
            self.forward_flexible(nsteps)

    def forward_full_parallel(self, nsteps:int) -> None:
        for i in range(nsteps):
            self.model_ens[0].field_p = step.shift.step(self.model_ens[0].field_p, self.pe, i+1, config.USE_PDAF)
            if config.USE_PDAF:
                self.assimilate_full_parallel()

    def forward_flexible(self, nsteps:int) -> None:
        # create directory to
        current_step = 0
        # full DA system integration loop
        while current_step < nsteps:
            if config.USE_PDAF:
                # model integration
                for _ in range(self.steps_for):
                    for i in range(self.pe.dim_ens_l):
                        self.model_ens[i].field_p = step.shift.step(self.model_ens[i].field_p, self.pe, current_step + 1, config.USE_PDAF)
                    current_step += 1

                self.assimilate_flexible()
            else:
                for i in range(self.pe.dim_ens_l):
                    self.model_ens[i].field_p = step.shift.step(self.model_ens[i].field_p, self.pe, current_step + 1, config.USE_PDAF)
            current_step += 1


    def assimilate_full_parallel(self) -> None:
        """Assimilation function for the full parallel implementation
        """
        doexit:int = 0
        status:int = 0
        cltor:collector.collector = collector.collector(self.model_ens[0])
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
                cltor = collector.collector(self.model_ens[i])
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
                    cltor = collector.collector(self.model_ens[i])
                    status = PDAF.omi_put_state_lenkf(cltor.collect_state,
                                  self.obs.init_dim_obs_pdafomi, self.obs.obs_op_pdafomi,
                                  prepost.prepostprocess,
                                  self.obs.localize_covar_pdafomi)
            else:
                for i in range(self.pe.dim_ens_l):
                    cltor = collector.collector(self.model_ens[i])
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
