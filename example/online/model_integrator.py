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
import log
from mpi4py import MPI
import numpy as np
import pyPDAF

import config
import model
import pdaf_system

class ModelIntegrator:
    """This class implements functions where the model ensemble is integrated"""
    def __init__(self, model_ens: list[model.Model]) -> None:
        self.model_ens: list[model.Model] = model_ens

    def forward(self, nsteps:int, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration"""
        # When each model task runs one ensemble member,
        # i.e. no need to run each ensemble member sequentially,
        # we call this full parallel implementation
        if da_system.pe.dim_ens_l == 1:
            self.forward_full_parallel(nsteps, da_system)
        else:
            self.forward_flexible(nsteps, da_system)

    def forward_full_parallel(self, nsteps:int, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration when each model task
        (a group of processes with at least one processes) runs one ensemble member,
        i.e. no need to run each ensemble member sequentially for this group of processes,
        we call this full parallel implementation"""
        for i in range(nsteps):
            self.model_ens[0].field_p = self.step(self.model_ens[0].field_p,
                                                  da_system.pe, i+1)
            da_system.assimilate(0)

    def forward_flexible(self, nsteps:int, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration when each model task
        runs multiple ensemble members sequentially, i.e. model tasks is smaller than
        the ensemble size, flexible implementation"""
        # create directory to
        current_step = 0
        timenow = 0.
        doexit = 0
        # full DA system integration loop
        while current_step < nsteps:
            nsteps, timenow, doexit = pyPDAF.get_fcst_info(da_system.steps_for,
                                                            timenow,
                                                            doexit)
            if doexit == 1:
                break

            for i in range(da_system.pe.dim_ens_l):
                pyPDAF.PDAF.set_debug_flag(10)
                # model integration
                for j in range(nsteps):
                    self.model_ens[i].field_p = self.step(
                        self.model_ens[i].field_p,
                        da_system.pe, current_step + j + 1)
                    da_system.assimilate(i)
            current_step += nsteps


    def step(self, field_p:np.ndarray, pe, step:int) -> np.ndarray:
        """shifting model forward 'integration'

        Parameters
        ----------
        field_p : `np.ndarray`
            model field
        pe : `parallelization.parallelization`
            parallelization object from example
        step : int
            current time step
        use_pdaf : bool
            whether PDAF is performed
        """
        if pe.task_id == 1 and pe.mype_model == 0:
            output_str = f'model step: {step}'
            log.logger.info(output_str)

        field_p = np.roll(field_p, 1, axis=0)

        return field_p
