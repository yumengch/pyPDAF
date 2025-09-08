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
import pyPDAF

import model
import pdaf_system

class ModelIntegrator:
    """This class implements functions where the model ensemble is integrated"""
    def __init__(self, model_t: model.Model) -> None:
        self.model_t = model_t

    def forward(self, nsteps:int, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration"""
        # When each model task runs one ensemble member,
        # i.e. no need to run each ensemble member sequentially,
        # we call this full parallel implementation
        if da_system.pe.n_modeltasks >= da_system.pe.dim_ens:
            self.forward_full_parallel(nsteps, da_system)
        else:
            self.forward_flexible(da_system)

    def forward_full_parallel(self, nsteps:int, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration when each model task
        (a group of processes with at least one processes) runs one ensemble member,
        i.e. no need to run each ensemble member sequentially for this group of processes,
        we call this full parallel implementation"""
        for j in range(nsteps):
            self.model_t.step(da_system.pe, j)
            da_system.assimilate()

    def forward_flexible(self, da_system:pdaf_system.PDAFsystem) -> None:
        """This function implements the model integration when each model task
        runs multiple ensemble members sequentially, i.e. model tasks is smaller than
        the ensemble size, flexible implementation"""
        timenow = 0.
        doexit = 0
        # full DA system integration loop
        nsteps_da = da_system.steps_for
        while True:
            # model integration
            for j in range(nsteps_da):
                self.model_t.step(da_system.pe, timenow + j)
                da_system.assimilate()
            nsteps_da, timenow, doexit = pyPDAF.get_fcst_info(nsteps_da,
                                                              timenow,
                                                              doexit)
            if doexit == 1:
                break
