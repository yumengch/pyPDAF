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

import log
import model
import sfields
import transforms
import io_nemo

class Prepost:
    """PDAF pre/post-processing callbacks around the analysis step.

    On the first call, the forecast ensemble is saved as background and optional
    state transforms/limits are applied before analysis. On the second call,
    transforms are reversed and NEMO ASM background/increment files are written.

    Attributes
    ----------
    model : `model.model`
        model object
    pe : `parallelisation.parallelisation`
        parallelisation object
    -------
    """
    def __init__(self, writer: io_nemo.Writer, nemo_grid: model.NemoDomain,
                 fields: dict[str, sfields.StateField]) -> None:
        self.first = True
        self.model = nemo_grid
        self.fields = fields
        self.writer = writer

    def initial_process(self, ens_p:np.ndarray) -> np.ndarray:
        """Save background ensemble and transform fields before assimilation."""
        self.first = False
        log.info('Analyze forecast state ensemble')
        self.model.bkg_ens = ens_p.copy()
        ens_p = transforms.transform_field_mv(1, ens_p, self.fields, 11)
        return ens_p

    def postprocess(self, dim_ens:int, ens_p:np.ndarray) -> np.ndarray:
        """Reverse transforms and write analysis/background increment files."""
        log.info('Analyze assimilated state ensemble')
        log.info('--- Write ensemble of increments')
        ens_p = transforms.transform_field_mv(2, ens_p, self.fields, 0)
        for i in range(1, dim_ens + 1):
            self.writer.write_asmdin_mv(i, ens_p[:, i-1], self.fields, self.model)
            self.writer.write_asminc_mv(i, ens_p[:, i-1] - self.model.bkg_ens[:, i-1],
                                        self.fields, self.model)
        return ens_p


    def prepostprocess(self, _step:int, _dim_p:int, dim_ens:int, _dim_ens_p:int,
                       _dim_obs_p:int, state_p:np.ndarray, uinv:np.ndarray,
                       ens_p:np.ndarray, _flag:int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Dispatch PDAF's pre/post callback to the first or second phase."""
        if self.first:
            ens_p = self.initial_process(ens_p)
        else:
            ens_p = self.postprocess(dim_ens, ens_p)
        return state_p, uinv, ens_p
