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
import io_nemo
import sfields
import transforms

class Collector:
    """Load the background/forecast ensemble from model restart files.

    PDAF calls :meth:`init_ens_pdaf` during initialisation. In this offline
    workflow the model is not advanced by Python, so the callback reads each
    ensemble member from NEMO restart files and packs configured fields into
    PDAF's compact wet-point state vector.
    """
    def __init__(self, reader: io_nemo.Reader, wet_grid: model.WetGrid,
                 fields: dict[str, sfields.StateField]) -> None:
        # initialise the model instance
        self.wet_grid = wet_grid
        self.reader = reader
        self.fields = fields

    def init_ens_pdaf(self, _filtertype:int, _dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray, ens_p:np.ndarray,
                      status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """read forecast ensemble for DA

        In offline DA, the ensemble is read from disk.
        One does not need collect/distribute ensemble to model.

        Parameters
        ----------
        filtertype: int
            filter type
        dim_p: int
            dimension of state vector on local processor
        dim_ens: int
            ensemble size
        state_p: np.ndarray
            state vector on local processor
        uinv: np.ndarray
            inverse of covariance matrix
        ens_p: np.ndarray
            ensemble matrix on local processor
        """
        # The ensemble for offline DA is read here.

        log.info("Initialize ensemble from list of restart files")

        for i in range(1, dim_ens + 1):
            # read model variables
            data_iterator = self.reader.read_restart(i, self.fields)
            # convert model variables to state vector
            for data, off, dim, variable in data_iterator:
                ens_p[off : off + dim, i - 1] = transforms.field2state(data, self.wet_grid)
                if dim > 0:
                    log.info(
                        f"PE:{self.reader.pe.mype_filter} Min and max for {variable} :  "
                        f"{ens_p[off : off + dim, i - 1].min()} "
                        f"{ens_p[off : off + dim, i - 1].max()}", ranks=self.reader.pe.mype_filter
                    )
        state_p[:] = ens_p.mean(axis=1)
        return state_p, uinv, ens_p, status_pdaf
