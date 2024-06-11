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
import copy

import model.shift


class Model:

    """Model information in PDAF

    Attributes
    ----------
    field_p : ndarray
        PE-local model field
    nx : ndarray
        integer array for grid size
    nx_p : ndarray
        integer array for PE-local grid size
    total_steps : int
        total number of time steps
    """

    def __init__(self, nx, nt, pe):
        """constructor

        Parameters
        ----------
        nx : ndarray
            integer array for grid size
        nt : int
            total number of time steps
        pe : `parallelization.parallelization`
            parallelization object
        """
        # model size
        self.nx = list(nx)
        # model size for each CPU
        self.get_nxp(pe)
        # model time steps
        self.total_steps = nt

    def get_nxp(self, pe):
        """Compute local-PE domain size/domain decomposition

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        """
        self.nx_p = copy.copy(self.nx)

        try:
            assert self.nx[-1] % pe.npes_model == 0
            self.nx_p[-1] = self.nx[-1]//pe.npes_model
        except AssertionError:
            print((f'...ERROR: Invalid number of'
                   f'processes: {pe.npes_model}...'))
            pe.abort_parallel()

    def init_field(self, filename, mype_model):
        """initialise PE-local model field

        Parameters
        ----------
        filename : string
            input filename
        mype_model : int
            rank of the process in model communicator
        """
        # model field
        self.field_p = np.zeros(self.nx_p)
        offset = self.nx_p[-1]*mype_model
        self.field_p = np.loadtxt(
                                    filename
                                    )[:, offset:self.nx_p[-1] + offset]

    def step(self, pe, step, USE_PDAF):
        """step model forward

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        step : int
            current time step
        USE_PDAF : bool
            whether PDAF is used
        """
        model.shift.step(self, pe, step, USE_PDAF)

    def printInfo(self, USE_PDAF, pe):
        """print model info

        Parameters
        ----------
        USE_PDAF : bool
            whether PDAF is used
        pe : `parallelization.parallelization`
            parallelization object
        """
        do_print = USE_PDAF and pe.mype_model == 0
        do_print = do_print or \
            (pe.task_id == 1 and pe.mype_model == 0 and not USE_PDAF)
        if do_print:
            print('MODEL-side: INITIALIZE PARALLELIZED Shifting model MODEL')
            print(f'Grid size: {self.nx}')
            print(f'Time steps {self.total_steps}')
            print(f'-- Domain decomposition over {pe.npes_model} PEs')
            print(f'-- local domain sizes (nx_p): {self.nx_p}')
