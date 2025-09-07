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
import numpy as np

import config
import parallelisation

class Model:
    """Model information in PDAF

    Attributes
    ----------
    field_p : ndarray
        PE-local model field
    nx : int
        number of grid points in x-direction
    nx_p : int
        number of grid points in x-direction on local PE
    ny : int
        number of grid points in y-direction
    ny_p : int
        number of grid points in y-direction on local PE
    current_step : int
        current time step
    total_steps : int
        total number of time steps
    """

    def __init__(self, pe:parallelisation.Parallelisation) -> None:
        """constructor

        Parameters
        ----------
        nx : ndarray
            integer array for grid size in x-dimension
        ny : ndarray
            integer array for grid size in y-dimension
        nt : int
            total number of time steps
        pe : `parallelization.parallelization`
            parallelization object
        """
        # model size
        self.nx:int = config.nx
        self.ny:int = config.ny
        # model domain for each CPU (model domain decomposition)
        self.nx_p:int
        self.ny_p:int
        self.nx_p, self.ny_p = self.get_local_domain(pe)
        # model time steps
        self.total_steps:int = config.nsteps
        # model field
        self.field_p:np.ndarray = np.zeros((self.ny_p, self.nx_p))

    def get_local_domain(self, pe:parallelisation.Parallelisation) -> tuple[int, int]:
        """Compute local-PE domain size/domain decomposition

        Parameters
        ----------
        pe : `parallelization.Parallelisation`
            parallelization object
        """
        nx_p:int = self.nx
        ny_p:int = self.ny

        assert self.nx % pe.npes_model == 0, f'...ERROR: Invalid number of' \
               f'processes: {pe.npes_model}...'
        # we parallelise the domain column by column
        nx_p = self.nx//pe.npes_model

        return nx_p, ny_p

    def init_field(self, mype_model:int) -> None:
        """initialise PE-local model field

        Parameters
        ----------
        mype_model : int
            rank of the process in model communicator
        """
        offset = self.nx_p*mype_model
        self.field_p[:] = np.loadtxt(config.init_truth_path)[:, offset:self.nx_p + offset]

    def step(self, pe:parallelisation.Parallelisation, timenow:float) -> None:
        """shifting model forward 'integration'

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object from example
        """
        if pe.task_id == 1 and pe.mype_model == 0:
            output_str = f'model step: {timenow}'
            log.logger.info(output_str)

        self.field_p = np.roll(self.field_p, 1, axis=0)


    def print_info(self, pe:parallelisation.Parallelisation) -> None:
        """print model info

        Parameters
        ----------
        pe : `parallelization.Parallelisation`
            parallelization object
        """
        if pe.task_id == 1 and pe.mype_model == 0:
            log.logger.info('MODEL-side: INITIALIZE PARALLELIZED Shifting model MODEL')
            output_str = f'Grid size: {self.nx} x {self.ny}'
            log.logger.info(output_str)
            output_str = f'Time steps {self.total_steps}'
            log.logger.info(output_str)
            output_str = f'-- Domain decomposition over {pe.npes_model} PEs'
            log.logger.info(output_str)
            output_str = f'-- local domain sizes: {self.nx_p} x {self.ny_p}'
            log.logger.info(output_str)
