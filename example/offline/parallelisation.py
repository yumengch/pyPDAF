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


class parallelisation:

    """Parallelisation for offline DA system

    Attributes
    ----------
    dim_ens : int
        total ensemble size

    comm_ens : `MPI.Comm`
        ensemble communicator
    npes_ens : int
        number of global PEs
    mype_ens : int
        rank of ens communicator
    local_npes_model : int
        number of model PEs used by each model task
    task_id : int
        model task id of the current processor
    comm_model : `MPI.Comm`
        model communicator
    npes_model : int
        number of model PEs
    mype_model : int
        rank of model communicator
    dim_ens_l  : int
        number of ensemble members per model task
        This is the same as dim_ens
        for each model task in parallel
    all_dim_ens_l : np.ndarray
        number of ensemble members across all PEs.
    filterpe : bool
        whether the PE is used for filter
    comm_filter : `MPI.Comm`
        filter communicator
    npes_filter : int
        number of filter PEs
    mype_filter : int
        rank of filter communicator
    comm_couple : `MPI.Comm`
        model and filter coupling communicator
    """

    def __init__(self, dim_ens:int) -> None:
        """Init the parallization required by PDAF

        Parameters
        ----------
        dim_ens : int
            number of ensemble members
        """
        # total number of ensemble members. This can be larger than
        # number of parallel model tasks where the 'flexible' implementation
        # is used
        self.dim_ens:int = dim_ens

        # Initialize model communicator, its size and the process rank
        # Here the same as for MPI_COMM_WORLD
        self.comm_ens:MPI.Comm
        self.npes_ens:int
        self.mype_ens:int
        self.comm_ens, self.npes_ens, self.mype_ens = \
            self.init_parallel()

        # Initialize communicators for ensemble evaluations
        if self.mype_ens == 0:
            log.logger.info('Initialize communicators '
                            'for assimilation with PDAF')

        # get model communicator
        self.comm_model:MPI.Comm = self.comm_ens
        self.npes_model:int = self.npes_ens
        self.mype_model:int = self.mype_ens

        log.logger.info(f'MODEL: mype(w)= {self.mype_ens}; '
                        f'mype(m)= {self.mype_model}; '
                        f'npes(m)= {self.npes_model}')

        # Generate communicator for filter
        self.filter_pe:bool = True
        self.comm_filter:MPI.Comm = self.comm_ens
        self.npes_filter:int = self.npes_ens
        self.mype_filter:int = self.mype_ens

        # Generate communicators for communication
        self.comm_couple:MPI.Comm = self.comm_ens

        self.print_info()

    def init_parallel(self) -> tuple[MPI.Comm, int, int]:
        """Initialize MPI

        Routine to initialize MPI.

        Here, we also return the communicator for the full ensemble
        as well as its size and the rank of the current processor in
        the communicator.

        In this simple example, the full ensemble communicator is the
        MPI_COMM_WORLD.
        """
        if not MPI.Is_initialized():
            MPI.Init()

        return MPI.COMM_WORLD, MPI.COMM_WORLD.Get_size(), \
            MPI.COMM_WORLD.Get_rank()

    def print_info(self) -> None:
        """print parallelization info
        """
        # *** local variables ***
        #  Rank and size in COMM_couple
        mype_couple = self.comm_couple.Get_rank()
        #  Variables for communicator-splitting
        color_couple = self.mype_model + 1

        if (self.mype_ens == 0):
            log.logger.info('PE configuration:')
            log.logger.info('ens     filter       model        couple   filterPE')
            log.logger.info('rank    rank  task   rank task       rank       T/F')
            log.logger.info('-----------------------------------------------------')
        MPI.COMM_WORLD.Barrier()
        log.logger.info(f'{self.mype_ens},      {self.mype_filter},'
                f'    1      {self.mype_model},'
                f'   {color_couple},         {mype_couple},       {self.filter_pe}')
        MPI.COMM_WORLD.Barrier()
        if (self.mype_ens == 0):
            log.logger.info('')

    def finalize_parallel(self) -> None:
        """Finalize MPI
        """
        MPI.COMM_WORLD.Barrier()
        MPI.Finalize()
