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

    """Summary

    Attributes
    ----------
    n_modeltasks : int
        number of model tasks run in parallel
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
        This is the same as dim_ens//n_modeltasks
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

    def __init__(self, dim_ens:int, n_modeltasks:int) -> None:
        """Init the parallization required by PDAF

        Parameters
        ----------
        dim_ens : TYPE
            Description
        n_modeltasks : int
            Number of model tasks/ ensemble size
            This parameter should be the same
            as the number of PEs
        """

        # number of model tasks run in parallel 
        # this is limited by the number of processors
        self.n_modeltasks:int = n_modeltasks
        # total number of ensemble members. This can be larger than 
        # number of parallel model tasks where the 'flexible' implementation 
        # is used
        self.dim_ens:int = dim_ens

        # Initialize model communicator, its size and the process rank
        # Here the same as for MPI_COMM_WORLD
        self.comm_ens:MPI.Comm
        self.npes_ens:int
        self.mype_ens:int
        self.comm_ens, self.npes_ens, self.mype_ens = self.init_parallel()

        self.is_task_consistent()
        self.is_CPU_consistent()

        # Initialize communicators for ensemble evaluations
        if self.mype_ens == 0:
            log.logger.info('Initialize communicators for assimilation with PDAF')
        
        # get the number of processors used by each model task
        self.local_npes_model:np.ndarray = self.get_processor_per_model()
        self.task_id:int = self.get_task_id()
        # get model communicator
        self.comm_model:MPI.Comm
        self.npes_model:int
        self.mype_model:int
        self.comm_model, self.npes_model, self.mype_model = self.get_model_communicator()

        # local ensemble member of each model task/processors.
        # When dim_ens > n_modeltasks, some processors/model tasks 
        # run dim_ens_l number of members sequentially. 
        # When n_modeltasks = dim_ens, dim_ens_l = 1
        # Here dim_ens_l is the ensemble size for current model task
        # all_dim_ens_l is the ensemble size for all task ids.
        self.dim_ens_l:int
        self.all_dim_ens_l:np.ndarray
        self.dim_ens_l, self.all_dim_ens_l = self.get_dim_ens_l()

        log.logger.info(f'MODEL: mype(w)= {self.mype_ens};'
                     f'; model task: {self.task_id}'
                     f'; mype(m)= {self.mype_model}'
                     f'; npes(m)= {self.npes_model}')

        # Generate communicator for filter
        self.filter_pe:bool
        self.comm_filter:MPI.Comm
        self.npes_filter:int
        self.mype_filter:int
        self.filter_pe, self.comm_filter = self.get_filter_communicator()
        self.npes_filter, self.mype_filter = self.get_filter_communicator_size_rank()

        # Generate communicators for communication
        self.comm_couple:MPI.Comm = self.get_couple_communicator()

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

        return MPI.COMM_WORLD, MPI.COMM_WORLD.Get_size(), MPI.COMM_WORLD.Get_rank()

    def get_processor_per_model(self) -> np.ndarray:
        """Store # PEs per ensemble
        used for info on PE 0 and for generation
        of model communicators on other Pes
        """

        local_npes_model:np.ndarray = np.zeros(self.n_modeltasks, dtype=int)

        local_npes_model[:] = np.floor(
            self.npes_ens/self.n_modeltasks)

        size:int = self.npes_ens \
            - self.n_modeltasks * local_npes_model[0]
        local_npes_model[:size] = local_npes_model[:size] + 1

        return local_npes_model

    def get_task_id(self) -> int:
        """get model communicator
        Generate communicators for model runs
        (Split COMM_ENSEMBLE)
        """
        pe_index:np.ndarray = np.cumsum(self.local_npes_model, dtype=int)
        # task id of the current processor
        task_id:int = np.arange(self.n_modeltasks)[pe_index <= self.mype_ens + self.local_npes_model[0]][-1] + 1
        
        return task_id

    def get_model_communicator(self) -> tuple[MPI.Comm, int, int]:
        """get model PE rank and size
        """
        comm_model = MPI.COMM_WORLD.Split(self.task_id, self.mype_ens)
        return comm_model, comm_model.Get_size(), comm_model.Get_rank()

    def get_dim_ens_l(self) -> tuple[int, np.ndarray]:
        """obtain number of ensemble members for each model task

        This means that each model task has to run dim_ens_l
        number of ensemble members sequentially
        """
        # number of ensemble for each model task
        dim_ens_l:int = self.dim_ens//self.n_modeltasks
        residual:int = self.dim_ens - dim_ens_l*self.n_modeltasks
        # number of ensmeble members across all PEs
        # e.g. self.all_dim_ens_l[0] is the ensemble size on the first PE
        # and self.all_dim_ens_l[-1] is the ensemble size on the last PE
        all_dim_ens_l:np.ndarray = dim_ens_l*np.ones(self.n_modeltasks, dtype=int)
        all_dim_ens_l[:residual] += 1
        # number of tasks on local PE
        dim_ens_l = all_dim_ens_l[self.task_id - 1]
        if self.mype_ens == 0:
            log.logger.debug (f'number of Ens per PE {all_dim_ens_l}')
        return dim_ens_l, all_dim_ens_l

    def get_filter_communicator(self) -> tuple[bool, MPI.Comm]:
        """Generate communicator for filter
        """
        # filter is only conducted in the first model task
        filterpe:bool = True if self.task_id == 1 else False
        my_color:int = self.task_id if filterpe else MPI.UNDEFINED
        comm_filter:MPI.Comm = MPI.COMM_WORLD.Split(my_color, self.mype_ens)
        return filterpe, comm_filter

    def get_filter_communicator_size_rank(self) -> tuple[int, int]:
        """get filter PE rank and size which should be same as model size and rank
        """
        return self.comm_model.Get_size(), self.comm_model.Get_rank()

    def get_couple_communicator(self) -> MPI.Comm:
        """Generate communicator for ensemble communications
        """
        return MPI.COMM_WORLD.Split(self.mype_model, self.mype_ens)

    def is_CPU_consistent(self) -> None:
        """Check consistency of number of parallel ensemble tasks
        """
        if self.n_modeltasks > self.npes_ens:
            # number of parallel tasks is set larger than available PEs ***
            self.n_modeltasks = self.npes_ens
            if self.mype_ens == 0:
                log.logger.warning('!!! Resetting number of parallel ensemble'
                                ' tasks to total number of PEs!')

    def is_task_consistent(self) -> None:
        """Check consistency of number of model tasks
        """
        assert self.dim_ens > 0, 'dim_ens (ensemble size) must be > 0'

        # Check consistency with ensemble size
        if self.n_modeltasks > self.dim_ens:
            # parallel ensemble tasks is set larger than ensemble size
            self.n_modeltasks = self.dim_ens

            if self.mype_ens == 0:
                log.logger.warning('!!! Resetting number of parallel'
                       'ensemble tasks to number of ensemble states!')

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
        if (self.task_id == 1):
            log.logger.info(f'{self.mype_ens},      {self.mype_filter},'
                   f'    {self.task_id}      {self.mype_model},'
                   f'   {color_couple},         {mype_couple},       {self.filter_pe}')
        MPI.COMM_WORLD.Barrier()
        if (self.task_id > 1):
            log.logger.info(f'{self.mype_ens},      {self.mype_filter},'
                   f'    {self.task_id}      {self.mype_model},'
                   f'   {color_couple},         {mype_couple},       {self.filter_pe}')
        MPI.COMM_WORLD.Barrier()
        if (self.mype_ens == 0):
            log.logger.info('')

    def finalize_parallel(self) -> None:
        """Finalize MPI
        """
        MPI.COMM_WORLD.Barrier()
        MPI.Finalize()
