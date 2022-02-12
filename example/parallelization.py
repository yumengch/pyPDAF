from mpi4py import MPI
import numpy as np


class parallelization:
    def __init__(self, dim_ens, n_modeltasks, screen):
        """Init the parallization required by PDAF

        Parameters
        ----------
        n_modeltasks: int
            Number of model tasks/ ensemble size
            This parameter should be the same
            as the number of PEs
        screen: int
            The verbose level of screen output.
            screen = 3 is the most verbose level.
        """

        self.n_modeltasks = n_modeltasks
        self.n_filterpes = 1

        self.init_parallel()
        self.getEnsembleSize()

        # Initialize communicators for ensemble evaluations
        if (self.mype_world == 0):
            print(('Initialize communicators for assimilation with PDAF'))

        self.isCPUConsistent()
        self.isTaskConsistent(n_modeltasks)

        # get ensemble communicator
        self.getModelComm()
        self.getModelPERankSize()

        if (screen > 1):
            print(('MODEL: mype(w)= ', self.mype_world, 
                 '; model task: ', self.task_id,
                 '; mype(m)= ', self.mype_model, 
                 '; npes(m)= ', self.npes_model))

        # Generate communicator for filter
        self.getFilterComm()
        self.getFilterPERankSize()

        # Generate communicators for communication
        self.getCoupleComm()

        self.printInfo(screen)

    def printInfo(self, screen):
        # *** local variables ***
        #  Rank and size in COMM_couple
        mype_couple = self.COMM_couple.Get_rank()
        #  Variables for communicator-splitting 
        color_couple = self.mype_model + 1

        if (screen > 0):
            if (self.mype_world == 0):
                print(('PE configuration:'))
                print(('world', 'filter', '   model    ', '   couple   ', 'filterPE'))
                print(('rank ', ' rank ', 'task', 'rank', 'task', 'rank', 'T/F'))
                print(('--------------------------------------'))
            MPI.COMM_WORLD.Barrier()
            if (self.task_id == 1):
                print((self.mype_world, self.mype_filter, self.task_id, 
                        self.mype_model, color_couple,
                        mype_couple, self.filterpe))
            if (self.task_id > 1):
                print((self.mype_world, ' ', self.task_id, self.mype_model, 
                    color_couple, mype_couple, self.filterpe))
            MPI.COMM_WORLD.Barrier()
            if (self.mype_world == 0):
                print('')

    def init_parallel(self):
        """!> Initialize MPI
       
        Routine to initialize MPI, the number of PEs
        (npes_world) and the rank of a PE (mype_world).
        The model is executed within the scope of the
        communicator Comm_model. It is also initialized
        here together with its size (npes_model) and 
        the rank of a PE (mype_model) within Comm_model."""

        if not MPI.Is_initialized():
            MPI.Init()

        # Initialize model communicator, its size and the process rank
        # Here the same as for MPI_COMM_WORLD
        self.COMM_model = MPI.COMM_WORLD
        self.npes_model = None
        self.npes_world = None
        self.mype_model = None
        self.mype_world = None

        self.npes_world = self.COMM_model.Get_size()
        self.mype_world = self.COMM_model.Get_rank()

    def getModelComm(self):
        """COMM_MODEL
          Generate communicators for model runs
          (Split COMM_ENSEMBLE)"""
        self.getPEperModel()

        pe_index = np.cumsum(self.local_npes_model, dtype=int)

        mype_ens = self.mype_world

        if mype_ens + 1 <= self.local_npes_model[0]:
            self.task_id = 1
        else:
            self.task_id = np.where(pe_index < mype_ens + 1)[0][-1] + 2
        
        self.COMM_model = MPI.COMM_WORLD.Split(self.task_id, mype_ens)

    def getFilterComm(self):
        """Generate communicator for filter"""
        # Init flag FILTERPE (all PEs of model task 1)
        self.filterpe = True if self.task_id == 1 else False

        my_color = self.task_id if self.filterpe \
                                else MPI.UNDEFINED

        self.COMM_filter = MPI.COMM_WORLD.Split(my_color, 
                                                self.mype_world)

    def getCoupleComm(self):
        """Generate communicator for filter"""
        # Init flag FILTERPE (all PEs of model task 1)
        color_couple = self.mype_model + 1

        self.COMM_couple = MPI.COMM_WORLD.Split(color_couple, 
                                                self.mype_world)

    def getModelPERankSize(self):
        self.npes_model = self.COMM_model.Get_size()
        self.mype_model = self.COMM_model.Get_rank()

    def getFilterPERankSize(self):
        self.npes_filter = self.COMM_model.Get_size()
        self.mype_filter = self.COMM_model.Get_rank()

    def getEnsembleSize(self):
        """Parse number of model tasks

           The module variable is N_MODELTASKS. 
           Since it has to be equal to the ensemble size 
           we parse dim_ens from the command line."""
        #  handle for command line parser
        # handle = 'dim_ens'
        # parse(handle, self.n_modeltasks)
        # self.n_modeltasks = 1
        pass

    def getPEperModel(self):
        """Store # PEs per ensemble 
           used for info on PE 0 and for generation
           of model communicators on other Pes"""

        self.local_npes_model = np.zeros(self.n_modeltasks, dtype=int)

        self.local_npes_model[:] = np.floor(
            self.npes_world/self.n_modeltasks)

        size = self.npes_world \
                - self.n_modeltasks * self.local_npes_model[0]
        self.local_npes_model[:size] = self.local_npes_model[:size] + 1

    def isCPUConsistent(self):
        """Check consistency of number of parallel ensemble tasks"""
        pass
        # try:
        #     assert self.n_modeltasks == self.npes_world
        #     # parallel tasks is not identical to available PEs
        #     # This needs to hold for a model without parallelization
        # except AssertionError:
        #     if (self.mype_world == 0):
        #         print((f'!!! ERROR: dim_ens {self.n_modeltasks} '
        #                 f'must equal number of processes {self.npes_world}!'))
        #     self.abort_parallel()

    def isTaskConsistent(self, dim_ens):
        """Check consistency of number of model tasks"""

        # For dim_ens=0 no consistency check
        # for the ensemble size with the 
        # number of model tasks is performed.
        if (dim_ens <= 0): return

        # Check consistency with ensemble size
        if (self.n_modeltasks > dim_ens):
            # parallel ensemble tasks is set larger than ensemble size
            self.n_modeltasks = dim_ens

            if (self.mype_world == 0):
                print(('!!! Resetting number of parallel'
                    'ensemble tasks to number of ensemble states!'))

    def finalize_parallel(self):
        """ Finalize MPI"""
        MPI.COMM_WORLD.Barrier()
        MPI.Finalize()

    def abort_parallel(self):
        """Abort MPI"""
        MPI.COMM_WORLD.Abort(1)