import log

from mpi4py import MPI
import numpy as np

import pyPDAF


class Parallelisation:
    def __init__(self, dim_ens:int, screen:int) -> None:
        """Init the parallization required by PDAF

        Parameters
        ----------
        dim_ens : int
            number of ensemble members
        screen : int
            verbosity of the PDAF screen output
        """
        result = pyPDAF.init_parallel(screen, 1, 0, dim_ens, 1, MPI.COMM_WORLD)

    def finalize_parallel(self) -> None:
        """Finalize MPI
        """
        MPI.COMM_WORLD.Barrier()
        MPI.Finalize()
