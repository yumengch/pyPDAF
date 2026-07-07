"""MPI/PDAF parallel initialisation helpers for the offline workflow."""
import log

from mpi4py import MPI
import numpy as np
import pyPDAF



class Parallel:
    """Store the MPI/PDAF communicator information used by this offline run.

    PDAF owns the filter communicator returned by ``pyPDAF.init_parallel``.
    For this offline example there is one model task group and the ensemble
    size is passed directly from ``[pdaf].dim_ens``.

    Attributes
    ----------
    comm_filter : MPI.Comm
        The MPI communicator for the filter
    mype_filter : int
        The rank of the current process in the filter communicator
    npes_filter : int
        The number of processes in the filter communicator
    """
    def __init__(self, dim_ens:int, screen:int) -> None:
        """Init the parallization required by PDAF

        Parameters
        ----------
        dim_ens : int
            number of ensemble members
        screen : int
            verbosity of the PDAF screen output
        """
        self.dim_ens = dim_ens
        result = pyPDAF.init_parallel(screen, 0, 0, dim_ens, 1, MPI.COMM_WORLD.py2f())
        pyPDAF.flush_fortran_stdout()
        self.comm_filter = MPI.Comm.f2py(result[4])
        self.mype_filter = result[5]
        self.npes_filter = result[6]

    def finalize_parallel(self) -> None:
        """Finalize MPI
        """
        MPI.COMM_WORLD.Barrier()
        MPI.Finalize()
