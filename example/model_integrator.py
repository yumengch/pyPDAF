import logging
from mpi4py import MPI
import numpy as np

import config
import model
import PDAF_system

class model_integrator:
    def __init__(self, model_ens: list[model.model]) -> None:
        self.model_ens: list[model.model] = model_ens

    def forward(self, nsteps:int, PDAF_system:PDAF_system.PDAF_system) -> None:
        # When each model task runs one ensemble member,
        # i.e. no need to run each ensemble member sequentially,
        # we call this full parallel implementation
        if PDAF_system.pe.dim_ens_l == 1:
            self.forward_full_parallel(nsteps, PDAF_system)
        else:
            self.forward_flexible(nsteps, PDAF_system)

    def forward_full_parallel(self, nsteps:int, PDAF_system:PDAF_system.PDAF_system) -> None:
        for i in range(nsteps):
            self.model_ens[0].field_p = self.step(self.model_ens[0].field_p, PDAF_system.pe, i+1, config.USE_PDAF)
            if config.USE_PDAF:
                PDAF_system.assimilate_full_parallel()

    def forward_flexible(self, nsteps:int, PDAF_system:PDAF_system.PDAF_system) -> None:
        # create directory to
        current_step = 0
        # full DA system integration loop
        while current_step < nsteps:
            if config.USE_PDAF:
                # model integration
                for _ in range(PDAF_system.steps_for):
                    for i in range(PDAF_system.pe.dim_ens_l):
                        self.model_ens[i].field_p = self.step(self.model_ens[i].field_p, PDAF_system.pe, current_step + 1, config.USE_PDAF)
                    current_step += 1

                PDAF_system.assimilate_flexible()
            else:
                for i in range(PDAF_system.pe.dim_ens_l):
                    self.model_ens[i].field_p = self.step(self.model_ens[i].field_p, PDAF_system.pe, current_step + 1, config.USE_PDAF)
            current_step += 1

    def step(self, field_p:np.ndarray, pe, step:int, USE_PDAF:bool) -> np.ndarray:
        """shifting model forward 'integration'

        Parameters
        ----------
        field_p : `np.ndarray`
            model field
        pe : `parallelization.parallelization`
            parallelization object from example
        step : int
            current time step
        USE_PDAF : bool
            whether PDAF is performed
        """
        if pe.task_id == 1 and pe.mype_model == 0:
            logging.info(f'model step: {step}')

        field_p = np.roll(field_p, 1, axis=0)

        if USE_PDAF:
            return field_p

        # collect the length of state vector on each processor
        ny, _ = field_p.shape
        n_gp = len(field_p.ravel()) 
        send_counts = np.array(pe.comm_model.gather(n_gp, root=0))
        if pe.mype_model == 0:
            # get the length of the full state vector
            n_gp_total = np.sum(send_counts)
            # get nx as domain decomposition only occurs at nx direction
            nx = n_gp_total//ny
            # declare the full ensemble
            field = np.zeros(n_gp_total)
            # displacement of each of the full ensemble
            displacements = np.insert(np.cumsum(send_counts), 0, 0)[0:-1]
        else:
            displacements = None
            field = None

        pe.comm_model.Gatherv([field_p.ravel(order='F'), MPI.DOUBLE],
                            [field,
                            send_counts, displacements, MPI.DOUBLE],
                            root=0)

        if pe.task_id == 1 and pe.mype_model == 0:
            assert field is not None, '...ERROR: field is None'
            field = field.reshape(ny, nx, order='F')
            np.savetxt(f'true_step{step}.txt', field)

        return field_p