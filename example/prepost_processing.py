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
import os
import typing

from mpi4py import MPI
import numpy as np

import model
import parallelisation

class prepost:
    """User-supplied functions for pre and post processing of the ensemble.

    Attributes
    ----------
    model : `model.model`
        model object
    pe : `parallelisation.parallelisation`
        parallelisation object
    -------
    """
    def __init__(self, model_t: model.model, pe:parallelisation.parallelisation) -> None:
        self.model:model.model = model_t
        self.pe:parallelisation.parallelisation = pe
        os.makedirs('outputs', exist_ok=True)

    def get_full_ens(self, dim_p:int, dim_ens:int, ens_p:np.ndarray) -> typing.Union[np.ndarray, None]:
        """Gather total ensemble from each local processors
        """
        if self.pe.npes_filter == 1: return ens_p
        # get total dim

        ## collect full ensemble from domain decomposed ensemble
        # collect the length of state vector on each processor (local domain)
        all_dim_p:np.ndarray = np.array(self.pe.comm_filter.gather(dim_p, root=0))

        displacements : np.ndarray | None
        send_counts : np.ndarray | None
        ens : np.ndarray | None
        if self.pe.mype_filter == 0:
            # number of elements of the array on each processor
            send_counts = all_dim_p*dim_ens
            # get the length of the full state vector
            dim:int = np.sum(all_dim_p)
            # declare the full ensemble
            ens = np.zeros(dim*dim_ens)
            # displacement of each of the full ensemble
            displacements = np.insert(np.cumsum(send_counts), 0, 0)[0:-1]
        else:
            displacements = None
            ens = None
            send_counts = None

        # using row-major C order to ensure that
        # MPI gathers a continuous row-major array of the model domain
        # that is dim_ens number of first element of the state vector
        # followed by dim_ens number of the second element of the state vector, etc.
        ens_p_send = ens_p.ravel()
        self.pe.comm_filter.Gatherv([ens_p_send, MPI.DOUBLE],
                               [ens,
                                send_counts, displacements, MPI.DOUBLE],
                                root=0)
        if self.pe.mype_filter == 0:
            assert isinstance(ens, np.ndarray), 'ens should be a numpy array'
            ens = ens.reshape(dim, dim_ens)
            # As a consequence of domain decomposition in nx instead of ny
            # (following the PDAF tutorial)
            # we need to reorder the array after merging from different processors
            displ = np.insert(np.cumsum(all_dim_p), 0, 0)[1:]
            ens_tmp = ens[:displ[0]].reshape(self.model.ny, self.model.nx_p, dim_ens)
            if len(displ) > 0:
                for c0, c1 in zip(displ[:-1], displ[1:]):
                    ens_tmp = np.concatenate([ens_tmp,
                            ens[c0:c1].reshape(self.model.ny, self.model.nx_p, dim_ens)], axis=1)
            ens = ens_tmp.reshape(dim, dim_ens)

        return ens

    def initial_process(self, step:int, dim_p:int, dim_ens:int, dim_ens_p:int,
                        dim_obs_p:int, state_p:np.ndarray, uinv:np.ndarray,
                        ens_p:np.ndarray, flag:int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """initial processing of the ensemble before it is distributed to model fields
        """
        ens = self.get_full_ens(dim_p, dim_ens, ens_p)
        if self.pe.mype_filter == 0:
            assert isinstance(ens, np.ndarray), 'ens should be a numpy array'
            log.logger.info (f'RMS error according to sampled variance: {np.sqrt(np.mean(np.var(ens, axis=1, ddof=1)))}')
        return state_p, uinv, ens_p

    def preprocess(self, step:int, dim_p:int, dim_ens:int, ens_p:np.ndarray) -> None:
        """preprocessing of the ensemble before it is used by DA algorithms
        """
        ens = self.get_full_ens(dim_p, dim_ens, ens_p)
        if self.pe.mype_filter == 0:
            assert isinstance(ens, np.ndarray), 'ens should be a numpy array'
            log.logger.info (f'Forecast RMS error according to sampled variance: {np.sqrt(np.mean(np.var(ens, axis=1, ddof=1)))}')
            os.makedirs('outputs', exist_ok=True)
            for i in range(dim_ens):
                np.savetxt(os.path.join('outputs', f'ens_{i+1}_step{-step}_for.txt') , ens[:, i].reshape(self.model.ny, self.model.nx) )

    def postprocess(self, step:int, dim_p:int, dim_ens:int, ens_p:np.ndarray) -> None:
        """initial processing of the ensemble before it is distributed to model fields
        """
        ens = self.get_full_ens(dim_p, dim_ens, ens_p)
        if self.pe.mype_filter == 0:
            assert isinstance(ens, np.ndarray), 'ens should be a numpy array'
            log.logger.info (f'Analysis RMS error according to sampled variance: {np.sqrt(np.mean(np.var(ens, axis=1, ddof=1)))}')
            os.makedirs('outputs', exist_ok=True)
            for i in range(dim_ens):
                np.savetxt(os.path.join('outputs', f'ens_{i+1}_step{step}_ana.txt') , ens[:, i].reshape(self.model.ny, self.model.nx) )

    def prepostprocess(self, step:int, dim_p:int, dim_ens:int, dim_ens_p:int,
                       dim_obs_p:int, state_p:np.ndarray, uinv:np.ndarray,
                       ens_p:np.ndarray, flag:int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """pre-/post-processing of the ensemble as user-supplied functions
        """
        if step < 0:
          self.preprocess(step, dim_p, dim_ens, ens_p)
        else:
          self.postprocess(step, dim_p, dim_ens, ens_p)
        return state_p, uinv, ens_p
