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

This is a translation of the shiting model in PDAF tutorials
"""
import logging
import numpy as np

from mpi4py import MPI

def step(field_p:np.ndarray, pe, step:int, USE_PDAF:bool) -> np.ndarray:
    """shifting model forward 'integration'

    Parameters
    ----------
    model : `Model.Model`
        Model object from example
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
        field = field.reshape(ny, nx, order='F')
        np.savetxt(f'true_step{step}.txt', field)

    return field_p
