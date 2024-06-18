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
import numpy as np


def step(model, pe, step, USE_PDAF):
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
        print(('step', step))

    model.field_p = np.roll(model.field_p, 1, axis=0)

    if USE_PDAF:
        return

    if pe.mype_model == 0:
        field_gather = np.empty((pe.npes_model,) + tuple(model.nx_p))
    else:
        field_gather = None

    pe.COMM_model.Gather(model.field_p, field_gather, 0)

    if pe.task_id == 1 and pe.mype_model == 0:
        field = np.zeros(model.nx)
        for i in range(pe.npes_model):
            offset = model.nx_p[-1]*i
            field[:, offset:offset+model.nx_p[-1]] = field_gather[i]
        np.savetxt(f'true_step{step}.txt', field)
