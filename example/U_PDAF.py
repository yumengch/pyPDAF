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

Attributes
----------
firsttime : bool
    global variable for prepoststep_ens_pdaf
"""
import numpy as np

import U_PDAFomi


def init_ens_pdaf(model, pe, assim_dim,
                  filtertype, state_p,
                  uinv, ens_p, status_pdaf):
    """user-defined init_ens_pdaf function

    Parameters
    ----------
    model : `Model.Model`
        model object
    pe : `parallelization.parallelization`
        parallelization object
    assim_dim : `AssimilationDimensions.AssimilationDimensions`
        an object of AssimilationDimensions
    filtertype : int
        type of filter
    state_p : ndarray
        1D state vector on local PE
    uinv : ndarray
        2D left eigenvector with shape (dim_ens - 1,dim_ens - 1)
    ens_p : ndarray
        ensemble state vector on local PE (dim_ens, dim_p)
    status_pdaf : int
        status of PDAF

    Returns
    -------
    status_pdaf : int
        status of PDAF
    """
    dim_ens, dim_p = ens_p.shape
    filename = 'inputs_online/ens_{i}.txt'
    off_nx = model.nx_p[-1]*pe.mype_model
    for i_ens in range(dim_ens):
        field_p = np.loadtxt(
                        filename.format(i=i_ens+1)
                            )[:, off_nx:off_nx+model.nx_p[-1]]
        ens_p[i_ens] = field_p.reshape(assim_dim.dim_state_p, order='F')
    return status_pdaf


def collect_state_pdaf(model, assim_dim, state_p):
    """Collect state vector in PDAF from model

    rely on Python's pass by reference

    Parameters
    ----------
    model : `Model.Model`
        model object
    assim_dim : `AssimilationDimensions.AssimilationDimensions`
        an object of AssimilationDimensions
    state_p : ndarray
        1D state vector on local PE
    """
    state_p[:] = model.field_p.reshape(assim_dim.dim_state_p, order='F')


def distribute_state_pdaf(model, state_p):
    """Distribute state vector to model field

    rely on Python's pass by reference

    Parameters
    ----------
    model : `Model.Model`
        model object
    state_p : ndarray
        1D state vector on local PE
    """
    model.field_p[:] = state_p.reshape(*model.nx_p, order='F')


def next_observation_pdaf(model, pe, delt_obs,
                          stepnow, nsteps, doexit, time):
    """The time for the next observation

    Parameters
    ----------
    model : `Model.Model`
        model object
    pe : `parallelization.parallelization`
        parallelization object
    delt_obs : int
        frequency of observations
    stepnow : int
        Current time step
    nsteps : int
        steps between assimilation
    doexit : int
        Whether exit PDAF assimilation
    time : double
        Current model time

    Returns
    -------
    nsteps : int
        steps between assimilation
    doexit : int
        Whether exit PDAF assimilation
    time : double
        Current model time
    """
    if (stepnow + nsteps <= model.total_steps):
        nsteps = delt_obs
        doexit = 0
        if (pe.mype_world == 0):
            print((stepnow, 'Next observation at time step',
                   stepnow + nsteps))
    else:
        nsteps = 0
        doexit = 1
        if (pe.mype_world == 0):
            print((stepnow, 'No more observations - end assimilation'))

    return nsteps, doexit, time


firsttime = True


def prepoststep_ens_pdaf(assim_dim, model, pe, obs,
                         step, state_p, uinv, ens_p):
    """pre- and post-processing of ensemble

    Parameters
    ----------
    assim_dim : `AssimilationDimensions.AssimilationDimensions`
        an object of AssimilationDimensions
    model : `Model.Model`
        model object
    pe : `parallelization.parallelization`
        parallelization object
    obs : `OBS.OBS`
        observation object
    step : int
        Current time step
    state_p : ndarray
        1D state vector on local PE
    uinv : ndarray
        2D left eigenvector with shape (dim_ens - 1,dim_ens - 1)
    ens_p : ndarray
        ensemble state vector on local PE (dim_ens, dim_p)
    """
    global firsttime
    if (firsttime):
        print('Analyze initial state ensemble')
        anastr = 'ini'
    else:
        if (step < 0):
            print('Analyze and write forecasted state ensemble')
            anastr = 'for'
        else:
            print('Analyze and write assimilated state ensemble')
            anastr = 'ana'

    dim_ens, dim_p = ens_p.shape
    variance = np.zeros(assim_dim.dim_state)

    state_p[:] = np.mean(ens_p, axis=0)[:]
    variance_p = np.var(ens_p, axis=0, ddof=1)
    if pe.mype_filter != 0:
        pe.COMM_filter.Send(variance_p, 0, pe.mype_filter)
    else:
        variance[:dim_p] = variance_p[:]
        for i in range(1, pe.npes_filter):
            pe.COMM_filter.Recv(
                              variance[i*dim_p:(i+1)*dim_p], i, i)

    rmserror_est = np.sqrt(np.sum(
                            variance
                                 )/assim_dim.dim_state)

    if pe.mype_filter == 0:
        print('RMS error: ', rmserror_est)

    if (not firsttime):
        ens = np.zeros((assim_dim.dim_ens, assim_dim.dim_state))
        state = np.zeros(assim_dim.dim_state)
        if pe.mype_filter != 0:
            pe.COMM_filter.Send(ens_p[:], 0, pe.mype_filter)
        else:
            ens[:, :dim_p] = ens_p
            ens_tmp = np.zeros(ens_p.shape)
            for i in range(1, pe.npes_filter):
                pe.COMM_filter.Recv(
                                ens_tmp, i, i)
                ens[:, i*dim_p:(i+1)*dim_p] = ens_tmp[:, :]
            print('--- write ensemble and state estimate')

            stepstr = step if step >= 0 else -step
            field = np.zeros(model.nx)
            for i in range(dim_ens):
                field = ens[i].reshape(*model.nx, order='F')
                filename = f'ens_{i}_step{stepstr}_{anastr}.txt'
                np.savetxt(filename, field, delimiter=';')

        if pe.mype_filter != 0:
            pe.COMM_filter.Send(state_p, 0, pe.mype_filter)
        else:
            state[:dim_p] = state_p[:]
            state_p_tmp = np.zeros(state_p.shape)
            for i in range(1, pe.npes_filter):
                pe.COMM_filter.Recv(
                            state_p_tmp, i, i)
                state[i*dim_p:(i+1)*dim_p] = state_p_tmp
            filename = f'state_step{stepstr}_{anastr}.txt'
            np.savetxt(filename,
                       state.reshape(*model.nx, order='F'), delimiter=';')

    U_PDAFomi.deallocate_obs_pdafomi(obs, step)

    firsttime = False
