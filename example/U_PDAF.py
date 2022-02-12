import numpy as np

import U_PDAFomi


firsttime = True

def init_ens_pdaf(model, pe, assim_dim, 
                        filtertype, state_p, 
                        uinv, ens_p, status_pdaf):
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
    state_p[:] = model.field_p.reshape(assim_dim.dim_state_p, order='F')

def distribute_state_pdaf(model, state_p):
    model.field_p[:] = state_p.reshape(*model.nx_p, order='F')

def next_observation_pdaf(model, pe, delt_obs, 
                            stepnow, nsteps, doexit, time):
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

def prepoststep_ens_pdaf(assim_dim, model, pe, obs,
                            step, state_p, uinv, ens_p):
    global firsttime
    print(('prepoststep_ens_pdaf:', firsttime))
    if (firsttime):
        print( 'Analyze initial state ensemble')
        anastr = 'ini'
    else:
        if (step<0):
            print('Analyze and write forecasted state ensemble')
            anastr = 'for'
        else:
            print('Analyze and write assimilated state ensemble')
            anastr = 'ana'

    dim_ens, dim_p = ens_p.shape
    variance = np.zeros(assim_dim.dim_state)

    state_p[:] = np.mean(ens_p, axis = 0)[:]
    variance_p = np.var(ens_p, axis = 0, ddof=1)
    if pe.mype_filter != 0:
        pe.COMM_filter.Send(variance_p, 0, 
                                    pe.mype_filter)
    else:
        variance[:dim_p] = variance_p[:]
        for i in range(1, pe.npes_filter):
            pe.COMM_filter.Recv(
                              variance[i*dim_p:(i+1)*dim_p], i, i)

    truth = np.random.random(dim_ens)
    rmserror_est = np.sqrt(np.sum(
                            variance
                                 )/assim_dim.dim_state)

    if pe.mype_filter == 0: print('RMS error: ', rmserror_est)

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
                print(np.isfortran(ens))
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
            state_p_tmp = np.zeros(state_p.shape, order='F')
            for i in range(1, pe.npes_filter):
                pe.COMM_filter.Recv(
                            state_p_tmp, i, i)
                state[i*dim_p:(i+1)*dim_p] = state_p_tmp
            filename = f'state_step{stepstr}_{anastr}.txt'
            np.savetxt(filename, state.reshape(*model.nx, order='F')
                                                , delimiter=';')

    U_PDAFomi.deallocate_obs_pdafomi(obs, step)

    firsttime = False