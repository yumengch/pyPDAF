import numpy as np


def step(model, pe, istep, USE_PDAF):
    if pe.task_id==1 and pe.mype_model==0:
    	print('step', istep)

    model.field_p = np.roll(model.field_p, 1, axis=0)

    if USE_PDAF: return

    if pe.mype_model == 0:
        field_gather = np.empty((pe.npes_model,)+ tuple(model.nx_p))
    else:
        field_gather = None

    pe.COMM_model.Gather(model.field_p, field_gather, 0)

    if pe.task_id==1 and pe.mype_model==0:
        field = np.zeros(model.nx)
        for i in range(pe.npes_model):
            offset = model.nx_p[-1]*i
            field[:, offset:offset+model.nx_p[-1]] = field_gather[i]
        np.savetxt(f'true_step{istep}.txt', field)