import numpy as np


def step(model, pe, istep):
    if pe.task_id==1 and pe.mype_model==0:
    	print('step', istep)

    model.field_p = np.roll(model.field_p, 1, axis=0)

    field = np.zeros(model.nx, order='F')


    pe.COMM_model.Gather(model.field_p, field, 0)


    if pe.task_id==1 and pe.mype_model==0:
        np.savetxt(f'true_step{istep}.txt', field)