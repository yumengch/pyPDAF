import numpy as np
import copy

import model.shift

class Model:
    def __init__(self, nx, nt, pe):
        # model size
        self.nx = list(nx)
        # model size for each CPU
        self.get_nxp(pe)
        # model time steps
        self.total_steps = nt

    def get_nxp(self, pe):
        self.nx_p = copy.copy(self.nx)

        try:
            assert self.nx[-1]%pe.npes_model == 0
            self.nx_p[-1] = self.nx[-1]//pe.npes_model
        except AssertionError:
            print((f'...ERROR: Invalid number of'
                     f'processes: {pe.npes_model}...'))
            pe.abort_parallel()
            
    def init_field(self, filename, mype_model):
        # model field
        self.field_p = np.zeros(self.nx_p, order='F')
        offset = self.nx_p[-1]*mype_model
        self.field_p[:] = np.loadtxt(
                                    filename
                                    )[:, offset:self.nx_p[-1] + offset]

    def step(self, pe, istep):
        model.shift.step(self, pe, istep)

    def printInfo(self, USE_PDAF, pe):
        do_print = USE_PDAF and pe.mype_model == 0
        do_print = do_print or \
            (pe.task_id == 1 and pe.mype_model == 0 and not USE_PDAF) 
        if do_print:
            print('MODEL-side: INITIALIZE PARALLELIZED Shifting model' 
                ' MODEL')
            print('Grid size:', self.nx)
            print('Time steps', self.total_steps)
            print('-- Domain decomposition over', 
                pe.npes_model, ' PEs')
            print('-- local domain sizes (nx_p):', self.nx_p)
              