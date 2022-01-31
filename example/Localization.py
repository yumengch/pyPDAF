import numpy as np
import pyPDAF.PDAF.PDAFomi as PDAFomi


class Localization:
    def __init__(self, loc_weight, local_range, srange):
        self.loc_weight = loc_weight
        self.local_range = local_range
        self.srange = srange
    
    def set_lim_coords(self, nx_p, pe):
        # Get offset of local domain in global domain in x-direction
        off_nx = nx_p[-1]*pe.mype_filter

        lim_coords = np.zeros((2, 2))
        lim_coords[0,0] = float(off_nx + 1)
        lim_coords[0,1] = float(off_nx + nx_p[-1])
        lim_coords[1] = 0

        PDAFomi.set_domain_limits(lim_coords)

    def init_dim_l_pdaf(self, nx_p, mype_filter, step, domain_p, dim_l):
        # initialize local state dimension
        dim_l = 1
        # initialize coordinates of local domain
        # we use grid point indices as coordinates,
        #  but could e.g. use meters
        self.coords_l = np.zeros(2)
        offset = mype_filter*np.prod(nx_p)
        self.coords_l[0] = domain_p + offset
        self.coords_l[0] = np.ceil(self.coords_l[0]/nx_p[0])
        self.coords_l[1] = domain_p + offset
        self.coords_l[1] = self.coords_l[1] \
                            - (self.coords_l[0] - 1)*nx_p[0]
        # initialize array of indices of the local state
        #  vector elements in the global state vector.

        # allocate array
        self.id_lstate_in_pstate = np.zeros(dim_l)

        # here the local domain is a single grid point 
        # and variable given by domain_p
        self.id_lstate_in_pstate[0] = domain_p

        return dim_l

    def init_n_domains_pdaf(self, assim_dim, step):
        return assim_dim.dim_state_p

    def g2l_state_pdaf(self, step, domain_p, state_p, state_l):
        # generic initialization 
        # using id_lstate_in_pstate set in init_dim_l_pdaf
        state_l[:] = state_p[self.id_lstate_in_pstate]

    def l2g_state_pdaf(self, step, domain_p, state_l, state_p):
        state_p[self.id_lstate_in_pstate] = state_l[:]