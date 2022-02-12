import numpy as np

def init(n_obs):
    cdef int c__n_obs = n_obs
    c__init_pdafomi(&c__n_obs)

def setOMIessential(int i_obs, int doassim, 
                    int disttype, int ncoord, id_obs_p):

    c__set_pdafomi_doassim(&i_obs, &doassim)
    c__set_pdafomi_disttype(&i_obs, &disttype)
    c__set_pdafomi_ncoord(&i_obs, &ncoord)
    shape = id_obs_p.shape
    cdef int nrows = shape[0]
    cdef int dim_obs_p = shape[1]
    cdef int[::1, ::] id_obs_p_view = np.array(id_obs_p, order='F', dtype=np.intc)
    c__set_pdafomi_id_obs_p(&i_obs, &nrows, &dim_obs_p, &id_obs_p_view[0][0])

def set_icoeff_p(int i_obs, icoeff_p):
    cdef int nrows = icoeff_p.shape[0]
    cdef int dim_obs_p = icoeff_p.shape[1]
    cdef double[::1, ::] icoeff_p_view = np.array(icoeff_p, order='F')
    c__set_pdafomi_icoeff_p(&i_obs, &nrows, &dim_obs_p, &icoeff_p_view[0][0])

def set_domainsize(int i_obs, domainsize):
    cdef double[::1] domainsize_view = domainsize
    cdef int ncoord = len(domainsize)
    c__set_pdafomi_domainsize(&i_obs, &ncoord, &domainsize_view[0])

def set_obs_err_type(int i_obs, int obs_err_type):
    c__set_pdafomi_obs_err_type(&i_obs, &obs_err_type)

def set_use_global_obs(int i_obs, int use_global_obs):
    c__set_pdafomi_use_global_obs(&i_obs, &use_global_obs)

def gather_obs(int i_obs, int dim_obs_p, int nrows, obs_p, ivar_obs_p, 
                ocoord_p, double local_range):
    cdef double[::1, ::] ocoord_p_view = np.array(ocoord_p, order='F')
    cdef double[::1] obs_p_view = obs_p
    cdef double[::1] ivar_obs_p_view = ivar_obs_p
    cdef int dim_obs
    c__pdafomi_gather_obs(&i_obs, &nrows, &dim_obs_p,
                           &obs_p_view[0], &ivar_obs_p_view[0], 
                           &ocoord_p_view[0][0], &local_range, &dim_obs)
    return dim_obs

def set_domain_limits(lim_coords):
    cdef double[::1, ::] lim_coords_view = np.array(
                                            lim_coords, order='F')
    c__pdafomi_set_domain_limits(&lim_coords_view[0][0])

def obs_op_gridpoint(int i_obs, state_p, ostate):
    cdef int dim_p, dim_obs
    dim_p = len(state_p) 
    dim_obs = len(ostate)
    cdef double[::1] state_p_view = state_p
    cdef double[::1] ostate_view = ostate
    c__pdafomi_obs_op_gridpoint(&i_obs, &dim_p, &dim_obs, 
                                &state_p_view[0], &ostate_view[0])

def init_dim_obs_l(int i_obs, coords_l, int loc_weight, 
                    double local_range, double srange):
    cdef double[::1, ::] coords_l_view = np.array(coords_l, order='F')
    cdef int dim_obs_l
    c__pdafomi_init_dim_obs_l(&i_obs, &coords_l_view[0][0], &loc_weight, 
                              &local_range, &srange, &dim_obs_l)
    return dim_obs_l

def localize_covar(int i_obs, int loc_weight, 
                   double local_range, double srange, 
                   coords_p, hp_p, hph):
    cdef double[::1, ::] coords_p_view = np.array(coords_p, order='F')
    cdef int dim_coords, dim_p, dim_obs
    dim_coords, dim_p = coords_p.shape
    dim_obs = hp_p.shape[0]
    cdef double[::1, ::] hp_p_view = np.array(hp_p, order='F')
    cdef double[::1, ::] hph_view = np.array(hph, order='F')
    c__pdafomi_localize_covar(&i_obs, &dim_p, &dim_obs, &dim_coords,
                              &loc_weight, &local_range, &srange, 
                              &coords_p_view[0][0], 
                              &hp_p_view[0][0], 
                              &hph_view[0][0]);

def deallocate_obs(int i_obs, int step):
    c__pdafomi_deallocate_obs(&i_obs, &step);
