cdef extern void c__init_pdafomi(int* n_obs);
cdef extern void c__set_pdafomi_doassim(int* i_obs, int* doassim);
cdef extern void c__set_pdafomi_disttype(int* i_obs, int* disttype);
cdef extern void c__set_pdafomi_ncoord(int* i_obs, int* ncoord);
cdef extern void c__set_pdafomi_id_obs_p(int* i_obs, int* nrows, 
                                         int* dim_obs_p, 
                                         int* id_obs_p);
cdef extern void c__set_pdafomi_icoeff_p(int* i_obs, int* nrows, 
                                         int* dim_obs_p, 
                                         double* icoeff_p);
cdef extern void c__set_pdafomi_domainsize(int* i_obs, int* ncoord, 
                                           double* domainsize);
cdef extern void c__set_pdafomi_obs_err_type(int* i_obs, 
                                             int* obs_err_type);
cdef extern void c__set_pdafomi_use_global_obs(int* i_obs, 
                                               int* use_global_obs);
cdef extern void c__pdafomi_gather_obs(int* i_obs, int* nrows, 
                                       int* dim_obs_p, double* obs_p, 
                                       double* ivar_obs_p, 
                                       double* ocoord_p, 
                                       double* local_range, 
                                       int* dim_obs);
cdef extern void c__pdafomi_set_domain_limits(double* lim_coords);
cdef extern void c__pdafomi_obs_op_gridpoint(int* i_obs, int* dim_p, 
                                             int* dim_obs, 
                                             double* state_p, 
                                             double* ostate);
cdef extern void c__pdafomi_init_dim_obs_l(int* i_obs, 
                                           double* coords_l, 
                                           int* locweight, 
                                           double* local_range, 
                                           double* srange,
                                           int* dim_obs_l);
cdef extern void c__pdafomi_localize_covar(int* i_obs, int* dim_p,
                                           int* dim_obs,
                                           int* dim_coords,
                                           int* locweight, 
                                           double* local_range, 
                                           double* srange, 
                                           double* coords_p, 
                                           double* hp_p, double* hph);
cdef extern void c__pdafomi_deallocate_obs(int* i_obs, int* step);
