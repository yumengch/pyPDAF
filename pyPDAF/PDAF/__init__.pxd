cdef extern void c__pdaf_init(int* filtertype, int* subtype,
     int filter_param_i[], int* dim_pint, 
     double filter_param_r[], int* dim_preal,
     int* comm_model, int* comm_filter, int* comm_couple,
     int* task_id, int* n_modeltasks, bint* filterpe, 
     void (*c__init_ens_pdaf) (int*, int*, int*, double*, double*, 
                               double*, int*),
     int* screen, int* status_pdaf);
cdef extern void c__pdaf_get_state(int* steps, double* timenow, 
        int* doexit, 
        void (*c__next_observation_pdaf)(int*, int*, int*, double*),
        void (*c__distribute_state_pdaf)(int*, double*),
        void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,  
                                     double*, double*, double*, int*),
        int* status_pdaf);
cdef extern void c__pdaf_get_localfilter(int* localfilter);

cdef extern void c__pdafomi_assimilate_local(
     void (*c__collect_state_pdaf)(int*, double*),
     void (*c__distribute_state_pdaf)(int*, double*),
     void (*c__init_dim_obs_PDAFomi)(int*, int*),
     void (*c__obs_op_PDAFomi)(int*, int*, int*, 
                               double*, double*),
     void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                     double*, double*, double*, int*),
     void (*c__init_n_domains_pdaf)(int*, int*),
     void (*c__init_dim_l_pdaf)(int*, int*, int*),
     void (*c__init_dim_obs_l_PDAFomi)(int*, int*, int*, int*),
     void (*c__g2l_state_pdaf)(int*, int*, int*, double*,int*,double*),
     void (*c__l2g_state_pdaf)(int*, int*, int*, double*, 
                               int*, double*),
     void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                                          int* status_pdaf);
cdef extern void c__pdafomi_assimilate_global(
     void (*c__collect_state_pdaf)(int*, double*),
     void (*c__distribute_state_pdaf)(int*, double*),
     void (*c__init_dim_obs_PDAFomi)(int*, int*),
     void (*c__obs_op_PDAFomi)(int*, int*, int*, 
                               double*, double*),
     void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                     double*, double*, double*, int*),
     void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                                          int* status_pdaf);
cdef extern void c__pdafomi_assimilate_lenkf(
     void (*c__collect_state_pdaf)(int*, double*),
     void (*c__distribute_state_pdaf)(int*, double*),
     void (*c__init_dim_obs_PDAFomi)(int*, int*),
     void (*c__obs_op_PDAFomi)(int*, int*, int*, 
                               double*, double*),
     void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                     double*, double*, double*, int*),
     void (*c__localize_covar_PDAFomi)(int*, int*, double*, double*),
     void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                                          int* status_pdaf);