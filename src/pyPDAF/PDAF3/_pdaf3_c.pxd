cdef extern void c__pdaf3_init(int* filtertype, int* subtype, int* stepnull,
    int* param_int, int* dim_pint, double* param_real, int* dim_preal,
    void (*c__init_ens_pdaf)(int* , int* , int* , double* , double* ,
                             double* , int* ),
    int* in_screen, int* outflag) noexcept nogil;

cdef extern void c__pdaf3_init_forecast(
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ),
    void (*c__distribute_state_pdaf)(int* , double* ),
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* ,
                                double* , double* , double* , int* ),
    int* outflag) noexcept nogil;

cdef extern void c__pdaf3_init_parallel(int* screen, int* type_parallel,
    int* online_coupling, int* dim_ens, int* n_modeltasks, int* COMM_model,
    int* mype_model, int* npes_model, int* COMM_assim, int* mype_assim,
    int* npes_assim, int* task_id) noexcept nogil;


cdef extern void c__pdaf3_set_parallel(int* in_comm_pdaf,
    int* in_comm_model, int* in_comm_filter, int* in_comm_couple, int* in_task_id,
    int* in_n_modeltasks, bint* in_filterpe, int* flag) noexcept nogil;
