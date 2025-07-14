from .cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_correlation_function(int* ctype, double* length, 
    double* distance, double* value) noexcept nogil;

cdef extern void c__pdaf_deallocate() noexcept nogil;

cdef extern void c__pdaf_eofcovar(int* dim, int* nstates, int* nfields, 
    int* dim_fields, int* offsets, int* remove_mstate, int* do_mv, 
    double* states, double* stddev, double* svals, double* svec, 
    double* meanstate, int* verbose, 
    int* status) noexcept nogil;

cdef extern void c__pdaf_force_analysis() noexcept nogil;

cdef extern void c__pdaf_gather_dim_obs_f(int* dim_obs_p, 
    int* dim_obs_f) noexcept nogil;

cdef extern void c__pdaf_gather_obs_f(double* obs_p, double* obs_f, 
    int* status) noexcept nogil;

cdef extern void c__pdaf_gather_obs_f2(double* coords_p, double* coords_f, 
    int* nrows, int* status) noexcept nogil;

cdef extern void c__pdaf_gather_obs_f_flex(int* dim_obs_p, int* dim_obs_f, 
    double* obs_p, double* obs_f, int* status) noexcept nogil;

cdef extern void c__pdaf_gather_obs_f2_flex(int* dim_obs_p, int* dim_obs_f, 
    double* coords_p, double* coords_f, int* nrows, 
    int* status) noexcept nogil;

cdef extern void c__pdaf_init(int* filtertype, int* subtype, int* stepnull, 
    int* param_int, int* dim_pint, double* param_real, int* dim_preal, 
    int* comm_model, int* comm_filter, int* comm_couple, int* task_id, 
    int* n_modeltasks, bint* in_filterpe, 
    void (*c__init_ens_pdaf)(int* , int* , int* , double* , double* , 
                             double* , int* ), 
    int* in_screen, int* outflag) noexcept nogil;

cdef extern void c__pdaf_init_forecast(
    void (*c__next_observation_pdaf)(int* , int* , int* , double* ), 
    void (*c__distribute_state_pdaf)(int* , double* ), 
    void (*c__prepoststep_pdaf)(int* , int* , int* , int* , int* , 
                                double* , double* , double* , int* ), 
    int* outflag) noexcept nogil;

cdef extern void c__pdaf_local_weight(int* wtype, int* rtype, 
    double* cradius, double* sradius, double* distance, int* nrows, 
    int* ncols, double* a, double* var_obs, double* weight, 
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_local_weights(int* wtype, double* cradius, 
    double* sradius, int* dim, double* distance, double* weight, 
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_print_filter_types(
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_print_da_types(
    int* verbose) noexcept nogil;

cdef extern void c__pdaf_print_info(
    int* printtype) noexcept nogil;

cdef extern void c__pdaf_reset_forget(
    double* forget_in) noexcept nogil;

cdef extern void c__pdaf_sampleens(int* dim, int* dim_ens, double* modes, 
    double* svals, double* state, double* ens, int* verbose, 
    int* flag) noexcept nogil;

