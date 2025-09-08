from pyPDAF.cfi_binding cimport CFI_cdesc_t
cdef extern void c__pdaf_diag_ensmean(int* dim, int* dim_ens,
    double* state, double* ens, int* status) noexcept nogil;

cdef extern void c__pdaf_diag_stddev_nompi(int* dim, int* dim_ens,
    double* state, double* ens, double* stddev, int* do_mean,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_stddev(int* dim_p, int* dim_ens,
    double* state_p, double* ens_p, double* stddev_g, int* do_mean,
    int* comm_filter, int* status) noexcept nogil;

cdef extern void c__pdaf_diag_variance_nompi(int* dim, int* dim_ens,
    double* state, double* ens, double* variance, double* stddev,
    int* do_mean, int* do_stddev, int* status) noexcept nogil;

cdef extern void c__pdaf_diag_variance(int* dim_p, int* dim_ens,
    double* state_p, double* ens_p, double* variance_p, double* stddev_g,
    int* do_mean, int* do_stddev, int* comm_filter,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_rmsd_nompi(int* dim_p, double* statea_p,
    double* stateb_p, double* rmsd_p,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_rmsd(int* dim_p, double* statea_p,
    double* stateb_p, double* rmsd_g, int* comm_filter,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_crps_mpi(int* dim_p, int* dim_ens,
    int* element, double* oens, double* obs, int* comm_filter,
    int* mype_filter, int* npes_filter, double* crps, double* reli,
    double* pot_crps, double* uncert,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_crps_nompi(int* dim, int* dim_ens,
    int* element, double* oens, double* obs, double* crps, double* reli,
    double* resol, double* uncert, int* status) noexcept nogil;

cdef extern void c__pdaf_diag_effsample(int* dim_sample, double* weights,
    double* n_eff) noexcept nogil;

cdef extern void c__pdaf_diag_ensstats(int* dim, int* dim_ens,
    int* element, double* state, double* ens, double* skewness,
    double* kurtosis, int* status) noexcept nogil;

cdef extern void c__pdaf_diag_compute_moments(int* dim_p, int* dim_ens,
    double* ens, int* kmax, double* moments,
    int* bias) noexcept nogil;

cdef extern void c__pdaf_diag_histogram(int* ncall, int* dim, int* dim_ens,
    int* element, double* state, double* ens, int* hist, double* delta,
    int* status) noexcept nogil;

cdef extern void c__pdaf_diag_reliability_budget(int* n_times,
    int* dim_ens, int* dim_p, double* ens_p, double* obsvar, double* obs_p,
    double* budget, double* bias_2) noexcept nogil;

