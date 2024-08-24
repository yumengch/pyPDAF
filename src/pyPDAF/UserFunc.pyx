import numpy as np
import warnings


cdef void c__add_obs_err_pdaf (int* step, int* dim_obs_p, double* C_p) noexcept:

    cdef double[::1,:] C_p_np = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> C_p, order='F')

    C_p_np = (<object>add_obs_err_pdaf)(step[0], dim_obs_p[0], C_p_np.base)


    cdef double[::1,:] C_p_new
    if C_p != &C_p_np[0,0]:
        C_p_new = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> C_p, order='F')
        C_p_new[...] = C_p_np
        warnings.warn("The memory address of C_p is changed in c__add_obs_err_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_ens_pdaf (int* filtertype, int* dim_p, int* dim_ens, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    assert dim_ens[0] > 1, "ensemble size must be > 1 for ensemble filters."

    cdef double[::1,:] uinv_np = np.asarray(<double[:dim_ens[0]-1:1,:dim_ens[0]-1]> uinv, order='F')
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')

    state_p_np, uinv_np, ens_p_np, flag[0] = (<object>init_ens_pdaf)(filtertype[0], dim_p[0], dim_ens[0], state_p_np.base, uinv_np.base, ens_p_np.base, flag[0])


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:(dim_ens[0]-1):1,:(dim_ens[0]-1)]> uinv, order='F')
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_ens_pdaf_single_member (int* filtertype, int* dim_p, int* dim_ens, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1,:] uinv_np = np.asarray(<double[:dim_ens[0]:1,:dim_ens[0]]> uinv, order='F')
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')

    state_p_np, uinv_np, ens_p_np, flag[0] = (<object>init_ens_pdaf_single_member)(filtertype[0], dim_p[0], dim_ens[0], state_p_np.base, uinv_np.base, ens_p_np.base, flag[0])


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:dim_ens[0]:1,:dim_ens[0]]> uinv, order='F')
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__next_observation_pdaf (int* stepnow, int* nsteps, int* doexit, double* time) noexcept:

    nsteps[0], doexit[0], time[0] = (<object>next_observation_pdaf)(stepnow[0], nsteps[0], doexit[0], time[0])



cdef void c__collect_state_pdaf (int* dim_p, double* state_p) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>collect_state_pdaf)(dim_p[0], state_p_np.base)


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__collect_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__distribute_state_pdaf (int* dim_p, double* state_p) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>distribute_state_pdaf)(dim_p[0], state_p_np.base)


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__distribute_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prepoststep_pdaf (int* step, int* dim_p, int* dim_ens, int* dim_ens_p, int* dim_obs_p, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    assert dim_ens[0] > 1, "ensemble size must be > 1 for ensemble filters."

    cdef double[::1,:] uinv_np = np.asarray(<double[:dim_ens[0]-1:1,:dim_ens[0]-1]> uinv, order='F')
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')

    state_p_np, uinv_np, ens_p_np = (<object>prepoststep_pdaf)(step[0], dim_p[0], dim_ens[0], dim_ens_p[0], dim_obs_p[0], state_p_np.base, uinv_np.base, ens_p_np.base, flag[0])


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:(dim_ens[0]-1):1,:(dim_ens[0]-1)]> uinv, order='F')
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_dim_obs_pdaf (int* step, int* dim_obs_p) noexcept:

    dim_obs_p[0] = (<object>init_dim_obs_pdaf)(step[0], dim_obs_p[0])



cdef void c__init_obs_pdaf (int* step, int* dim_obs_p, double* observation_p) noexcept:

    cdef double[::1] observation_p_np = np.asarray(<double[:dim_obs_p[0]]> observation_p)

    observation_p_np = (<object>init_obs_pdaf)(step[0], dim_obs_p[0], observation_p_np.base)


    cdef double[::1] observation_p_new
    if observation_p != &observation_p_np[0]:
        observation_p_new = np.asarray(<double[:dim_obs_p[0]]> observation_p)
        observation_p_new[...] = observation_p_np
        warnings.warn("The memory address of observation_p is changed in c__init_obs_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obs_covar_pdaf (int* step, int* dim_obs, int* dim_obs_p, double* covar, double* obs_p, bint* isdiag) noexcept:

    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)

    covar[0], isdiag[0] = (<object>init_obs_covar_pdaf)(step[0], dim_obs[0], dim_obs_p[0], covar[0], obs_p_np.base, isdiag[0])



cdef void c__init_obsvar_pdaf (int* step, int* dim_obs_p, double* obs_p, double* meanvar) noexcept:

    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)

    meanvar[0] = (<object>init_obsvar_pdaf)(step[0], dim_obs_p[0], obs_p_np.base, meanvar[0])



cdef void c__prodRinvA_pdaf (int* step, int* dim_obs_p, int* rank, double* obs_p, double* A_p, double* C_p) noexcept:

    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)
    cdef double[::1,:] A_p_np = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> A_p, order='F')
    cdef double[::1,:] C_p_np = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> C_p, order='F')

    C_p_np = (<object>prodRinvA_pdaf)(step[0], dim_obs_p[0], rank[0], obs_p_np.base, A_p_np.base, C_p_np.base)


    cdef double[::1,:] C_p_new
    if C_p != &C_p_np[0,0]:
        C_p_new = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> C_p, order='F')
        C_p_new[...] = C_p_np
        warnings.warn("The memory address of C_p is changed in c__prodRinvA_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    m_state_p_np = (<object>obs_op_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np.base, m_state_p_np.base)


    cdef double[::1] m_state_p_new
    if m_state_p != &m_state_p_np[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]]> m_state_p)
        m_state_p_new[...] = m_state_p_np
        warnings.warn("The memory address of m_state_p is changed in c__obs_op_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__g2l_obs_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l, int* mstate_f, int* dim_p, int* mstate_l, int* dim_l) noexcept:

    cdef int[::1] mstate_f_np = np.asarray(<int[:dim_p[0]]> mstate_f)
    cdef int[::1] mstate_l_np = np.asarray(<int[:dim_l[0]]> mstate_l)

    mstate_l_np = (<object>g2l_obs_pdaf)(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0], mstate_f_np.base, dim_p[0], mstate_l_np.base, dim_l[0])


    cdef int[::1] mstate_l_new
    if mstate_l != &mstate_l_np[0]:
        mstate_l_new = np.asarray(<int[:dim_l[0]]> mstate_l)
        mstate_l_new[...] = mstate_l_np
        warnings.warn("The memory address of mstate_l is changed in c__g2l_obs_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__g2l_state_pdaf (int* step, int* domain_p, int* dim_p, double* state_p, int* dim_l, double* state_l) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1] state_l_np = np.asarray(<double[:dim_l[0]]> state_l)

    state_l_np = (<object>g2l_state_pdaf)(step[0], domain_p[0], dim_p[0], state_p_np.base, dim_l[0], state_l_np.base)


    cdef double[::1] state_l_new
    if state_l != &state_l_np[0]:
        state_l_new = np.asarray(<double[:dim_l[0]]> state_l)
        state_l_new[...] = state_l_np
        warnings.warn("The memory address of state_l is changed in c__g2l_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_dim_l_pdaf (int* step, int* domain_p, int* dim_l) noexcept:

    dim_l[0] = (<object>init_dim_l_pdaf)(step[0], domain_p[0], dim_l[0])



cdef void c__init_dim_obs_f_pdaf (int* step, int* dim_obs_f) noexcept:

    dim_obs_f[0] = (<object>init_dim_obs_f_pdaf)(step[0], dim_obs_f[0])



cdef void c__init_dim_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l) noexcept:

    dim_obs_l[0] = (<object>init_dim_obs_l_pdaf)(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0])



cdef void c__init_n_domains_p_pdaf (int* step, int* n_domains_p) noexcept:

    n_domains_p[0] = (<object>init_n_domains_p_pdaf)(step[0], n_domains_p[0])



cdef void c__init_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f) noexcept:

    cdef double[::1] observation_f_np = np.asarray(<double[:dim_obs_f[0]]> observation_f)

    observation_f_np = (<object>init_obs_f_pdaf)(step[0], dim_obs_f[0], observation_f_np.base)


    cdef double[::1] observation_f_new
    if observation_f != &observation_f_np[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]]> observation_f)
        observation_f_new[...] = observation_f_np
        warnings.warn("The memory address of observation_f is changed in c__init_obs_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* observation_l) noexcept:

    cdef double[::1] observation_l_np = np.asarray(<double[:dim_obs_l[0]]> observation_l)

    observation_l_np = (<object>init_obs_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], observation_l_np.base)


    cdef double[::1] observation_l_new
    if observation_l != &observation_l_np[0]:
        observation_l_new = np.asarray(<double[:dim_obs_l[0]]> observation_l)
        observation_l_new[...] = observation_l_np
        warnings.warn("The memory address of observation_l is changed in c__init_obs_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obsvar_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, int* dim_obs_p, double* meanvar_l) noexcept:

    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_p[0]]> obs_l)

    meanvar_l[0] = (<object>init_obsvar_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np.base, dim_obs_p[0], meanvar_l[0])



cdef void c__init_obserr_f_pdaf (int* step, int* dim_obs_f, double* obs_f, double* obserr_f) noexcept:

    cdef double[::1] obs_f_np = np.asarray(<double[:dim_obs_f[0]]> obs_f)
    cdef double[::1] obserr_f_np = np.asarray(<double[:dim_obs_f[0]]> obserr_f)

    obserr_f_np = (<object>init_obserr_f_pdaf)(step[0], dim_obs_f[0], obs_f_np.base, obserr_f_np.base)


    cdef double[::1] obserr_f_new
    if obserr_f != &obserr_f_np[0]:
        obserr_f_new = np.asarray(<double[:dim_obs_f[0]]> obserr_f)
        obserr_f_new[...] = obserr_f_np
        warnings.warn("The memory address of obserr_f is changed in c__init_obserr_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__l2g_state_pdaf (int* step, int* domain_p, int* dim_l, double* state_l, int* dim_p, double* state_p) noexcept:

    cdef double[::1] state_l_np = np.asarray(<double[:dim_l[0]]> state_l)
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>l2g_state_pdaf)(step[0], domain_p[0], dim_l[0], state_l_np.base, dim_p[0], state_p_np.base)


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__l2g_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_f_pdaf (int* step, int* dim_p, int* dim_obs_f, double* state_p, double* m_state_f) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1] m_state_f_np = np.asarray(<double[:dim_obs_f[0]]> m_state_f)

    m_state_f_np = (<object>obs_op_f_pdaf)(step[0], dim_p[0], dim_obs_f[0], state_p_np.base, m_state_f_np.base)


    cdef double[::1] m_state_f_new
    if m_state_f != &m_state_f_np[0]:
        m_state_f_new = np.asarray(<double[:dim_obs_f[0]]> m_state_f)
        m_state_f_new[...] = m_state_f_np
        warnings.warn("The memory address of m_state_f is changed in c__obs_op_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prodRinvA_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* A_l, double* C_l) noexcept:

    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    cdef double[::1,:] A_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> A_l, order='F')
    cdef double[::1,:] C_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> C_l, order='F')

    C_l_np = (<object>prodRinvA_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np.base, A_l_np.base, C_l_np.base)


    cdef double[::1,:] C_l_new
    if C_l != &C_l_np[0,0]:
        C_l_new = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> C_l, order='F')
        C_l_new[...] = C_l_np
        warnings.warn("The memory address of C_l is changed in c__prodRinvA_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__localize_covar_pdaf (int* dim_p, int* dim_obs, double* hp_p, double* hph) noexcept:

    cdef double[::1,:] hp_p_np = np.asarray(<double[:dim_obs[0]:1,:dim_p[0]]> hp_p, order='F')
    cdef double[::1,:] hph_np = np.asarray(<double[:dim_obs[0]:1,:dim_obs[0]]> hph, order='F')

    hp_p_np, hph_np = (<object>localize_covar_pdaf)(dim_p[0], dim_obs[0], hp_p_np.base, hph_np.base)


    cdef double[::1,:] hp_p_new
    if hp_p != &hp_p_np[0,0]:
        hp_p_new = np.asarray(<double[:dim_obs[0]:1,:dim_p[0]]> hp_p, order='F')
        hp_p_new[...] = hp_p_np
        warnings.warn("The memory address of hp_p is changed in c__localize_covar_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] hph_new
    if hph != &hph_np[0,0]:
        hph_new = np.asarray(<double[:dim_obs[0]:1,:dim_obs[0]]> hph, order='F')
        hph_new[...] = hph_np
        warnings.warn("The memory address of hph is changed in c__localize_covar_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__likelihood_pdaf (int* step, int* dim_obs_p, double* obs_p, double* resid, double* likely) noexcept:

    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)
    cdef double[::1] resid_np = np.asarray(<double[:dim_obs_p[0]]> resid)

    likely[0] = (<object>likelihood_pdaf)(step[0], dim_obs_p[0], obs_p_np.base, resid_np.base, likely[0])



cdef void c__likelihood_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* likely_l) noexcept:

    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    cdef double[::1] resid_l_np = np.asarray(<double[:dim_obs_l[0]]> resid_l)

    likely_l[0] = (<object>likelihood_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np.base, resid_l_np.base, likely_l[0])



cdef void c__get_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f) noexcept:

    cdef double[::1] observation_f_np = np.asarray(<double[:dim_obs_f[0]]> observation_f)

    observation_f_np = (<object>get_obs_f_pdaf)(step[0], dim_obs_f[0], observation_f_np.base)


    cdef double[::1] observation_f_new
    if observation_f != &observation_f_np[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]]> observation_f)
        observation_f_new[...] = observation_f_np
        warnings.warn("The memory address of observation_f is changed in c__get_obs_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_adj_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cv_ens_p, double* ens_p, double* Vcv_p, double* cv_p) noexcept:

    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
    cdef double[::1] Vcv_p_np = np.asarray(<double[:dim_p[0]]> Vcv_p)
    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cv_ens_p[0]]> cv_p)

    cv_p_np = (<object>cvt_adj_ens_pdaf)(iter[0], dim_p[0], dim_ens[0], dim_cv_ens_p[0], ens_p_np.base, Vcv_p_np.base, cv_p_np.base)


    cdef double[::1] cv_p_new
    if cv_p != &cv_p_np[0]:
        cv_p_new = np.asarray(<double[:dim_cv_ens_p[0]]> cv_p)
        cv_p_new[...] = cv_p_np
        warnings.warn("The memory address of cv_p is changed in c__cvt_adj_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_adj_pdaf (int* iter, int* dim_p, int* dim_cvec, double* Vcv_p, double* cv_p) noexcept:

    cdef double[::1] Vcv_p_np = np.asarray(<double[:dim_p[0]]> Vcv_p)
    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cvec[0]]> cv_p)

    cv_p_np = (<object>cvt_adj_pdaf)(iter[0], dim_p[0], dim_cvec[0], Vcv_p_np.base, cv_p_np.base)


    cdef double[::1] cv_p_new
    if cv_p != &cv_p_np[0]:
        cv_p_new = np.asarray(<double[:dim_cvec[0]]> cv_p)
        cv_p_new[...] = cv_p_np
        warnings.warn("The memory address of cv_p is changed in c__cvt_adj_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_pdaf (int* iter, int* dim_p, int* dim_cvec, double* cv_p, double* Vv_p) noexcept:

    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cvec[0]]> cv_p)
    cdef double[::1] Vv_p_np = np.asarray(<double[:dim_p[0]]> Vv_p)

    Vv_p_np = (<object>cvt_pdaf)(iter[0], dim_p[0], dim_cvec[0], cv_p_np.base, Vv_p_np.base)


    cdef double[::1] Vv_p_new
    if Vv_p != &Vv_p_np[0]:
        Vv_p_new = np.asarray(<double[:dim_p[0]]> Vv_p)
        Vv_p_new[...] = Vv_p_np
        warnings.warn("The memory address of Vv_p is changed in c__cvt_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cvec_ens, double* ens_p, double* v_p, double* Vv_p) noexcept:

    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
    cdef double[::1] v_p_np = np.asarray(<double[:dim_cvec_ens[0]]> v_p)
    cdef double[::1] Vv_p_np = np.asarray(<double[:dim_p[0]]> Vv_p)

    Vv_p_np = (<object>cvt_ens_pdaf)(iter[0], dim_p[0], dim_ens[0], dim_cvec_ens[0], ens_p_np.base, v_p_np.base, Vv_p_np.base)


    cdef double[::1] Vv_p_new
    if Vv_p != &Vv_p_np[0]:
        Vv_p_new = np.asarray(<double[:dim_p[0]]> Vv_p)
        Vv_p_new[...] = Vv_p_np
        warnings.warn("The memory address of Vv_p is changed in c__cvt_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_adj_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    state_p_np = (<object>obs_op_adj_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np.base, m_state_p_np.base)


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__obs_op_adj_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_lin_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    m_state_p_np = (<object>obs_op_lin_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np.base, m_state_p_np.base)


    cdef double[::1] m_state_p_new
    if m_state_p != &m_state_p_np[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]]> m_state_p)
        m_state_p_new[...] = m_state_p_np
        warnings.warn("The memory address of m_state_p is changed in c__obs_op_lin_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__dist_stateinc_pdaf (int* dim_p, double* state_inc_p, int* first, int* steps) noexcept:

    cdef double[::1] state_inc_p_np = np.asarray(<double[:dim_p[0]]> state_inc_p)

    (<object>dist_stateinc_pdaf)(dim_p[0], state_inc_p_np.base, first[0], steps[0])



cdef void c__prodRinvA_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* gamma, double* likely_l) noexcept:

    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    cdef double[::1] resid_l_np = np.asarray(<double[:dim_obs_l[0]]> resid_l)

    likely_l[0] = (<object>prodRinvA_hyb_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np.base, resid_l_np.base, gamma[0], likely_l[0])



cdef void c__likelihood_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* gamma, double* A_l, double* C_l) noexcept:

    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    cdef double[::1,:] A_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> A_l, order='F')
    cdef double[::1,:] C_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> C_l, order='F')

    C_l_np = (<object>likelihood_hyb_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np.base, gamma[0], A_l_np.base, C_l_np.base)


    cdef double[::1,:] C_l_new
    if C_l != &C_l_np[0,0]:
        C_l_new = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> C_l, order='F')
        C_l_new[...] = C_l_np
        warnings.warn("The memory address of C_l is changed in c__likelihood_hyb_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


