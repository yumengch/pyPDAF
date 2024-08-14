import numpy as np
import warnings


cdef void c__add_obs_err_pdaf (int* step, int* dim_obs_p, double* c_p) noexcept:


    cdef int c_p_dim = dim_obs_p[0]*dim_obs_p[0]
    c_p_np = np.asarray(<double[:c_p_dim]> c_p).reshape(dim_obs_p[0],dim_obs_p[0], order='F')

    c_p_np = (<object>add_obs_err_pdaf)(step[0], dim_obs_p[0], c_p_np)


    cdef double[::1] c_p_view = c_p_np.ravel(order='F')
    if c_p != &c_p_view[0]:
        c_p_new = np.asarray(<double[:c_p_dim]> c_p).reshape(dim_obs_p[0],dim_obs_p[0], order='F')
        c_p_new[:] = c_p_np[:]
        warnings.RuntimeWarning("The memory address of c_p is changed in c__add_obs_err_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_ens_pdaf (int* filtertype, int* dim_p, int* dim_ens, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    assert dim_ens[0] > 1, "ensemble size must be > 1 for ensemble filters."


    cdef int uinv_dim = (dim_ens[0]-1)*(dim_ens[0]-1)
    uinv_np = np.asarray(<double[:uinv_dim]> uinv).reshape(dim_ens[0]-1,dim_ens[0]-1, order='F')

    cdef int ens_p_dim = dim_p[0]*dim_ens[0]
    ens_p_np = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')

    state_p_np, uinv_np, ens_p_np, flag[0] = (<object>init_ens_pdaf)(filtertype[0], dim_p[0], dim_ens[0], state_p_np, uinv_np, ens_p_np, flag[0])


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] uinv_view = uinv_np.ravel(order='F')
    if uinv != &uinv_view[0]:
        uinv_new = np.asarray(<double[:uinv_dim]> uinv).reshape((dim_ens[0]-1),(dim_ens[0]-1), order='F')
        uinv_new[:] = uinv_np[:]
        warnings.RuntimeWarning("The memory address of uinv is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] ens_p_view = ens_p_np.ravel(order='F')
    if ens_p != &ens_p_view[0]:
        ens_p_new = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')
        ens_p_new[:] = ens_p_np[:]
        warnings.RuntimeWarning("The memory address of ens_p is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_ens_pdaf_single_member (int* filtertype, int* dim_p, int* dim_ens, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    cdef int uinv_dim = dim_ens[0]*dim_ens[0]
    uinv_np = np.asarray(<double[:uinv_dim]> uinv).reshape(dim_ens[0],dim_ens[0], order='F')

    cdef int ens_p_dim = dim_p[0]*dim_ens[0]
    ens_p_np = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')

    state_p_np, uinv_np, ens_p_np, flag[0] = (<object>init_ens_pdaf_single_member)(filtertype[0], dim_p[0], dim_ens[0], state_p_np, uinv_np, ens_p_np, flag[0])


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] uinv_view = uinv_np.ravel(order='F')
    if uinv != &uinv_view[0]:
        uinv_new = np.asarray(<double[:uinv_dim]> uinv).reshape(dim_ens[0],dim_ens[0], order='F')
        uinv_new[:] = uinv_np[:]
        warnings.RuntimeWarning("The memory address of uinv is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] ens_p_view = ens_p_np.ravel(order='F')
    if ens_p != &ens_p_view[0]:
        ens_p_new = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')
        ens_p_new[:] = ens_p_np[:]
        warnings.RuntimeWarning("The memory address of ens_p is changed in c__init_ens_pdaf_single_member." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__next_observation_pdaf (int* stepnow, int* nsteps, int* doexit, double* time) noexcept:

    nsteps[0], doexit[0], time[0] = (<object>next_observation_pdaf)(stepnow[0], nsteps[0], doexit[0], time[0])



cdef void c__collect_state_pdaf (int* dim_p, double* state_p) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>collect_state_pdaf)(dim_p[0], state_p_np)


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__collect_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__distribute_state_pdaf (int* dim_p, double* state_p) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>distribute_state_pdaf)(dim_p[0], state_p_np)


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__distribute_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__prepoststep_pdaf (int* step, int* dim_p, int* dim_ens, int* dim_ens_p, int* dim_obs_p, double* state_p, double* uinv, double* ens_p, int* flag) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    assert dim_ens[0] > 1, "ensemble size must be > 1 for ensemble filters."


    cdef int uinv_dim = (dim_ens[0]-1)*(dim_ens[0]-1)
    uinv_np = np.asarray(<double[:uinv_dim]> uinv).reshape(dim_ens[0]-1,dim_ens[0]-1, order='F')

    cdef int ens_p_dim = dim_p[0]*dim_ens[0]
    ens_p_np = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')

    state_p_np, uinv_np, ens_p_np = (<object>prepoststep_pdaf)(step[0], dim_p[0], dim_ens[0], dim_ens_p[0], dim_obs_p[0], state_p_np, uinv_np, ens_p_np, flag[0])


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] uinv_view = uinv_np.ravel(order='F')
    if uinv != &uinv_view[0]:
        uinv_new = np.asarray(<double[:uinv_dim]> uinv).reshape((dim_ens[0]-1),(dim_ens[0]-1), order='F')
        uinv_new[:] = uinv_np[:]
        warnings.RuntimeWarning("The memory address of uinv is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] ens_p_view = ens_p_np.ravel(order='F')
    if ens_p != &ens_p_view[0]:
        ens_p_new = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')
        ens_p_new[:] = ens_p_np[:]
        warnings.RuntimeWarning("The memory address of ens_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_dim_obs_pdaf (int* step, int* dim_obs_p) noexcept:

    dim_obs_p[0] = (<object>init_dim_obs_pdaf)(step[0], dim_obs_p[0])



cdef void c__init_obs_pdaf (int* step, int* dim_obs_p, double* observation_p) noexcept:

    observation_p_np = np.asarray(<double[:dim_obs_p[0]]> observation_p)

    observation_p_np = (<object>init_obs_pdaf)(step[0], dim_obs_p[0], observation_p_np)


    cdef double[::1] observation_p_view = observation_p_np
    if observation_p != &observation_p_view[0]:
        observation_p_new = np.asarray(<double[:dim_obs_p[0]]> observation_p)
        observation_p_new[:] = observation_p_np[:]
        warnings.RuntimeWarning("The memory address of observation_p is changed in c__init_obs_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_obs_covar_pdaf (int* step, int* dim_obs, int* dim_obs_p, double* covar, double* obs_p, bint* isdiag) noexcept:

    obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)

    covar[0], isdiag[0] = (<object>init_obs_covar_pdaf)(step[0], dim_obs[0], dim_obs_p[0], covar[0], obs_p_np, isdiag[0])



cdef void c__init_obsvar_pdaf (int* step, int* dim_obs_p, double* obs_p, double* meanvar) noexcept:

    obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)

    meanvar[0] = (<object>init_obsvar_pdaf)(step[0], dim_obs_p[0], obs_p_np, meanvar[0])



cdef void c__prodrinva_pdaf (int* step, int* dim_obs_p, int* rank, double* obs_p, double* a_p, double* c_p) noexcept:

    obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)

    cdef int a_p_dim = dim_obs_p[0]*rank[0]
    a_p_np = np.asarray(<double[:a_p_dim]> a_p).reshape(dim_obs_p[0],rank[0], order='F')
    c_p_np = np.asarray(<double[:a_p_dim]> c_p).reshape(dim_obs_p[0],rank[0], order='F')

    c_p_np = (<object>prodrinva_pdaf)(step[0], dim_obs_p[0], rank[0], obs_p_np, a_p_np, c_p_np)


    cdef double[::1] c_p_view = c_p_np.ravel(order='F')
    if c_p != &c_p_view[0]:
        c_p_new = np.asarray(<double[:a_p_dim]> c_p).reshape(dim_obs_p[0],rank[0], order='F')
        c_p_new[:] = c_p_np[:]
        warnings.RuntimeWarning("The memory address of c_p is changed in c__prodrinva_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__obs_op_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    m_state_p_np = (<object>obs_op_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)


    cdef double[::1] m_state_p_view = m_state_p_np
    if m_state_p != &m_state_p_view[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]]> m_state_p)
        m_state_p_new[:] = m_state_p_np[:]
        warnings.RuntimeWarning("The memory address of m_state_p is changed in c__obs_op_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__g2l_obs_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l, int* mstate_f, int* dim_p, int* mstate_l, int* dim_l) noexcept:

    mstate_f_np = np.asarray(<int[:dim_p[0]]> mstate_f)
    mstate_l_np = np.asarray(<int[:dim_l[0]]> mstate_l)

    mstate_l_np = (<object>g2l_obs_pdaf)(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0], mstate_f_np, dim_p[0], mstate_l_np, dim_l[0])


    cdef int[::1] mstate_l_view = mstate_l_np
    if mstate_l != &mstate_l_view[0]:
        mstate_l_new = np.asarray(<int[:dim_l[0]]> mstate_l)
        mstate_l_new[:] = mstate_l_np[:]
        warnings.RuntimeWarning("The memory address of mstate_l is changed in c__g2l_obs_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__g2l_state_pdaf (int* step, int* domain_p, int* dim_p, double* state_p, int* dim_l, double* state_l) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    state_l_np = np.asarray(<double[:dim_l[0]]> state_l)

    state_l_np = (<object>g2l_state_pdaf)(step[0], domain_p[0], dim_p[0], state_p_np, dim_l[0], state_l_np)


    cdef double[::1] state_l_view = state_l_np
    if state_l != &state_l_view[0]:
        state_l_new = np.asarray(<double[:dim_l[0]]> state_l)
        state_l_new[:] = state_l_np[:]
        warnings.RuntimeWarning("The memory address of state_l is changed in c__g2l_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_dim_l_pdaf (int* step, int* domain_p, int* dim_l) noexcept:

    dim_l[0] = (<object>init_dim_l_pdaf)(step[0], domain_p[0], dim_l[0])



cdef void c__init_dim_obs_f_pdaf (int* step, int* dim_obs_f) noexcept:

    dim_obs_f[0] = (<object>init_dim_obs_f_pdaf)(step[0], dim_obs_f[0])



cdef void c__init_dim_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_f, int* dim_obs_l) noexcept:

    dim_obs_l[0] = (<object>init_dim_obs_l_pdaf)(domain_p[0], step[0], dim_obs_f[0], dim_obs_l[0])



cdef void c__init_n_domains_p_pdaf (int* step, int* n_domains_p) noexcept:

    n_domains_p[0] = (<object>init_n_domains_p_pdaf)(step[0], n_domains_p[0])



cdef void c__init_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f) noexcept:

    observation_f_np = np.asarray(<double[:dim_obs_f[0]]> observation_f)

    observation_f_np = (<object>init_obs_f_pdaf)(step[0], dim_obs_f[0], observation_f_np)


    cdef double[::1] observation_f_view = observation_f_np
    if observation_f != &observation_f_view[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]]> observation_f)
        observation_f_new[:] = observation_f_np[:]
        warnings.RuntimeWarning("The memory address of observation_f is changed in c__init_obs_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_obs_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* observation_l) noexcept:

    observation_l_np = np.asarray(<double[:dim_obs_l[0]]> observation_l)

    observation_l_np = (<object>init_obs_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], observation_l_np)


    cdef double[::1] observation_l_view = observation_l_np
    if observation_l != &observation_l_view[0]:
        observation_l_new = np.asarray(<double[:dim_obs_l[0]]> observation_l)
        observation_l_new[:] = observation_l_np[:]
        warnings.RuntimeWarning("The memory address of observation_l is changed in c__init_obs_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__init_obsvar_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, int* dim_obs_p, double* meanvar_l) noexcept:

    obs_l_np = np.asarray(<double[:dim_obs_p[0]]> obs_l)

    meanvar_l[0] = (<object>init_obsvar_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np, dim_obs_p[0], meanvar_l[0])



cdef void c__init_obserr_f_pdaf (int* step, int* dim_obs_f, double* obs_f, double* obserr_f) noexcept:

    obs_f_np = np.asarray(<double[:dim_obs_f[0]]> obs_f)
    obserr_f_np = np.asarray(<double[:dim_obs_f[0]]> obserr_f)

    obserr_f_np = (<object>init_obserr_f_pdaf)(step[0], dim_obs_f[0], obs_f_np, obserr_f_np)


    cdef double[::1] obserr_f_view = obserr_f_np
    if obserr_f != &obserr_f_view[0]:
        obserr_f_new = np.asarray(<double[:dim_obs_f[0]]> obserr_f)
        obserr_f_new[:] = obserr_f_np[:]
        warnings.RuntimeWarning("The memory address of obserr_f is changed in c__init_obserr_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__l2g_state_pdaf (int* step, int* domain_p, int* dim_l, double* state_l, int* dim_p, double* state_p) noexcept:

    state_l_np = np.asarray(<double[:dim_l[0]]> state_l)
    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    state_p_np = (<object>l2g_state_pdaf)(step[0], domain_p[0], dim_l[0], state_l_np, dim_p[0], state_p_np)


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__l2g_state_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__obs_op_f_pdaf (int* step, int* dim_p, int* dim_obs_f, double* state_p, double* m_state_f) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    m_state_f_np = np.asarray(<double[:dim_obs_f[0]]> m_state_f)

    m_state_f_np = (<object>obs_op_f_pdaf)(step[0], dim_p[0], dim_obs_f[0], state_p_np, m_state_f_np)


    cdef double[::1] m_state_f_view = m_state_f_np
    if m_state_f != &m_state_f_view[0]:
        m_state_f_new = np.asarray(<double[:dim_obs_f[0]]> m_state_f)
        m_state_f_new[:] = m_state_f_np[:]
        warnings.RuntimeWarning("The memory address of m_state_f is changed in c__obs_op_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__prodrinva_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* a_l, double* c_l) noexcept:

    obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)

    cdef int a_l_dim = dim_obs_l[0]*rank[0]
    a_l_np = np.asarray(<double[:a_l_dim]> a_l).reshape(dim_obs_l[0],rank[0], order='F')
    c_l_np = np.asarray(<double[:a_l_dim]> c_l).reshape(dim_obs_l[0],rank[0], order='F')

    c_l_np = (<object>prodrinva_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np, a_l_np, c_l_np)


    cdef double[::1] c_l_view = c_l_np.ravel(order='F')
    if c_l != &c_l_view[0]:
        c_l_new = np.asarray(<double[:a_l_dim]> c_l).reshape(dim_obs_l[0],rank[0], order='F')
        c_l_new[:] = c_l_np[:]
        warnings.RuntimeWarning("The memory address of c_l is changed in c__prodrinva_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__localize_covar_pdaf (int* dim_p, int* dim_obs, double* hp_p, double* hph) noexcept:


    cdef int hp_p_dim = dim_obs[0]*dim_p[0]
    hp_p_np = np.asarray(<double[:hp_p_dim]> hp_p).reshape(dim_obs[0],dim_p[0], order='F')

    cdef int hph_dim = dim_obs[0]*dim_obs[0]
    hph_np = np.asarray(<double[:hph_dim]> hph).reshape(dim_obs[0],dim_obs[0], order='F')

    hp_p_np, hph_np = (<object>localize_covar_pdaf)(dim_p[0], dim_obs[0], hp_p_np, hph_np)


    cdef double[::1] hp_p_view = hp_p_np.ravel(order='F')
    if hp_p != &hp_p_view[0]:
        hp_p_new = np.asarray(<double[:hp_p_dim]> hp_p).reshape(dim_obs[0],dim_p[0], order='F')
        hp_p_new[:] = hp_p_np[:]
        warnings.RuntimeWarning("The memory address of hp_p is changed in c__localize_covar_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")

    cdef double[::1] hph_view = hph_np.ravel(order='F')
    if hph != &hph_view[0]:
        hph_new = np.asarray(<double[:hph_dim]> hph).reshape(dim_obs[0],dim_obs[0], order='F')
        hph_new[:] = hph_np[:]
        warnings.RuntimeWarning("The memory address of hph is changed in c__localize_covar_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__likelihood_pdaf (int* step, int* dim_obs_p, double* obs_p, double* resid, double* likely) noexcept:

    obs_p_np = np.asarray(<double[:dim_obs_p[0]]> obs_p)
    resid_np = np.asarray(<double[:dim_obs_p[0]]> resid)

    likely[0] = (<object>likelihood_pdaf)(step[0], dim_obs_p[0], obs_p_np, resid_np, likely[0])



cdef void c__likelihood_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* likely_l) noexcept:

    obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    resid_l_np = np.asarray(<double[:dim_obs_l[0]]> resid_l)

    likely_l[0] = (<object>likelihood_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np, resid_l_np, likely_l[0])



cdef void c__get_obs_f_pdaf (int* step, int* dim_obs_f, double* observation_f) noexcept:

    observation_f_np = np.asarray(<double[:dim_obs_f[0]]> observation_f)

    observation_f_np = (<object>get_obs_f_pdaf)(step[0], dim_obs_f[0], observation_f_np)


    cdef double[::1] observation_f_view = observation_f_np
    if observation_f != &observation_f_view[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]]> observation_f)
        observation_f_new[:] = observation_f_np[:]
        warnings.RuntimeWarning("The memory address of observation_f is changed in c__get_obs_f_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__cvt_adj_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cv_ens_p, double* ens_p, double* vcv_p, double* cv_p) noexcept:


    cdef int ens_p_dim = dim_p[0]*dim_ens[0]
    ens_p_np = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')
    vcv_p_np = np.asarray(<double[:dim_p[0]]> vcv_p)
    cv_p_np = np.asarray(<double[:dim_cv_ens_p[0]]> cv_p)

    cv_p_np = (<object>cvt_adj_ens_pdaf)(iter[0], dim_p[0], dim_ens[0], dim_cv_ens_p[0], ens_p_np, vcv_p_np, cv_p_np)


    cdef double[::1] cv_p_view = cv_p_np
    if cv_p != &cv_p_view[0]:
        cv_p_new = np.asarray(<double[:dim_cv_ens_p[0]]> cv_p)
        cv_p_new[:] = cv_p_np[:]
        warnings.RuntimeWarning("The memory address of cv_p is changed in c__cvt_adj_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__cvt_adj_pdaf (int* iter, int* dim_p, int* dim_cvec, double* vcv_p, double* cv_p) noexcept:

    vcv_p_np = np.asarray(<double[:dim_p[0]]> vcv_p)
    cv_p_np = np.asarray(<double[:dim_cvec[0]]> cv_p)

    cv_p_np = (<object>cvt_adj_pdaf)(iter[0], dim_p[0], dim_cvec[0], vcv_p_np, cv_p_np)


    cdef double[::1] cv_p_view = cv_p_np
    if cv_p != &cv_p_view[0]:
        cv_p_new = np.asarray(<double[:dim_cvec[0]]> cv_p)
        cv_p_new[:] = cv_p_np[:]
        warnings.RuntimeWarning("The memory address of cv_p is changed in c__cvt_adj_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__cvt_pdaf (int* iter, int* dim_p, int* dim_cvec, double* cv_p, double* vv_p) noexcept:

    cv_p_np = np.asarray(<double[:dim_cvec[0]]> cv_p)
    vv_p_np = np.asarray(<double[:dim_p[0]]> vv_p)

    vv_p_np = (<object>cvt_pdaf)(iter[0], dim_p[0], dim_cvec[0], cv_p_np, vv_p_np)


    cdef double[::1] vv_p_view = vv_p_np
    if vv_p != &vv_p_view[0]:
        vv_p_new = np.asarray(<double[:dim_p[0]]> vv_p)
        vv_p_new[:] = vv_p_np[:]
        warnings.RuntimeWarning("The memory address of vv_p is changed in c__cvt_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__cvt_ens_pdaf (int* iter, int* dim_p, int* dim_ens, int* dim_cvec_ens, double* ens_p, double* v_p, double* vv_p) noexcept:


    cdef int ens_p_dim = dim_p[0]*dim_ens[0]
    ens_p_np = np.asarray(<double[:ens_p_dim]> ens_p).reshape(dim_p[0],dim_ens[0], order='F')
    v_p_np = np.asarray(<double[:dim_cvec_ens[0]]> v_p)
    vv_p_np = np.asarray(<double[:dim_p[0]]> vv_p)

    vv_p_np = (<object>cvt_ens_pdaf)(iter[0], dim_p[0], dim_ens[0], dim_cvec_ens[0], ens_p_np, v_p_np, vv_p_np)


    cdef double[::1] vv_p_view = vv_p_np
    if vv_p != &vv_p_view[0]:
        vv_p_new = np.asarray(<double[:dim_p[0]]> vv_p)
        vv_p_new[:] = vv_p_np[:]
        warnings.RuntimeWarning("The memory address of vv_p is changed in c__cvt_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__obs_op_adj_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    state_p_np = (<object>obs_op_adj_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)


    cdef double[::1] state_p_view = state_p_np
    if state_p != &state_p_view[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[:] = state_p_np[:]
        warnings.RuntimeWarning("The memory address of state_p is changed in c__obs_op_adj_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__obs_op_lin_pdaf (int* step, int* dim_p, int* dim_obs_p, double* state_p, double* m_state_p) noexcept:

    state_p_np = np.asarray(<double[:dim_p[0]]> state_p)
    m_state_p_np = np.asarray(<double[:dim_obs_p[0]]> m_state_p)

    m_state_p_np = (<object>obs_op_lin_pdaf)(step[0], dim_p[0], dim_obs_p[0], state_p_np, m_state_p_np)


    cdef double[::1] m_state_p_view = m_state_p_np
    if m_state_p != &m_state_p_view[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]]> m_state_p)
        m_state_p_new[:] = m_state_p_np[:]
        warnings.RuntimeWarning("The memory address of m_state_p is changed in c__obs_op_lin_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


cdef void c__dist_stateinc_pdaf (int* dim_p, double* state_inc_p, int* first, int* steps) noexcept:

    state_inc_p_np = np.asarray(<double[:dim_p[0]]> state_inc_p)

    (<object>dist_stateinc_pdaf)(dim_p[0], state_inc_p_np, first[0], steps[0])



cdef void c__prodrinva_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, double* obs_l, double* resid_l, double* gamma, double* likely_l) noexcept:

    obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)
    resid_l_np = np.asarray(<double[:dim_obs_l[0]]> resid_l)

    likely_l[0] = (<object>prodrinva_hyb_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], obs_l_np, resid_l_np, gamma[0], likely_l[0])



cdef void c__likelihood_hyb_l_pdaf (int* domain_p, int* step, int* dim_obs_l, int* rank, double* obs_l, double* gamma, double* a_l, double* c_l) noexcept:

    obs_l_np = np.asarray(<double[:dim_obs_l[0]]> obs_l)

    cdef int a_l_dim = dim_obs_l[0]*rank[0]
    a_l_np = np.asarray(<double[:a_l_dim]> a_l).reshape(dim_obs_l[0],rank[0], order='F')
    c_l_np = np.asarray(<double[:a_l_dim]> c_l).reshape(dim_obs_l[0],rank[0], order='F')

    c_l_np = (<object>likelihood_hyb_l_pdaf)(domain_p[0], step[0], dim_obs_l[0], rank[0], obs_l_np, gamma[0], a_l_np, c_l_np)


    cdef double[::1] c_l_view = c_l_np.ravel(order='F')
    if c_l != &c_l_view[0]:
        c_l_new = np.asarray(<double[:a_l_dim]> c_l).reshape(dim_obs_l[0],rank[0], order='F')
        c_l_new[:] = c_l_np[:]
        warnings.RuntimeWarning("The memory address of c_l is changed in c__likelihood_hyb_l_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.")


