cdef void c__prepoststep_pdaf (int* step, int* dim_p, int* dim_ens, int* dim_ens_l, int* dim_obs_p, double* state_p, double* uinv, double* ens_p, int* flag) noexcept with gil:

    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]]> state_p)

    cdef int uinv_size = max(dim_ens[0]-1, 1)
    # if filtertype in [1, 3, 6, 7, 200]:
    #     uinv_size = max(dim_ens[0]-1, 1)
    # elif filtertype in [0, 4, 5, 9, 10, 11, 12]:
    #     uinv_size = dim_ens[0]
    # else:
    #     uinv_size = 1

    cdef double[::1,:] uinv_np = np.asarray(<double[:uinv_size:1,:uinv_size]> uinv, order='F')
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')

    state_p_np, uinv_np, ens_p_np = (<object>prepoststep_pdaf)(step[0], dim_p[0], dim_ens[0], dim_ens_l[0], dim_obs_p[0], state_p_np.base, uinv_np.base, ens_p_np.base, flag[0])


    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double [: dim_p [0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:uinv_size:1,:uinv_size]> uinv, order='F')
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)

    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order='F')
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__prepoststep_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
