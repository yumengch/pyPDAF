import sys
import numpy as np
import warnings

cdef void c__add_obs_err_pdaf(int* step, int* dim_obs_p,
    double* c_p) noexcept with gil:
    """Add the observation error covariance matrix to the matrix C.

    The input matrix is the projection of the ensemble
    covariance matrix onto the observation space that is computed
    during the analysis step of the EnKF.

    The operation is for the global observation space.
    Thus, it is independent of whether the filter is executed with or
    without parallelization.

    Parameters
    ----------
    step : int
        Current time step
    dim_obs_p : int
        Size of observation vector
    C_p : ndarray[np.float64, ndim=2]
        Matrix to which the observation error covariance matrix is added
        shape: (dim_obs_p, dim_obs_p)

    Returns
    -------
    C_p : ndarray[np.float64, ndim=2]
        Matrix with added obs error covariance
        shape: (dim_obs_p, dim_obs_p)
    """
    cdef double[::1,:] c_p_np = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> c_p, order="F")

    c_p_np = (<object>add_obs_err_pdaf)(step[0], dim_obs_p[0], c_p_np.base)

    cdef double[::1,:] c_p_new
    if c_p != &c_p_np[0,0]:
        c_p_new = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> c_p, order="F")
        c_p_new[...] = c_p_np
        warnings.warn("The memory address of c_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_ens_pdaf(int* filtertype, int* dim_p, int* dim_ens,
    double* state_p, double* uinv, double* ens_p, int* flag) noexcept with gil:
    """Fill the ensemble array that is provided by PDAF with an initial ensemble of model states.

    This function is called by :func:`pyPDAF.PDAF.init`. The initialised
    ensemble array will be distributed to model by :func:`pyPDAF.PDAF.init_forecast`.

    Parameters
    ----------
    filtertype : int
            filter type given in PDAF_init
    dim_p : int
            PE-local state dimension given by PDAF_init
    dim_ens : int
            number of ensemble members
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag
    """
    cdef size_t uinv_len = max(dim_ens[0]-1, 1)
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1,:] uinv_np = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")

    state_p_np,uinv_np,ens_p_np,flag[0] = (<object>init_ens_pdaf)(
                                                                  filtertype[0],
                                                                  dim_p[0],
                                                                  dim_ens[0],
                                                                  state_p_np.base,
                                                                  uinv_np.base,
                                                                  ens_p_np.base,
                                                                  flag[0])

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__next_observation_pdaf(int* stepnow, int* nsteps, int* doexit,
    double* time) noexcept with gil:
    """Get the number of time steps to be computed in the forecast phase.

    At the beginning of a forecast phase, this is called once by
        * :func:`pyPDAF.PDAF.init_forecast`
        * :func:`pyPDAF.PDAF3.assimilate_X`
        * ...

    Parameters
    ----------
    stepnow : int
            the current time step given by PDAF

    Returns
    -------
    nsteps : int
            number of forecast time steps until next assimilation;
            this can also be interpreted as
            number of assimilation function calls
            to perform a new assimilation
    doexit : int
            whether to exit forecasting (1 for exit)
    time : double
            current model (physical) time
    """
    nsteps[0],doexit[0],time[0] = (<object>next_observation_pdaf)(
                                                                  stepnow[0],
                                                                  nsteps[0],
                                                                  doexit[0],
                                                                  time[0])



cdef void c__collect_state_pdaf(int* dim_p, double* state_p) noexcept with gil:
    """Collect state vector from model/any arrays to pdaf arrays

    Parameters
    ----------
    dim_p : int
        pe-local state dimension

    state_p : ndarray[tuple[dim_p, ...], np.float64]
        local state vector

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        local state vector
    """
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")

    state_p_np = (<object>collect_state_pdaf)(dim_p[0], state_p_np.base)

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__distribute_state_pdaf(int* dim_p,
    double* state_p) noexcept with gil:
    """Distribute a state vector from pdaf to the model/any arrays

    Parameters
    ----------
    dim_p : int
            PE-local state dimension
    state_p : ndarray[np.float64, ndim=1]
            PE-local state vector
            Array shape: (dim_p)

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            PE-local state vector
            Array shape: (dim_p)
    """
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")

    state_p_np = (<object>distribute_state_pdaf)(dim_p[0], state_p_np.base)

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prepoststep_pdaf(int* step, int* dim_p, int* dim_ens,
    int* dim_ens_l, int* dim_obs_p, double* state_p, double* uinv,
    double* ens_p, int* flag) noexcept with gil:
    """Process ensemble before or after DA.

    Parameters
    ----------
    step : int
            current time step
            (negative for call before analysis/preprocessing)
    dim_p : int
            PE-local state vector dimension
    dim_ens : int
            number of ensemble members
    dim_ens_l : int
            number of ensemble members run serially
            on each model task
    dim_obs_p : int
            PE-local dimension of observation vector
    state_p : ndarray[np.float64, ndim=1]
            pe-local forecast/analysis state
            (the array 'state_p' is generally not
            initialised in the case of ESTKF/ETKF/EnKF/SEIK,
            so it can be used freely here.)
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            Inverse of the transformation matrix in ETKF and ESKTF;
            inverse of matrix formed by right singular vectors of error
            covariance matrix of ensemble perturbations in SEIK/SEEK.
            not used in EnKF.
            Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            pe-local forecast/analysis state
            (the array 'state_p' is generally not
            initialised in the case of ESTKF/ETKF/EnKF/SEIK,
            so it can be used freely here.)
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            Inverse of the transformation matrix in ETKF and ESKTF;
            inverse of matrix formed by right singular vectors of error
            covariance matrix of ensemble perturbations in SEIK/SEEK.
            not used in EnKF.
            Array shape: (dim_ens-1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    """
    cdef size_t uinv_len = max(dim_ens[0]-1, 1)
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1,:] uinv_np = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")

    state_p_np,uinv_np,ens_p_np = (<object>prepoststep_pdaf)(step[0],
                                                             dim_p[0],
                                                             dim_ens[0],
                                                             dim_ens_l[0],
                                                             dim_obs_p[0],
                                                             state_p_np.base,
                                                             uinv_np.base,
                                                             ens_p_np.base,
                                                             flag[0])

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_dim_obs_pdaf(int* step, int* dim_obs_p) noexcept with gil:
    """Determine the size of the vector of observations

    The primary purpose of this function is to
    obtain the dimension of the observation vector.
    In OMI, in this function, one also sets the properties
    of `obs_f`, read the observation vector from
    files, setting the observation error variance
    when diagonal observation error covariance matrix
    is used. The `pyPDAF.PDAF.omi_gather_obs` function
    is also called here.

    Furthermore, in this user-supplied function, one also sets the interpolation
    coefficients used by observation operators.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.

    Returns
    -------
    dim_obs_p : int
        Dimension of the observation vector.
    """
    dim_obs_p[0] = (<object>init_dim_obs_pdaf)(step[0], dim_obs_p[0])



cdef void c__init_dim_obs_f_pdaf(int* step, int* dim_obs_p) noexcept with gil:
    """Determine the size of the vector of observations

    The primary purpose of this function is to
    obtain the dimension of the observation vector.
    In OMI, in this function, one also sets the properties
    of `obs_f`, read the observation vector from
    files, setting the observation error variance
    when diagonal observation error covariance matrix
    is used. The `pyPDAF.PDAF.omi_gather_obs` function
    is also called here.

    Furthermore, in this user-supplied function, one also sets the interpolation
    coefficients used by observation operators.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.

    Returns
    -------
    dim_obs_p : int
        Dimension of the observation vector.
    """
    dim_obs_p[0] = (<object>init_dim_obs_f_pdaf)(step[0], dim_obs_p[0])



cdef void c__init_obs_pdaf(int* step, int* dim_obs_p,
    double* observation_p) noexcept with gil:
    """Provide the observation vector for the current time step.

    The routine is used with all global filters.

    Parameters
    ----------
    step : int
        Current time step.
    dim_obs_p : int
        Dimension of the observation vector.
    observation_p : ndarray[np.float64, ndim=1]
        Observation vector.

    Returns
    -------
    observation_p : ndarray[np.float64, ndim=1]
        Filled observation vector.
    """
    cdef double[::1] observation_p_np = np.asarray(<double[:dim_obs_p[0]:1]> observation_p, order="F")

    observation_p_np = (<object>init_obs_pdaf)(step[0], dim_obs_p[0],
                                               observation_p_np.base)

    cdef double[::1] observation_p_new
    if observation_p != &observation_p_np[0]:
        observation_p_new = np.asarray(<double[:dim_obs_p[0]:1]> observation_p, order="F")
        observation_p_new[...] = observation_p_np
        warnings.warn("The memory address of observation_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obs_covar_pdaf(int* step, int* dim_obs, int* dim_obs_p,
    double* covar, double* obs_p, bint* isdiag) noexcept with gil:
    cdef double[::1,:] covar_np = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> covar, order="F")
    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]:1]> obs_p, order="F")

    covar_np,isdiag[0] = (<object>init_obs_covar_pdaf)(step[0], dim_obs[0],
                                                       dim_obs_p[0],
                                                       covar_np.base,
                                                       obs_p_np.base, isdiag[0])

    cdef double[::1,:] covar_new
    if covar != &covar_np[0,0]:
        covar_new = np.asarray(<double[:dim_obs_p[0]:1,:dim_obs_p[0]]> covar, order="F")
        covar_new[...] = covar_np
        warnings.warn("The memory address of covar is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obsvar_pdaf(int* step, int* dim_obs_p, double* obs_p,
    double* meanvar) noexcept with gil:
    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]:1]> obs_p, order="F")

    meanvar[0] = (<object>init_obsvar_pdaf)(step[0], dim_obs_p[0],
                                            obs_p_np.base, meanvar[0])



cdef void c__init_obsvars_pdaf(int* step, int* dim_obs_f,
    double* var_f) noexcept with gil:
    cdef double[::1] var_f_np = np.asarray(<double[:dim_obs_f[0]:1]> var_f, order="F")

    var_f_np = (<object>init_obsvars_pdaf)(step[0], dim_obs_f[0], var_f_np.base)

    cdef double[::1] var_f_new
    if var_f != &var_f_np[0]:
        var_f_new = np.asarray(<double[:dim_obs_f[0]:1]> var_f, order="F")
        var_f_new[...] = var_f_np
        warnings.warn("The memory address of var_f is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prodrinva_pdaf(int* step, int* dim_obs_p, int* rank,
    double* obs_p, double* a_p, double* c_p) noexcept with gil:
    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]:1]> obs_p, order="F")
    cdef double[::1,:] a_p_np = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> a_p, order="F")
    cdef double[::1,:] c_p_np = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> c_p, order="F")

    c_p_np = (<object>prodrinva_pdaf)(step[0], dim_obs_p[0], rank[0],
                                      obs_p_np.base, a_p_np.base, c_p_np.base)

    cdef double[::1,:] c_p_new
    if c_p != &c_p_np[0,0]:
        c_p_new = np.asarray(<double[:dim_obs_p[0]:1,:rank[0]]> c_p, order="F")
        c_p_new[...] = c_p_np
        warnings.warn("The memory address of c_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_pdaf(int* step, int* dim_p, int* dim_obs_p,
    double* state_p, double* m_state_p) noexcept with gil:
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")

    m_state_p_np = (<object>obs_op_pdaf)(step[0], dim_p[0], dim_obs_p[0],
                                         state_p_np.base, m_state_p_np.base)

    cdef double[::1] m_state_p_new
    if m_state_p != &m_state_p_np[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")
        m_state_p_new[...] = m_state_p_np
        warnings.warn("The memory address of m_state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_f_pdaf(int* step, int* dim_p, int* dim_obs_p,
    double* state_p, double* m_state_p) noexcept with gil:
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")

    m_state_p_np = (<object>obs_op_f_pdaf)(step[0], dim_p[0], dim_obs_p[0],
                                           state_p_np.base, m_state_p_np.base)

    cdef double[::1] m_state_p_new
    if m_state_p != &m_state_p_np[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")
        m_state_p_new[...] = m_state_p_np
        warnings.warn("The memory address of m_state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__g2l_obs_pdaf(int* domain_p, int* step, int* dim_obs_f,
    int* dim_obs_l, int* mstate_f, int* dim_p, int* mstate_l,
    int* dim_l) noexcept with gil:
    cdef int[::1] mstate_f_np = np.asarray(<int[:dim_p[0]:1]> mstate_f, order="F")
    cdef int[::1] mstate_l_np = np.asarray(<int[:dim_l[0]:1]> mstate_l, order="F")

    mstate_l_np = (<object>g2l_obs_pdaf)(domain_p[0], step[0],
                                         dim_obs_f[0], dim_obs_l[0],
                                         mstate_f_np.base, dim_p[0],
                                         mstate_l_np.base, dim_l[0])

    cdef int[::1] mstate_l_new
    if mstate_l != &mstate_l_np[0]:
        mstate_l_new = np.asarray(<int[:dim_l[0]:1]> mstate_l, order="F")
        mstate_l_new[...] = mstate_l_np
        warnings.warn("The memory address of mstate_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__g2l_state_pdaf(int* step, int* domain_p, int* dim_p,
    double* state_p, int* dim_l, double* state_l) noexcept with gil:
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1] state_l_np = np.asarray(<double[:dim_l[0]:1]> state_l, order="F")

    state_l_np = (<object>g2l_state_pdaf)(step[0], domain_p[0], dim_p[0],
                                          state_p_np.base, dim_l[0],
                                          state_l_np.base)

    cdef double[::1] state_l_new
    if state_l != &state_l_np[0]:
        state_l_new = np.asarray(<double[:dim_l[0]:1]> state_l, order="F")
        state_l_new[...] = state_l_np
        warnings.warn("The memory address of state_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_dim_l_pdaf(int* step, int* domain_p,
    int* dim_l) noexcept with gil:

    dim_l[0] = (<object>init_dim_l_pdaf)(step[0], domain_p[0], dim_l[0])



cdef void c__init_dim_obs_l_pdaf(int* domain_p, int* step, int* dim_obs_f,
    int* dim_obs_l) noexcept with gil:

    dim_obs_l[0] = (<object>init_dim_obs_l_pdaf)(domain_p[0], step[0],
                                                 dim_obs_f[0], dim_obs_l[0])



cdef void c__init_n_domains_p_pdaf(int* step,
    int* n_domains_p) noexcept with gil:

    n_domains_p[0] = (<object>init_n_domains_p_pdaf)(step[0], n_domains_p[0])



cdef void c__init_obs_f_pdaf(int* step, int* dim_obs_f,
    double* observation_f) noexcept with gil:
    cdef double[::1] observation_f_np = np.asarray(<double[:dim_obs_f[0]:1]> observation_f, order="F")

    observation_f_np = (<object>init_obs_f_pdaf)(step[0], dim_obs_f[0],
                                                 observation_f_np.base)

    cdef double[::1] observation_f_new
    if observation_f != &observation_f_np[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]:1]> observation_f, order="F")
        observation_f_new[...] = observation_f_np
        warnings.warn("The memory address of observation_f is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obs_l_pdaf(int* domain_p, int* step, int* dim_obs_l,
    double* observation_l) noexcept with gil:
    cdef double[::1] observation_l_np = np.asarray(<double[:dim_obs_l[0]:1]> observation_l, order="F")

    observation_l_np = (<object>init_obs_l_pdaf)(domain_p[0], step[0],
                                                 dim_obs_l[0],
                                                 observation_l_np.base)

    cdef double[::1] observation_l_new
    if observation_l != &observation_l_np[0]:
        observation_l_new = np.asarray(<double[:dim_obs_l[0]:1]> observation_l, order="F")
        observation_l_new[...] = observation_l_np
        warnings.warn("The memory address of observation_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__init_obsvar_l_pdaf(int* domain_p, int* step, int* dim_obs_l,
    double* obs_l, int* dim_obs_p, double* meanvar_l) noexcept with gil:
    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_p[0]:1]> obs_l, order="F")

    meanvar_l[0] = (<object>init_obsvar_l_pdaf)(domain_p[0], step[0],
                                                dim_obs_l[0],
                                                obs_l_np.base,
                                                dim_obs_p[0], meanvar_l[0])



cdef void c__init_obserr_f_pdaf(int* step, int* dim_obs_f, double* obs_f,
    double* obserr_f) noexcept with gil:
    cdef double[::1] obs_f_np = np.asarray(<double[:dim_obs_f[0]:1]> obs_f, order="F")
    cdef double[::1] obserr_f_np = np.asarray(<double[:dim_obs_f[0]:1]> obserr_f, order="F")

    obserr_f_np = (<object>init_obserr_f_pdaf)(step[0], dim_obs_f[0],
                                               obs_f_np.base, obserr_f_np.base)

    cdef double[::1] obserr_f_new
    if obserr_f != &obserr_f_np[0]:
        obserr_f_new = np.asarray(<double[:dim_obs_f[0]:1]> obserr_f, order="F")
        obserr_f_new[...] = obserr_f_np
        warnings.warn("The memory address of obserr_f is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__l2g_state_pdaf(int* step, int* domain_p, int* dim_l,
    double* state_l, int* dim_p, double* state_p) noexcept with gil:
    cdef double[::1] state_l_np = np.asarray(<double[:dim_l[0]:1]> state_l, order="F")
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")

    state_p_np = (<object>l2g_state_pdaf)(step[0], domain_p[0], dim_l[0],
                                          state_l_np.base, dim_p[0],
                                          state_p_np.base)

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prodrinva_l_pdaf(int* domain_p, int* step, int* dim_obs_l,
    int* rank, double* obs_l, double* a_l, double* c_l) noexcept with gil:
    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]:1]> obs_l, order="F")
    cdef double[::1,:] a_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> a_l, order="F")
    cdef double[::1,:] c_l_np = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> c_l, order="F")

    a_l_np,c_l_np = (<object>prodrinva_l_pdaf)(domain_p[0], step[0],
                                               dim_obs_l[0], rank[0],
                                               obs_l_np.base, a_l_np.base,
                                               c_l_np.base)

    cdef double[::1,:] a_l_new
    if a_l != &a_l_np[0,0]:
        a_l_new = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> a_l, order="F")
        a_l_new[...] = a_l_np
        warnings.warn("The memory address of a_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] c_l_new
    if c_l != &c_l_np[0,0]:
        c_l_new = np.asarray(<double[:dim_obs_l[0]:1,:rank[0]]> c_l, order="F")
        c_l_new[...] = c_l_np
        warnings.warn("The memory address of c_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__localize_covar_pdaf(int* dim_p, int* dim_obs, double* hp_p,
    double* hph) noexcept with gil:
    cdef double[::1,:] hp_p_np = np.asarray(<double[:dim_obs[0]:1,:dim_p[0]]> hp_p, order="F")
    cdef double[::1,:] hph_np = np.asarray(<double[:dim_obs[0]:1,:dim_obs[0]]> hph, order="F")

    hp_p_np,hph_np = (<object>localize_covar_pdaf)(dim_p[0], dim_obs[0],
                                                   hp_p_np.base, hph_np.base)

    cdef double[::1,:] hp_p_new
    if hp_p != &hp_p_np[0,0]:
        hp_p_new = np.asarray(<double[:dim_obs[0]:1,:dim_p[0]]> hp_p, order="F")
        hp_p_new[...] = hp_p_np
        warnings.warn("The memory address of hp_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] hph_new
    if hph != &hph_np[0,0]:
        hph_new = np.asarray(<double[:dim_obs[0]:1,:dim_obs[0]]> hph, order="F")
        hph_new[...] = hph_np
        warnings.warn("The memory address of hph is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__localize_covar_serial_pdaf(int* iobs, int* dim_p,
    int* dim_obs, double* hp_p, double* hxy_p) noexcept with gil:
    cdef double[::1] hp_p_np = np.asarray(<double[:dim_p[0]:1]> hp_p, order="F")
    cdef double[::1] hxy_p_np = np.asarray(<double[:dim_obs[0]:1]> hxy_p, order="F")

    hp_p_np,hxy_p_np = (<object>localize_covar_serial_pdaf)(iobs[0],
                                                            dim_p[0],
                                                            dim_obs[0],
                                                            hp_p_np.base,
                                                            hxy_p_np.base)

    cdef double[::1] hp_p_new
    if hp_p != &hp_p_np[0]:
        hp_p_new = np.asarray(<double[:dim_p[0]:1]> hp_p, order="F")
        hp_p_new[...] = hp_p_np
        warnings.warn("The memory address of hp_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1] hxy_p_new
    if hxy_p != &hxy_p_np[0]:
        hxy_p_new = np.asarray(<double[:dim_obs[0]:1]> hxy_p, order="F")
        hxy_p_new[...] = hxy_p_np
        warnings.warn("The memory address of hxy_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__likelihood_pdaf(int* step, int* dim_obs_p, double* obs_p,
    double* resid, double* likely) noexcept with gil:
    cdef double[::1] obs_p_np = np.asarray(<double[:dim_obs_p[0]:1]> obs_p, order="F")
    cdef double[::1] resid_np = np.asarray(<double[:dim_obs_p[0]:1]> resid, order="F")

    likely[0] = (<object>likelihood_pdaf)(step[0], dim_obs_p[0],
                                          obs_p_np.base, resid_np.base,
                                          likely[0])



cdef void c__likelihood_l_pdaf(int* domain_p, int* step, int* dim_obs_l,
    double* obs_l, double* resid_l, double* likely_l) noexcept with gil:
    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]:1]> obs_l, order="F")
    cdef double[::1] resid_l_np = np.asarray(<double[:dim_obs_l[0]:1]> resid_l, order="F")

    resid_l_np,likely_l[0] = (<object>likelihood_l_pdaf)(domain_p[0],
                                                         step[0],
                                                         dim_obs_l[0],
                                                         obs_l_np.base,
                                                         resid_l_np.base,
                                                         likely_l[0])

    cdef double[::1] resid_l_new
    if resid_l != &resid_l_np[0]:
        resid_l_new = np.asarray(<double[:dim_obs_l[0]:1]> resid_l, order="F")
        resid_l_new[...] = resid_l_np
        warnings.warn("The memory address of resid_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__get_obs_f_pdaf(int* step, int* dim_obs_f,
    double* observation_f) noexcept with gil:
    cdef double[::1] observation_f_np = np.asarray(<double[:dim_obs_f[0]:1]> observation_f, order="F")

    observation_f_np = (<object>get_obs_f_pdaf)(step[0], dim_obs_f[0],
                                                observation_f_np.base)

    cdef double[::1] observation_f_new
    if observation_f != &observation_f_np[0]:
        observation_f_new = np.asarray(<double[:dim_obs_f[0]:1]> observation_f, order="F")
        observation_f_new[...] = observation_f_np
        warnings.warn("The memory address of observation_f is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_adj_ens_pdaf(int* iter, int* dim_p, int* dim_ens,
    int* dim_cv_ens_p, double* ens_p, double* vcv_p,
    double* cv_p) noexcept with gil:
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")
    cdef double[::1] vcv_p_np = np.asarray(<double[:dim_p[0]:1]> vcv_p, order="F")
    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cv_ens_p[0]:1]> cv_p, order="F")

    cv_p_np = (<object>cvt_adj_ens_pdaf)(iter[0], dim_p[0], dim_ens[0],
                                         dim_cv_ens_p[0], ens_p_np.base,
                                         vcv_p_np.base, cv_p_np.base)

    cdef double[::1] cv_p_new
    if cv_p != &cv_p_np[0]:
        cv_p_new = np.asarray(<double[:dim_cv_ens_p[0]:1]> cv_p, order="F")
        cv_p_new[...] = cv_p_np
        warnings.warn("The memory address of cv_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_adj_pdaf(int* iter, int* dim_p, int* dim_cvec,
    double* vcv_p, double* cv_p) noexcept with gil:
    cdef double[::1] vcv_p_np = np.asarray(<double[:dim_p[0]:1]> vcv_p, order="F")
    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cvec[0]:1]> cv_p, order="F")

    cv_p_np = (<object>cvt_adj_pdaf)(iter[0], dim_p[0], dim_cvec[0],
                                     vcv_p_np.base, cv_p_np.base)

    cdef double[::1] cv_p_new
    if cv_p != &cv_p_np[0]:
        cv_p_new = np.asarray(<double[:dim_cvec[0]:1]> cv_p, order="F")
        cv_p_new[...] = cv_p_np
        warnings.warn("The memory address of cv_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_pdaf(int* iter, int* dim_p, int* dim_cvec, double* cv_p,
    double* vv_p) noexcept with gil:
    cdef double[::1] cv_p_np = np.asarray(<double[:dim_cvec[0]:1]> cv_p, order="F")
    cdef double[::1] vv_p_np = np.asarray(<double[:dim_p[0]:1]> vv_p, order="F")

    vv_p_np = (<object>cvt_pdaf)(iter[0], dim_p[0], dim_cvec[0],
                                 cv_p_np.base, vv_p_np.base)

    cdef double[::1] vv_p_new
    if vv_p != &vv_p_np[0]:
        vv_p_new = np.asarray(<double[:dim_p[0]:1]> vv_p, order="F")
        vv_p_new[...] = vv_p_np
        warnings.warn("The memory address of vv_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__cvt_ens_pdaf(int* iter, int* dim_p, int* dim_ens,
    int* dim_cvec_ens, double* ens_p, double* v_p,
    double* vv_p) noexcept with gil:
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")
    cdef double[::1] v_p_np = np.asarray(<double[:dim_cvec_ens[0]:1]> v_p, order="F")
    cdef double[::1] vv_p_np = np.asarray(<double[:dim_p[0]:1]> vv_p, order="F")

    vv_p_np = (<object>cvt_ens_pdaf)(iter[0], dim_p[0], dim_ens[0],
                                     dim_cvec_ens[0], ens_p_np.base,
                                     v_p_np.base, vv_p_np.base)

    cdef double[::1] vv_p_new
    if vv_p != &vv_p_np[0]:
        vv_p_new = np.asarray(<double[:dim_p[0]:1]> vv_p, order="F")
        vv_p_new[...] = vv_p_np
        warnings.warn("The memory address of vv_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_adj_pdaf(int* step, int* dim_p, int* dim_obs_p,
    double* m_state_p, double* state_p) noexcept with gil:
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")

    state_p_np = (<object>obs_op_adj_pdaf)(step[0], dim_p[0], dim_obs_p[0],
                                           m_state_p_np.base, state_p_np.base)

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__obs_op_lin_pdaf(int* step, int* dim_p, int* dim_obs_p,
    double* state_p, double* m_state_p) noexcept with gil:
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1] m_state_p_np = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")

    m_state_p_np = (<object>obs_op_lin_pdaf)(step[0], dim_p[0],
                                             dim_obs_p[0], state_p_np.base,
                                             m_state_p_np.base)

    cdef double[::1] m_state_p_new
    if m_state_p != &m_state_p_np[0]:
        m_state_p_new = np.asarray(<double[:dim_obs_p[0]:1]> m_state_p, order="F")
        m_state_p_new[...] = m_state_p_np
        warnings.warn("The memory address of m_state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__dist_stateinc_pdaf(int* dim_p, double* state_inc_p,
    int* first, int* steps) noexcept with gil:
    cdef double[::1] state_inc_p_np = np.asarray(<double[:dim_p[0]:1]> state_inc_p, order="F")

    (<object>dist_stateinc_pdaf)(dim_p[0], state_inc_p_np.base, first[0],
                                 steps[0])



cdef void c__likelihood_hyb_l_pdaf(int* domain_p, int* step,
    int* dim_obs_l, double* obs_l, double* resid_l, double* gamma,
    double* likely_l) noexcept with gil:
    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]:1]> obs_l, order="F")
    cdef double[::1] resid_l_np = np.asarray(<double[:dim_obs_l[0]:1]> resid_l, order="F")

    resid_l_np,likely_l[0] = (<object>likelihood_hyb_l_pdaf)(domain_p[0],
                                                             step[0],
                                                             dim_obs_l[0],
                                                             obs_l_np.base,
                                                             resid_l_np.base,
                                                             gamma[0],
                                                             likely_l[0])

    cdef double[::1] resid_l_new
    if resid_l != &resid_l_np[0]:
        resid_l_new = np.asarray(<double[:dim_obs_l[0]:1]> resid_l, order="F")
        resid_l_new[...] = resid_l_np
        warnings.warn("The memory address of resid_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


cdef void c__prodrinva_hyb_l_pdaf(int* domain_p, int* step, int* dim_obs_l,
    int* dim_ens, double* obs_l, double* gamma, double* a_l,
    double* c_l) noexcept with gil:
    cdef double[::1] obs_l_np = np.asarray(<double[:dim_obs_l[0]:1]> obs_l, order="F")
    cdef double[::1,:] a_l_np = np.asarray(<double[:dim_obs_l[0]:1,:dim_ens[0]]> a_l, order="F")
    cdef double[::1,:] c_l_np = np.asarray(<double[:dim_obs_l[0]:1,:dim_ens[0]]> c_l, order="F")

    a_l_np,c_l_np = (<object>prodrinva_hyb_l_pdaf)(domain_p[0], step[0],
                                                   dim_obs_l[0],
                                                   dim_ens[0],
                                                   obs_l_np.base, gamma[0],
                                                   a_l_np.base, c_l_np.base)

    cdef double[::1,:] a_l_new
    if a_l != &a_l_np[0,0]:
        a_l_new = np.asarray(<double[:dim_obs_l[0]:1,:dim_ens[0]]> a_l, order="F")
        a_l_new[...] = a_l_np
        warnings.warn("The memory address of a_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] c_l_new
    if c_l != &c_l_np[0,0]:
        c_l_new = np.asarray(<double[:dim_obs_l[0]:1,:dim_ens[0]]> c_l, order="F")
        c_l_new[...] = c_l_np
        warnings.warn("The memory address of c_l is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)


