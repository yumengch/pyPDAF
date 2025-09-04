import numpy as np
cimport numpy as cnp

def g2l_cb(int  step, int  domain_p, int  dim_p, double [::1] state_p,
    int  dim_l):
    r"""Project a global to a local state vector for the localized filters.

    This is the full callback function to be used internally.
    The mapping is done using the index vector id_lstate_in_pstate that is
    initialised in `pyPDAF.PDAFlocal.set_indices`.

    Parameters
    ----------
    step : int
        Current time step
    domain_p : int
        Current local analysis domain
    dim_p : int
        PE-local full state dimension
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local full state vector
        The array dimension `dim_p` is PE-local full state dimension
    dim_l : int
        Local state dimension

    Returns
    -------
    state_l : ndarray[tuple[dim_l, ...], np.float64]
         State vector on local analysis domain

        The array dimension `dim_l` is Local state dimension
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_l_np = np.zeros((dim_l), dtype=np.float64, order="F")
    cdef double [::1] state_l = state_l_np
    with nogil:
        c__pdaflocal_g2l_cb(&step, &domain_p, &dim_p, &state_p[0], &dim_l,
                            &state_l[0])

    return state_l_np


def l2g_cb(int  step, int  domain_p, int  dim_l, double [::1] state_l,
    int  dim_p, double [::1] state_p):
    r"""Initialise elements of a global state vector from a local state vector.

    This is the full callback function to be used internally.
    The mapping is done using the index vector `id_lstate_in_pstate` that is
    initialised in :func:`pyPDAF.PDAFlocal.set_indices`.

    To exclude any element of the local state vector from the initialisationone
    can set the corresponding index value to 0.

    Parameters
    ----------
    step : int
        Current time step
    domain_p : int
        Current local analysis domain
    dim_l : int
        Local state dimension
    state_l : ndarray[tuple[dim_l, ...], np.float64]
        State vector on local analysis domain
        The array dimension `dim_l` is Local state dimension
    dim_p : int
        PE-local full state dimension
    state_p : ndarray[tuple[dim_p, ...], np.float64]
        PE-local full state vector
        The array dimension `dim_p` is PE-local full state dimension

    Returns
    -------
    state_p : ndarray[tuple[dim_p, ...], np.float64]
         PE-local full state vector

        The array dimension `dim_p` is PE-local full state dimension
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] state_p_np = np.asarray(state_p, dtype=np.float64, order="F")
    with nogil:
        c__pdaflocal_l2g_cb(&step, &domain_p, &dim_l, &state_l[0], &dim_p,
                            &state_p[0])

    return state_p_np