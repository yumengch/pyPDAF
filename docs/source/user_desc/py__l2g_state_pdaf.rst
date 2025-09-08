.. function:: py__l2g_state_pdaf

    Assign local state vector to process-local global state vector.

    Parameters
    ----------
    step: int
        Current time step
    domain_p: int
        Current local domain index
    dim_l: int
        Local analysis domain state vector dimension.
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector.
    dim_p: int
        Process-local global state vector dimension.
    state_p: np.ndarray[np.float, dim=1]
        Process-local global state vector.

    Returns
    -------
    state_p: np.ndarray[np.float, dim=1]
        Process-local global state vector.
