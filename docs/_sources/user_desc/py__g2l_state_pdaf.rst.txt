py__g2l_state_pdaf
==================

.. py:function:: py__g2l_state_pdaf(step: int, domain_p: int, dim_p: int, state_p: np.ndarray, dim_l: int, state_l: np.ndarray) -> np.ndarray

    Get local state vector.

    Get the state vector for analysis local domain.

    Parameters
    ----------
    step: int
        Current time step
    domain_p: int
        Current local domain index
    dim_p: int
        Process-local state vector dimension.
    state_p: np.ndarray[np.float, dim=1]
        Process-local state vector.
    dim_l: int
        Local analysis domain state vector dimension.
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector. Shape: (dim_l,)

    Returns
    -------
    state_l: np.ndarray[np.float, dim=1]
        Local analysis domain state vector. Shape: (dim_l,)
