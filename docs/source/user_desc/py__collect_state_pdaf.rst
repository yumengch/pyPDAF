py__collect_state_pdaf
======================

.. py:function:: py__collect_state_pdaf(dim_p: int, state_p: np.ndarray) -> np.ndarray

    Collect state vector from model/any arrays to pdaf arrays

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
