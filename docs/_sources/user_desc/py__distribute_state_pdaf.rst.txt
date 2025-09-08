.. function:: py__distribute_state_pdaf
    Distribute a state vector from pdaf to the model/any arrays

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