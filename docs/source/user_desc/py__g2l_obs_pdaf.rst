py__g2l_obs_pdaf
================

.. py:function:: py__g2l_obs_pdaf(domain_p: int, step: int, dim_obs_f: int, dim_obs_l: int, mstate_f: np.ndarray, mstate_l: np.ndarray) -> np.ndarray

    Convert global observed state vector to local vector.

    This is used by domain localisation methods. In these methods, each local
    domain has their own observation vector and observed state vector.

    Parameters
    ----------
    domain_p:int
        Current local domain index
    step: int
        Current time step
    dim_obs_f: int
        Global observation vector dimension.
    dim_obs_l: int
        Local observation vector dimension.
    mstate_f: np.ndarray[np.float64, dim=1]
        Global observed state vector. shape: (dim_obs_f,)
    mstate_l: np.ndarray[np.float64, dim=1]
        Local observed state vector. shape: (dim_obs_l,)

    Returns
    -------
    mstate_l: np.ndarray[np.float64, dim=1]
        Local observed state vector. shape: (dim_obs_l,)
