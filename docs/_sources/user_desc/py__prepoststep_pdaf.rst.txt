py__prepoststep_pdaf
====================

.. py:function:: py__prepoststep_pdaf(step: int, dim_p: int, dim_ens: int, dim_ens_l: int, dim_obs_p: int, state_p: np.ndarray, uinv: np.ndarray, ens_p: np.ndarray, flag: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]

    Process ensemble before or after DA.

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
