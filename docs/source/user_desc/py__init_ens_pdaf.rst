py__init_ens_pdaf
=================

.. py:function:: py__init_ens_pdaf(filtertype: int, dim_p: int, dim_ens: int, state_p: np.ndarray, uinv: np.ndarray, ens_p: np.ndarray, flag: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]

    Fill the ensemble array that is provided by PDAF with an initial ensemble of model states.

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
