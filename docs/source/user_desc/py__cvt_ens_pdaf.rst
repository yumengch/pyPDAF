py__cvt_ens_pdaf
================

.. py:function:: py__cvt_ens_pdaf(iter: int, dim_p: int, dim_ens: int, dim_cv_ens_p: int, v_p: np.ndarray, vv_p: np.ndarray) -> np.ndarray

    The  control variable transformation involving ensembles.

    Here, this function performs :math:`\mathbf{U} \mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DEnVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_ens: int
        Dimension of the ensemble.
    dim_cv_ens_p: int
        Dimension of the control vector.
    ens_p: np.ndarray[np.float, dim=2]
        Ensemble matrix. shape: (dim_p, dim_ens)
    v_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}`. shape: (dim_cv_ens_p, )
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )

    Returns
    -------
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )
