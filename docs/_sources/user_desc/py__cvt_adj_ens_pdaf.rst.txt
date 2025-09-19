py__cvt_adj_ens_pdaf
=====================

.. py:function:: py__cvt_adj_ens_pdaf(iter: int, dim_p: int, dim_ens: int, dim_cv_ens_p:int, ens: np.ndarray, vcv_p: np.ndarray, cv_p: np.ndarray) -> np.ndarray:

    The adjoint control variable transformation involving ensembles.

    Here, this function performs
    :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. Here, the input vector can be
    :math:`\mathbf{v}_h = \mathbf{H}^\mathrm{T}\mathbf{R}^{-1}\delta \mathbf{x}`
    with :math:`\delta \mathbf{x}` the innovation vector.

    The control vector transform is given by :math:`\mathbf{U}\mathbf{v}` with
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
    vcv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}_h`. shape: (dim_p, )
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cv_ens_p, )

    Returns
    -------
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cv_ens_p, )
