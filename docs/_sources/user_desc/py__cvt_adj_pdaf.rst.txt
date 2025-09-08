.. function:: py__cvt_adj_pdaf

    The adjoint control variable transformation.

    Here, this function performs
    :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. Here, the input vector can be
    :math:`\mathbf{v}_h = \mathbf{H}^\mathrm{T}\mathbf{R}^{-1}\delta \mathbf{x}`
    with :math:`\delta \mathbf{x}` the innovation vector.

    The control vector transform is given by :math:`\mathbf{U}\mathbf{v}` with
    :math:`\mathbf{v}` the control vector and :math:`\mathbf{U}` the transformation
    matrix, which can be :math:`\mathbf{B}^\frac{1}{2}`.

    This function is used in 3DVar and hybrid 3DVar.

    Parameters
    ----------
    iter: int
        Current optimisation iteration number.
    dim_p: int
        Dimension of the state vector.
    dim_cvec: int
        Dimension of the control vector.
    vcv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}_h`. shape: (dim_p, )
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cvec, )

    Returns
    -------
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U}^\mathrm{T} \mathbf{v}_h`. shape: (dim_cvec, )
