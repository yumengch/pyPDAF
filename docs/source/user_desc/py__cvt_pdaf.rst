py__cvt_pdaf
============

.. py:function:: py__cvt_pdaf(iter: int, dim_p: int, dim_cvec: int, cv_p: np.ndarray, vv_p: np.ndarray) -> np.ndarray

    The control variable transformation.

    Here, this function performs :math:`\mathbf{U} \mathbf{v}` with
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
    cv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{v}`. shape: (dim_cvec, )
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )

    Returns
    -------
    vv_p: np.ndarray[np.float, dim=1]
        :math:`\mathbf{U} \mathbf{v}`. shape: (dim_p, )
