.. function:: py__add_obs_err_pdaf
   Add the observation error covariance matrix to the matrix ``C``.

   The input matrix is the projection of the ensemble covariance
   matrix onto the observation space that is computed during the
   analysis step of the stochastic EnKF, i.e. :math:`HPH^T`.
   The function returns :math:`HPH^T + R`.

   The operation is for the global observation space, thus it is
   independent of whether the filter is executed with or without
   parallelization.

   Parameters
   ----------
   step : int
      Current time step
   dim_obs_p : int
      Size of observation vector
   C_p : ndarray[np.float64, ndim=2]
      Matrix to which the observation error covariance matrix is added
      shape: (dim_obs_p, dim_obs_p)

   Returns
   -------
   C_p : ndarray[np.float64, ndim=2]
      Matrix with added obs error covariance, i.e., HPH.T + R
      shape: (dim_obs_p, dim_obs_p)