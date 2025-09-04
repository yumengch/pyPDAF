import numpy as np
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb


def assimilate(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Online ensemble filters and smoothers except for 3DVars for a single DA step
    using diagnoal observation error covariance matrix.

    Here, this function call is used for
    global stochastic EnKF [1]_, E(S)TKF [2]_, EAKF, EnSRF,
    SEEK [2]_, SEIK [2]_, NETF [3]_, and particle filter [4]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. core DA algorithm
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    References
    ----------
    .. [1] Evensen, G. (1994),
           Sequential data assimilation with
           a nonlinear quasi-geostrophic model
           using Monte Carlo methods to forecast error statistics,
           J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572.
    .. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1
    .. [3] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.
    .. [4] Van Leeuwen, P. J., Künsch, H. R., Nerger, L.,
           Potthast, R., & Reich, S. (2019).
           Particle filters for high‐dimensional geoscience applications:
           A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365.


    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate(pdaf_cb.c__collect_state_pdaf,
                            pdaf_cb.c__distribute_state_pdaf,
                            pdaf_cb.c__init_dim_obs_pdaf,
                            pdaf_cb.c__obs_op_pdaf,
                            pdaf_cb.c__init_n_domains_p_pdaf,
                            pdaf_cb.c__init_dim_l_pdaf,
                            pdaf_cb.c__init_dim_obs_l_pdaf,
                            pdaf_cb.c__prepoststep_pdaf,
                            pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assim_offline(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf, int  outflag):
    r"""Offline ensemble filters and smoothers except for 3DVars for a single DA step
    using diagnoal observation error covariance matrix.

    Here, this function call is used for
    global stochastic EnKF [1]_, E(S)TKF [2]_, EAKF, EnSRF,
    SEEK [2]_, SEIK [2]_, NETF [3]_, and particle filter [4]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_n_domains_p_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. core DA algorithm
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    References
    ----------
    .. [1] Evensen, G. (1994),
           Sequential data assimilation with
           a nonlinear quasi-geostrophic model
           using Monte Carlo methods to forecast error statistics,
           J. Geophys. Res., 99(C5), 10143–10162, doi:10.1029/94JC00572.
    .. [2] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1
    .. [3] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.
    .. [4] Van Leeuwen, P. J., Künsch, H. R., Nerger, L.,
           Potthast, R., & Reich, S. (2019).
           Particle filters for high‐dimensional geoscience applications:
           A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365.


    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_assim_offline(pdaf_cb.c__init_dim_obs_pdaf,
                               pdaf_cb.c__obs_op_pdaf,
                               pdaf_cb.c__init_n_domains_p_pdaf,
                               pdaf_cb.c__init_dim_l_pdaf,
                               pdaf_cb.c__init_dim_obs_l_pdaf,
                               pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assimilate_3dvar_all(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Online assimilation for all types of 3DVar DA for a single DA step
    using diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate_3dvar_nondiagR`,
    :func:`pyPDAF.PDAF3.assimilate_en3dvar_lestkf_nondiagR`,
    :func:`pyPDAF.PDAF3.assimilate_en3dvar_estkf_nondiagR`,
    :func:`pyPDAF.PDAF3.assimilate_hyb3dvar_lestkf_nondiagR`,
    :func:`pyPDAF.PDAF3.assimilate_hyb3dvar_estkf_nondiagR`
    for non-diagonal observation error covariance matrix.

    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable transformation.
    This is a deterministic filtering scheme
    so no ensemble and parallelisation is needed.
    This function should be called at each model time step.

    For parametrised 3DVar, user-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__obs_op_adj_pdaf
            4. py__cvt_adj_pdaf
            5. core DA algorithm
        6. py__cvt_pdaf
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    For 3DEnVar, user-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__obs_op_adj_pdaf
            4. py__cvt_adj_ens_pdaf
            5. core DA algorithm
        6. py__cvt_ens_pdaf
        7. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf


    For hybrid 3DVar, user-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. py__cvt_adj_ens_pdaf
            7. core DA algorithm
        6. py__cvt_pdaf
        7. py__cvt_ens_pdaf
        8. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. core DA algorithm
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf


    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_3dvar_all(pdaf_cb.c__collect_state_pdaf,
                                      pdaf_cb.c__distribute_state_pdaf,
                                      pdaf_cb.c__init_dim_obs_pdaf,
                                      pdaf_cb.c__obs_op_pdaf,
                                      pdaf_cb.c__cvt_ens_pdaf,
                                      pdaf_cb.c__cvt_adj_ens_pdaf,
                                      pdaf_cb.c__cvt_pdaf,
                                      pdaf_cb.c__cvt_adj_pdaf,
                                      pdaf_cb.c__obs_op_lin_pdaf,
                                      pdaf_cb.c__obs_op_adj_pdaf,
                                      pdaf_cb.c__init_n_domains_p_pdaf,
                                      pdaf_cb.c__init_dim_l_pdaf,
                                      pdaf_cb.c__init_dim_obs_l_pdaf,
                                      pdaf_cb.c__prepoststep_pdaf,
                                      pdaf_cb.c__next_observation_pdaf,
                                      &outflag)

    return outflag

def assim_offline_3dvar_all(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation for all types of 3DVar DA for a single DA step
    using diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.put_state_3dvar_nondiagR`,
    :func:`pyPDAF.PDAF3.put_state_en3dvar_lestkf_nondiagR`,
    :func:`pyPDAF.PDAF3.put_state_en3dvar_estkf_nondiagR`,
    :func:`pyPDAF.PDAF3.put_state_hyb3dvar_lestkf_nondiagR`,
    :func:`pyPDAF.PDAF3.put_state_hyb3dvar_estkf_nondiagR`
    for non-diagonal observation error covariance matrix.

    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable transformation.
    This is a deterministic filtering scheme
    so no ensemble and parallelisation is needed.
    This function should be called at each model time step.

    For parametrised 3DVar, user-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__obs_op_adj_pdaf
            4. py__cvt_adj_pdaf
            5. core DA algorithm
        5. py__cvt_pdaf
        6. py__prepoststep_state_pdaf

    For 3DEnVar, user-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__obs_op_adj_pdaf
            4. py__cvt_adj_ens_pdaf
            5. core DA algorithm
        5. py__cvt_ens_pdaf
        6. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. core DA algorithm
        7. py__prepoststep_state_pdaf


    For hybrid 3DVar, user-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. py__cvt_adj_ens_pdaf
            7. core DA algorithm
        5. py__cvt_pdaf
        6. py__cvt_ens_pdaf
        7. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. core DA algorithm
        8. py__prepoststep_state_pdaf


    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_3dvar_all(pdaf_cb.c__init_dim_obs_pdaf,
                                         pdaf_cb.c__obs_op_pdaf,
                                         pdaf_cb.c__cvt_ens_pdaf,
                                         pdaf_cb.c__cvt_adj_ens_pdaf,
                                         pdaf_cb.c__cvt_pdaf,
                                         pdaf_cb.c__cvt_adj_pdaf,
                                         pdaf_cb.c__obs_op_lin_pdaf,
                                         pdaf_cb.c__obs_op_adj_pdaf,
                                         pdaf_cb.c__init_n_domains_p_pdaf,
                                         pdaf_cb.c__init_dim_l_pdaf,
                                         pdaf_cb.c__init_dim_obs_l_pdaf,
                                         pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag

def assimilate_local_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Online assimilation of domain local filters for a single DA step
    using non-diagnoal observation error covariance matrix.

    Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__init_obs_l_pdaf
            4. py__prodRinvA_l_pdaf
            5. core DA algorithm
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_local_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__init_n_domains_p_pdaf,
                                           pdaf_cb.c__init_dim_l_pdaf,
                                           pdaf_cb.c__init_dim_obs_l_pdaf,
                                           pdaf_cb.c__prodrinva_l_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_global_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    r"""Online assimilation of global filters except for 3DVar and stochastic EnKF
    for a single DA step using non-diagnoal observation
    error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate`
    for diagonal observation error covariance matrix.

    Here, this function call is used for global, E(S)TKF [1]_,
    SEIK [1]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__prodRinvA_pdaf
        7. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product of inverse of R with matrix A
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_global_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                            pdaf_cb.c__distribute_state_pdaf,
                                            pdaf_cb.c__init_dim_obs_pdaf,
                                            pdaf_cb.c__obs_op_pdaf,
                                            pdaf_cb.c__prodrinva_pdaf,
                                            pdaf_cb.c__prepoststep_pdaf,
                                            pdaf_cb.c__next_observation_pdaf,
                                            &outflag)

    return outflag


def assimilate_lnetf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__likelihood_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Online assimilation of LNETF for a single DA step using
    non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate` for
    using diagnoal observation error covariance matrix.
    The non-linear filter is proposed in [1]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__likelihood_l_pdaf
            4. core DA algorithm
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_lnetf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__init_n_domains_p_pdaf,
                                           pdaf_cb.c__init_dim_l_pdaf,
                                           pdaf_cb.c__init_dim_obs_l_pdaf,
                                           pdaf_cb.c__likelihood_l_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_lknetf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf,
    py__prodrinva_hyb_l_pdaf, py__likelihood_l_pdaf,
    py__likelihood_hyb_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Online assimilation of LKNETF for a single DA step using
    non-diagonal observation error covariance matrix.


    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    LKNETF [1]_ for a single DA step using non-diagnoal
    observation error covariance matrix.
    See :func:`pyPDAF.PDAF3.assimilate`
    for using diagnoal observation error covariance matrix.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_n_domains_p_pdaf
        4. py__init_dim_obs_pdaf
        5. py__obs_op_pdaf (for each ensemble member)
        6. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__prodRinvA_pdaf
            4. py__likelihood_l_pdaf
            5. core DA algorithm
            6. py__obs_op_pdaf
               (only called with `HKN` and `HNK` options
               called for each ensemble member)
            7. py__likelihood_hyb_l_pdaf
            8. py__prodRinvA_hyb_l_pdaf
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A
    py__prodrinva_hyb_l_pdaf : Callable
        Product R^-1 A on local analysis domain with hybrid weight
    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization
    py__likelihood_hyb_l_pdaf : Callable
        Compute likelihood and apply localization with tempering
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_lknetf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                            pdaf_cb.c__distribute_state_pdaf,
                                            pdaf_cb.c__init_dim_obs_pdaf,
                                            pdaf_cb.c__obs_op_pdaf,
                                            pdaf_cb.c__init_n_domains_p_pdaf,
                                            pdaf_cb.c__init_dim_l_pdaf,
                                            pdaf_cb.c__init_dim_obs_l_pdaf,
                                            pdaf_cb.c__prodrinva_l_pdaf,
                                            pdaf_cb.c__prodrinva_hyb_l_pdaf,
                                            pdaf_cb.c__likelihood_l_pdaf,
                                            pdaf_cb.c__likelihood_hyb_l_pdaf,
                                            pdaf_cb.c__prepoststep_pdaf,
                                            pdaf_cb.c__next_observation_pdaf,
                                            &outflag)

    return outflag


def assimilate_enkf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    r"""Online assimilation of global or Covariance localised stochastic EnKF
    for a single DA step using non-diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate`
    for diagonal observation error covariance matrix.

    This stochastic EnKF is implemented based on [1]_

    This is the only scheme for covariance localisation with non-diagonal
    observation error covariance matrix in PDAF.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__localize_pdaf
        6. py__add_obs_err_pdaf
        7. py__init_obscovar_pdaf
        8. py__obs_op_pdaf (repeated to reduce storage)
        9. core DA algorith
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix
    py__init_obs_covar_pdaf : Callable
        Initialize mean observation error variance
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_enkf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__distribute_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__add_obs_err_pdaf,
                                          pdaf_cb.c__init_obs_covar_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__next_observation_pdaf,
                                          &outflag)

    return outflag


def assimilate_lenkf_nondiagr(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__localize_covar_pdaf, py__add_obs_err_pdaf, py__init_obs_covar_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    r"""Covariance localised stochastic EnKF
    for a single DA step using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate` or :func:`pyPDAF.PDAF3.assim_offline`
    for diagnoal observation error covariance matrix.

    This stochastic EnKF is implemented based on [1]_

    This is the only scheme for covariance localisation with non-diagonal
    observation error covariance matrix in PDAF.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__localize_pdaf
        6. py__add_obs_err_pdaf
        7. py__init_obscovar_pdaf
        8. py__obs_op_pdaf (repeated to reduce storage)
        9. core DA algorith
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
   py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__localize_covar_pdaf : Callable
        Apply covariance localization
    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix
    py__init_obs_covar_pdaf : Callable
        Initialize mean observation error variance
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_lenkf_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__localize_covar_pdaf,
                                           pdaf_cb.c__add_obs_err_pdaf,
                                           pdaf_cb.c__init_obs_covar_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_nonlin_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__likelihood_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    r"""Online assimilation of global nonlinear filters for a single DA step
    using non-diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate`
    for diagonal observation error covariance matrix.

    Here, this function call is used for global NETF [1]_,
    and particle filter [2]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR` and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for ensemble mean)
        5. py__obs_op_pdaf (for each ensemble member)
        6. py__likelihood_pdaf
        7. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.
    .. [2] Van Leeuwen, P. J., Künsch, H. R., Nerger, L.,
           Potthast, R., & Reich, S. (2019).
           Particle filters for high‐dimensional geoscience applications:
           A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__likelihood_pdaf : Callable
        Compute likelihood
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_nonlin_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                            pdaf_cb.c__distribute_state_pdaf,
                                            pdaf_cb.c__init_dim_obs_pdaf,
                                            pdaf_cb.c__obs_op_pdaf,
                                            pdaf_cb.c__likelihood_pdaf,
                                            pdaf_cb.c__prepoststep_pdaf,
                                            pdaf_cb.c__next_observation_pdaf,
                                            &outflag)

    return outflag








def assimilate_3dvar_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    r"""3DVar DA for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate_3dvar_all`
    for diagonal observation error covariance matrix.

    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable transformation.
    This is a deterministic filtering scheme
    so no ensemble and parallelisation is needed.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. core DA algorithm
        6. py__cvt_pdaf
        7. py__prepoststep_state_pdaf
        8. py__distribute_state_pdaf
        9. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_3dvar_nondiagr(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__prodrinva_pdaf,
                                           pdaf_cb.c__cvt_pdaf,
                                           pdaf_cb.c__cvt_adj_pdaf,
                                           pdaf_cb.c__obs_op_lin_pdaf,
                                           pdaf_cb.c__obs_op_adj_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_en3dvar_estkf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    r"""3DEnVar for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate_3dvar_all`
    for diagonal observation error covariance matirx.

    Here, the background error covariance matrix is
    estimated by an ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An ESTKF is used along with 3DEnVar to generate ensemble perturbations.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core 3DEnVar algorithm
        6. py__cvt_ens_pdaf
        7. ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf (for ensemble mean)
            3. py__obs_op_pdaf (for each ensemble member)
            4. py__prodRinvA_pdaf
            5. core ESTKF algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_en3dvar_estkf_nondiagr(
                                                   pdaf_cb.c__collect_state_pdaf,
                                                   pdaf_cb.c__distribute_state_pdaf,
                                                   pdaf_cb.c__init_dim_obs_pdaf,
                                                   pdaf_cb.c__obs_op_pdaf,
                                                   pdaf_cb.c__prodrinva_pdaf,
                                                   pdaf_cb.c__cvt_ens_pdaf,
                                                   pdaf_cb.c__cvt_adj_ens_pdaf,
                                                   pdaf_cb.c__obs_op_lin_pdaf,
                                                   pdaf_cb.c__obs_op_adj_pdaf,
                                                   pdaf_cb.c__prepoststep_pdaf,
                                                   pdaf_cb.c__next_observation_pdaf,
                                                   &outflag)

    return outflag


def assimilate_en3dvar_lestkf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF using non-diagonal observation
    error covariance matrix.

    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF
    using non-diagnoal observation error covariance matrix.
    The background error covariance matrix is estimated by ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An LESTKF is used to generate ensemble perturbations.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core DA algorithm
        6. py__cvt_ens_pdaf
        7. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__prodRinvA_l_pdaf
                4. core DA algorithm
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf


    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_en3dvar_lestkf_nondiagr(
                                                    pdaf_cb.c__collect_state_pdaf,
                                                    pdaf_cb.c__distribute_state_pdaf,
                                                    pdaf_cb.c__init_dim_obs_pdaf,
                                                    pdaf_cb.c__obs_op_pdaf,
                                                    pdaf_cb.c__prodrinva_pdaf,
                                                    pdaf_cb.c__cvt_ens_pdaf,
                                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                                    pdaf_cb.c__obs_op_lin_pdaf,
                                                    pdaf_cb.c__obs_op_adj_pdaf,
                                                    pdaf_cb.c__prodrinva_l_pdaf,
                                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                                    pdaf_cb.c__init_dim_l_pdaf,
                                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                                    pdaf_cb.c__prepoststep_pdaf,
                                                    pdaf_cb.c__next_observation_pdaf,
                                                    &outflag)

    return outflag


def assimilate_hyb3dvar_estkf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    r"""Hybrid 3DEnVar for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate_3dvar_all`
    for diagonal observation error covariance matrix.

    Here the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    ESTKF in this implementation.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. the iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core 3DEnVar algorithm
        6. py__cvt_pdaf
        7. py__cvt_ens_pdaf
        8. Perform ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf
               (for ensemble mean)
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__prodRinvA_pdaf
            5. core ESTKF algorithm
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_hyb3dvar_estkf_nondiagr(
                                                    pdaf_cb.c__collect_state_pdaf,
                                                    pdaf_cb.c__distribute_state_pdaf,
                                                    pdaf_cb.c__init_dim_obs_pdaf,
                                                    pdaf_cb.c__obs_op_pdaf,
                                                    pdaf_cb.c__prodrinva_pdaf,
                                                    pdaf_cb.c__cvt_ens_pdaf,
                                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                                    pdaf_cb.c__cvt_pdaf,
                                                    pdaf_cb.c__cvt_adj_pdaf,
                                                    pdaf_cb.c__obs_op_lin_pdaf,
                                                    pdaf_cb.c__obs_op_adj_pdaf,
                                                    pdaf_cb.c__prepoststep_pdaf,
                                                    pdaf_cb.c__next_observation_pdaf,
                                                    &outflag)

    return outflag


def assimilate_hyb3dvar_lestkf_nondiagr(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    r"""Hybrid 3DEnVar for a single DA step
    using non-diagonal observation error covariance matrix.

    Here, the background error covariance is
    hybridised by a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    LESTKF in this implementation.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core DA algorithm
        6. py__cvt_pdaf
        7. py__cvt_ens_pdaf
        8. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__prodRinvA_l_pdaf
                4. core DA algorithm
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_hyb3dvar_lestkf_nondiagr(
                                                     pdaf_cb.c__collect_state_pdaf,
                                                     pdaf_cb.c__distribute_state_pdaf,
                                                     pdaf_cb.c__init_dim_obs_pdaf,
                                                     pdaf_cb.c__obs_op_pdaf,
                                                     pdaf_cb.c__prodrinva_pdaf,
                                                     pdaf_cb.c__cvt_ens_pdaf,
                                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                                     pdaf_cb.c__cvt_pdaf,
                                                     pdaf_cb.c__cvt_adj_pdaf,
                                                     pdaf_cb.c__obs_op_lin_pdaf,
                                                     pdaf_cb.c__obs_op_adj_pdaf,
                                                     pdaf_cb.c__prodrinva_l_pdaf,
                                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                                     pdaf_cb.c__init_dim_l_pdaf,
                                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                                     pdaf_cb.c__prepoststep_pdaf,
                                                     pdaf_cb.c__next_observation_pdaf,
                                                     &outflag)

    return outflag







def assim_offline_local_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation of domain local filters for a single DA step
    using non-diagnoal observation error covariance matrix.

    Here, this function call is used for LE(S)TKF [1]_ and LSEIK [1]_
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_n_domains_p_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__init_obs_l_pdaf
            4. py__prodRinvA_l_pdaf
            5. core DA algorithm
        6. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_local_nondiagr(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__init_n_domains_p_pdaf,
                                              pdaf_cb.c__init_dim_l_pdaf,
                                              pdaf_cb.c__init_dim_obs_l_pdaf,
                                              pdaf_cb.c__prodrinva_l_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_global_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation of global filters except for 3DVar and stochastic EnKF
    for a single DA step using non-diagnoal observation
    error covariance matrix.

    See :func:`pyPDAF.PDAF3.assimilate`
    for diagonal observation error covariance matrix.

    Here, this function call is used for global, E(S)TKF [1]_,
    SEIK [1]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf (for ensemble mean)
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__prodRinvA_pdaf
        6. core DA algorithm
        7. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Nerger, L., Janjić, T., Schröter, J., Hiller, W. (2012).
           A unification of ensemble square root Kalman filters.
           Monthly Weather Review, 140, 2335-2345.
           doi:10.1175/MWR-D-11-00102.1
    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prodrinva_pdaf : Callable
        Provide product of inverse of R with matrix A

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_global_nondiagr(
                                               pdaf_cb.c__init_dim_obs_pdaf,
                                               pdaf_cb.c__obs_op_pdaf,
                                               pdaf_cb.c__prodrinva_pdaf,
                                               pdaf_cb.c__prepoststep_pdaf,
                                               &outflag)

    return outflag


def assim_offline_lnetf_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__likelihood_l_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation of LNETF for a single DA step using
    non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline` for
    using diagnoal observation error covariance matrix.
    The non-linear filter is proposed in [1]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_n_domains_p_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__likelihood_l_pdaf
            4. core DA algorithm
        6. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_lnetf_nondiagr(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__init_n_domains_p_pdaf,
                                              pdaf_cb.c__init_dim_l_pdaf,
                                              pdaf_cb.c__init_dim_obs_l_pdaf,
                                              pdaf_cb.c__likelihood_l_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_lknetf_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prodrinva_l_pdaf,
    py__prodrinva_hyb_l_pdaf, py__likelihood_l_pdaf, py__likelihood_hyb_l_pdaf,
    py__prepoststep_pdaf):
    r"""Offline assimilation of LKNETF for a single DA step using
    non-diagonal observation error covariance matrix.


    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    LKNETF [1]_ for a single DA step using non-diagnoal
    observation error covariance matrix.
    See :func:`pyPDAF.PDAF3.assimilate`
    for using diagnoal observation error covariance matrix.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_n_domains_p_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. loop over each local domain:
            1. py__init_dim_l_pdaf
            2. py__init_dim_obs_l_pdaf
            3. py__prodRinvA_pdaf
            4. py__likelihood_l_pdaf
            5. core DA algorithm
            6. py__obs_op_pdaf
               (only called with `HKN` and `HNK` options
               called for each ensemble member)
            7. py__likelihood_hyb_l_pdaf
            8. py__prodRinvA_hyb_l_pdaf
        6. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prodrinva_l_pdaf : Callable
        Provide product of inverse of R with matrix A

    py__prodrinva_hyb_l_pdaf : Callable
        Product R^-1 A on local analysis domain with hybrid weight

    py__likelihood_l_pdaf : Callable
        Compute likelihood and apply localization

    py__likelihood_hyb_l_pdaf : Callable
        Compute likelihood and apply localization with tempering

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.prodrinva_hyb_l_pdaf = <void*>py__prodrinva_hyb_l_pdaf
    pdaf_cb.likelihood_l_pdaf = <void*>py__likelihood_l_pdaf
    pdaf_cb.likelihood_hyb_l_pdaf = <void*>py__likelihood_hyb_l_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_lknetf_nondiagr(
                                               pdaf_cb.c__init_dim_obs_pdaf,
                                               pdaf_cb.c__obs_op_pdaf,
                                               pdaf_cb.c__init_n_domains_p_pdaf,
                                               pdaf_cb.c__init_dim_l_pdaf,
                                               pdaf_cb.c__init_dim_obs_l_pdaf,
                                               pdaf_cb.c__prodrinva_l_pdaf,
                                               pdaf_cb.c__prodrinva_hyb_l_pdaf,
                                               pdaf_cb.c__likelihood_l_pdaf,
                                               pdaf_cb.c__likelihood_hyb_l_pdaf,
                                               pdaf_cb.c__prepoststep_pdaf,
                                               &outflag)

    return outflag


def assim_offline_enkf_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__add_obs_err_pdaf, py__init_obs_covar_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation of global or Covariance localised stochastic EnKF
    for a single DA step using non-diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline`
    for diagonal observation error covariance matrix.

    This stochastic EnKF is implemented based on [1]_

    This is the only scheme for covariance localisation with non-diagonal
    observation error covariance matrix in PDAF.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf (for each ensemble member)
        4. py__localize_pdaf
        5. py__add_obs_err_pdaf
        6. py__init_obscovar_pdaf
        7. py__obs_op_pdaf (repeated to reduce storage)
        8. core DA algorithm
        9. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

    py__init_obs_covar_pdaf : Callable
        Initialize mean observation error variance

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_enkf_nondiagr(pdaf_cb.c__init_dim_obs_pdaf,
                                             pdaf_cb.c__obs_op_pdaf,
                                             pdaf_cb.c__add_obs_err_pdaf,
                                             pdaf_cb.c__init_obs_covar_pdaf,
                                             pdaf_cb.c__prepoststep_pdaf,
                                             &outflag)

    return outflag


def assim_offline_lenkf_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__localize_covar_pdaf, py__add_obs_err_pdaf,
    py__init_obs_covar_pdaf, py__prepoststep_pdaf):
    r"""Online assimilation of covariance localised stochastic EnKF
    for a single DA step using non-diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline`
    for diagonal observation error covariance matrix.

    This stochastic EnKF is implemented based on [1]_

    This is the only scheme for covariance localisation with non-diagonal
    observation error covariance matrix in PDAF.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf (for each ensemble member)
        4. py__localize_pdaf
        5. py__add_obs_err_pdaf
        6. py__init_obscovar_pdaf
        7. py__obs_op_pdaf (repeated to reduce storage)
        8. core DA algorithm
        9. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__localize_covar_pdaf : Callable
        Apply covariance localization

    py__add_obs_err_pdaf : Callable
        Add observation error covariance matrix

    py__init_obs_covar_pdaf : Callable
        Initialize mean observation error variance

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.add_obs_err_pdaf = <void*>py__add_obs_err_pdaf
    pdaf_cb.init_obs_covar_pdaf = <void*>py__init_obs_covar_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_lenkf_nondiagr(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__localize_covar_pdaf,
                                              pdaf_cb.c__add_obs_err_pdaf,
                                              pdaf_cb.c__init_obs_covar_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_nonlin_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__likelihood_pdaf, py__prepoststep_pdaf):
    r"""Offline assimilation of global nonlinear filters for a single DA step
    using non-diagonal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline`
    for diagonal observation error covariance matrix.

    Here, this function call is used for global NETF [1]_,
    and particle filter [2]_.
    The filter type is set in :func:`pyPDAF.PDAF.init`.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.omi_put_state_global_nondiagR` and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf (for ensemble mean)
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__likelihood_pdaf
        6. core DA algorithm
        7. py__prepoststep_state_pdaf

    References
    ----------
    .. [1] Tödter, J., and B. Ahrens, 2015:
           A second-order exact ensemble square root filter
           for nonlinear data assimilation. Mon. Wea. Rev.,
           143, 1347–1367, doi:10.1175/MWR-D-14-00108.1.
    .. [2] Van Leeuwen, P. J., Künsch, H. R., Nerger, L.,
           Potthast, R., & Reich, S. (2019).
           Particle filters for high‐dimensional geoscience applications:
           A review. Quarterly Journal of the Royal Meteorological Society, 145(723), 2335-2365.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__likelihood_pdaf : Callable
        Compute likelihood

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.likelihood_pdaf = <void*>py__likelihood_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_nonlin_nondiagr(
                                               pdaf_cb.c__init_dim_obs_pdaf,
                                               pdaf_cb.c__obs_op_pdaf,
                                               pdaf_cb.c__likelihood_pdaf,
                                               pdaf_cb.c__prepoststep_pdaf,
                                               &outflag)

    return outflag





def assimilate_3dvar(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DVar DA for a single step without OMI.
    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable
    transformation. This is a deterministic filtering
    scheme so no ensemble and
    parallelisation is needed.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_3dvar`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. core DA algorithm
        7. py__cvt_pdaf
        8. py__prepoststep_state_pdaf
        9. py__distribute_state_pdaf
        10. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by :func:`pyPDAF.PDAF.omi_assimilate_3dvar`
       and :func:`pyPDAF.PDAF.omi_assimilate_3dvar_nondiagR`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_3dvar(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__cvt_pdaf,
                                  pdaf_cb.c__cvt_adj_pdaf,
                                  pdaf_cb.c__obs_op_lin_pdaf,
                                  pdaf_cb.c__obs_op_adj_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_en3dvar(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_en3dvar(pdaf_cb.c__collect_state_pdaf,
                                    pdaf_cb.c__distribute_state_pdaf,
                                    pdaf_cb.c__init_dim_obs_pdaf,
                                    pdaf_cb.c__obs_op_pdaf,
                                    pdaf_cb.c__cvt_ens_pdaf,
                                    pdaf_cb.c__cvt_adj_ens_pdaf,
                                    pdaf_cb.c__obs_op_lin_pdaf,
                                    pdaf_cb.c__obs_op_adj_pdaf,
                                    pdaf_cb.c__init_n_domains_p_pdaf,
                                    pdaf_cb.c__init_dim_l_pdaf,
                                    pdaf_cb.c__init_dim_obs_l_pdaf,
                                    pdaf_cb.c__prepoststep_pdaf,
                                    pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_en3dvar_estkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step.
    The background error covariance matrix is estimated
    by an ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An ESTKF is used along with 3DEnVar
    to generate ensemble perturbations.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_en3dvar_estkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core 3DEnVar algorithm
        7. py__cvt_ens_pdaf
        8. ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf (for ensemble mean)
            3. py__init_obs_pdaf
            4. py__obs_op_pdaf (for each ensemble member)
            5. py__init_obsvar_pdaf
               (only relevant for adaptive
               forgetting factor schemes)
            6. py__prodRinvA_pdaf
            7. core ESTKF algorithm
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf`
       and :func:`pyPDAF.PDAF.omi_assimilate_en3dvar_estkf_nondiagR`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_en3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                          pdaf_cb.c__distribute_state_pdaf,
                                          pdaf_cb.c__init_dim_obs_pdaf,
                                          pdaf_cb.c__obs_op_pdaf,
                                          pdaf_cb.c__cvt_ens_pdaf,
                                          pdaf_cb.c__cvt_adj_ens_pdaf,
                                          pdaf_cb.c__obs_op_lin_pdaf,
                                          pdaf_cb.c__obs_op_adj_pdaf,
                                          pdaf_cb.c__prepoststep_pdaf,
                                          pdaf_cb.c__next_observation_pdaf,
                                          &outflag)

    return outflag


def assimilate_en3dvar_lestkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF.
    The background error covariance matrix is
    estimated by ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An LESTKF is used to generate ensemble perturbations.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_en3dvar_lestkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core DA algorithm
        7. py__cvt_ens_pdaf
        8. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__init_obs_pdaf
               (if global adaptive forgetting factor is used
               `type_forget=1` in :func:`pyPDAF.PDAF.init`)
            5. py__init_obsvar_pdaf
               (if global adaptive forgetting factor is used)
            6. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__g2l_state_pdaf
                4. py__g2l_obs_pdaf
                   (localise mean ensemble in observation space)
                5. py__init_obs_l_pdaf
                6. py__g2l_obs_pdaf
                   (localise each ensemble member
                   in observation space)
                7. py__init_obsvar_l_pdaf
                   (only called if local
                   adaptive forgetting factor
                   `type_forget=2` is used)
                8. py__prodRinvA_l_pdaf
                9. core DA algorithm
                10. py__l2g_state_pdaf
        9. py__prepoststep_state_pdaf
        10. py__distribute_state_pdaf
        11. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf`
       and
       :func:`pyPDAF.PDAF.localomi_assimilate_en3dvar_lestkf_nondiagR`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_en3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__cvt_ens_pdaf,
                                           pdaf_cb.c__cvt_adj_ens_pdaf,
                                           pdaf_cb.c__obs_op_lin_pdaf,
                                           pdaf_cb.c__obs_op_adj_pdaf,
                                           pdaf_cb.c__init_n_domains_p_pdaf,
                                           pdaf_cb.c__init_dim_l_pdaf,
                                           pdaf_cb.c__init_dim_obs_l_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_hyb3dvar(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_hyb3dvar(pdaf_cb.c__collect_state_pdaf,
                                     pdaf_cb.c__distribute_state_pdaf,
                                     pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__cvt_ens_pdaf,
                                     pdaf_cb.c__cvt_adj_ens_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                     pdaf_cb.c__init_dim_l_pdaf,
                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf,
                                     pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_hyb3dvar_estkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    ESTKF in this implementation.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_hyb3dvar_estkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. the iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core 3DEnVar algorithm
        7. py__cvt_pdaf
        8. py__cvt_ens_pdaf
        9. Perform ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf
               (for ensemble mean)
            3. py__init_obs_pdaf
            4. py__obs_op_pdaf
               (for each ensemble member)
            5. py__init_obsvar_pdaf
               (only relevant for adaptive
               forgetting factor schemes)
            6. py__prodRinvA_pdaf
            7. core ESTKF algorithm
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf` and
       :func:`pyPDAF.PDAF.omi_assimilate_hyb3dvar_estkf_nondiagR`.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_hyb3dvar_estkf(pdaf_cb.c__collect_state_pdaf,
                                           pdaf_cb.c__distribute_state_pdaf,
                                           pdaf_cb.c__init_dim_obs_pdaf,
                                           pdaf_cb.c__obs_op_pdaf,
                                           pdaf_cb.c__cvt_ens_pdaf,
                                           pdaf_cb.c__cvt_adj_ens_pdaf,
                                           pdaf_cb.c__cvt_pdaf,
                                           pdaf_cb.c__cvt_adj_pdaf,
                                           pdaf_cb.c__obs_op_lin_pdaf,
                                           pdaf_cb.c__obs_op_adj_pdaf,
                                           pdaf_cb.c__prepoststep_pdaf,
                                           pdaf_cb.c__next_observation_pdaf,
                                           &outflag)

    return outflag


def assimilate_hyb3dvar_lestkf(py__collect_state_pdaf,
    py__distribute_state_pdaf, py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf, int  outflag):
    """It is recommended to use :func:`pyPDAF.PDAF3.assimilate_3dvar_all`

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Hybrid 3DEnVar for a single DA step where
    the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    LESTKF in this implementation.
    This function should be called at each model time step.

    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_hyb3dvar_lestkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__init_obs_pdaf
        6. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core DA algorithm
        7. py__cvt_pdaf
        8. py__cvt_ens_pdaf
        9. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__init_obs_pdaf
               (if global adaptive forgetting factor
               `type_forget=1` in :func:`pyPDAF.PDAF.init`)
            5. py__init_obsvar_pdaf
               (if global adaptive forgetting factor is used)
            6. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__g2l_state_pdaf
                4. py__g2l_obs_pdaf
                   (localise mean ensemble in observation space)
                5. py__init_obs_l_pdaf
                6. py__g2l_obs_pdaf
                   (localise each ensemble member
                   in observation space)
                7. py__init_obsvar_l_pdaf
                   (only called if local adaptive forgetting
                   factor `type_forget=2` is used)
                8. py__prodRinvA_l_pdaf
                9. core DA algorithm
                10. py__l2g_state_pdaf
        10. py__prepoststep_state_pdaf
        11. py__distribute_state_pdaf
        12. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf`
       and
       :func:`pyPDAF.PDAF.localomi_assimilate_hyb3dvar_lestkf_nondiagR`

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_hyb3dvar_lestkf(pdaf_cb.c__collect_state_pdaf,
                                            pdaf_cb.c__distribute_state_pdaf,
                                            pdaf_cb.c__init_dim_obs_pdaf,
                                            pdaf_cb.c__obs_op_pdaf,
                                            pdaf_cb.c__cvt_ens_pdaf,
                                            pdaf_cb.c__cvt_adj_ens_pdaf,
                                            pdaf_cb.c__cvt_pdaf,
                                            pdaf_cb.c__cvt_adj_pdaf,
                                            pdaf_cb.c__obs_op_lin_pdaf,
                                            pdaf_cb.c__obs_op_adj_pdaf,
                                            pdaf_cb.c__init_n_domains_p_pdaf,
                                            pdaf_cb.c__init_dim_l_pdaf,
                                            pdaf_cb.c__init_dim_obs_l_pdaf,
                                            pdaf_cb.c__prepoststep_pdaf,
                                            pdaf_cb.c__next_observation_pdaf,
                                            &outflag)

    return outflag


def generate_obs(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__get_obs_f_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf, int  outflag):
    """Generation of synthetic observations based on
    given error statistics and observation operator.

    The generated synthetic observations are based on
    each member of model forecast.
    Therefore, an ensemble of observations can be obtained.
    In a typical experiment,
    one may only need one ensemble member.
    The implementation strategy is similar to
    an assimilation step. This means that,
    one can reuse many user-supplied functions for
    assimilation and observation generation.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf
        5. py__get_obs_f_pdaf
        6. py__prepoststep_state_pdaf
        7. py__distribute_state_pdaf
        8. py__next_observation_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__get_obs_f_pdaf : Callable
        Initialize observation vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_generate_obs(pdaf_cb.c__collect_state_pdaf,
                              pdaf_cb.c__distribute_state_pdaf,
                              pdaf_cb.c__init_dim_obs_pdaf,
                              pdaf_cb.c__obs_op_pdaf,
                              pdaf_cb.c__get_obs_f_pdaf,
                              pdaf_cb.c__prepoststep_pdaf,
                              pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def generate_obs_offline(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__get_obs_f_pdaf, py__prepoststep_pdaf):
    """Generation of synthetic observations based on
    given error statistics and observation operator in offline setup.

    The generated synthetic observations are based on
    each member of model forecast.
    Therefore, an ensemble of observations can be obtained.
    In a typical experiment,
    one may only need one ensemble member.
    The implementation strategy is similar to
    an assimilation step. This means that,
    one can reuse many user-supplied functions for
    assimilation and observation generation.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. py__get_obs_f_pdaf
        5. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__get_obs_f_pdaf : Callable
        Initialize observation vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.get_obs_f_pdaf = <void*>py__get_obs_f_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_generate_obs_offline(pdaf_cb.c__init_dim_obs_pdaf,
                                      pdaf_cb.c__obs_op_pdaf,
                                      pdaf_cb.c__get_obs_f_pdaf,
                                      pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_3dvar(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_pdaf, py__cvt_adj_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_3dvar(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__cvt_pdaf,
                                     pdaf_cb.c__cvt_adj_pdaf,
                                     pdaf_cb.c__obs_op_lin_pdaf,
                                     pdaf_cb.c__obs_op_adj_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_en3dvar(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_en3dvar(pdaf_cb.c__init_dim_obs_pdaf,
                                       pdaf_cb.c__obs_op_pdaf,
                                       pdaf_cb.c__cvt_ens_pdaf,
                                       pdaf_cb.c__cvt_adj_ens_pdaf,
                                       pdaf_cb.c__obs_op_lin_pdaf,
                                       pdaf_cb.c__obs_op_adj_pdaf,
                                       pdaf_cb.c__init_n_domains_p_pdaf,
                                       pdaf_cb.c__init_dim_l_pdaf,
                                       pdaf_cb.c__init_dim_obs_l_pdaf,
                                       pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_en3dvar_estkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_en3dvar_estkf(pdaf_cb.c__init_dim_obs_pdaf,
                                             pdaf_cb.c__obs_op_pdaf,
                                             pdaf_cb.c__cvt_ens_pdaf,
                                             pdaf_cb.c__cvt_adj_ens_pdaf,
                                             pdaf_cb.c__obs_op_lin_pdaf,
                                             pdaf_cb.c__obs_op_adj_pdaf,
                                             pdaf_cb.c__prepoststep_pdaf,
                                             &outflag)

    return outflag


def assim_offline_en3dvar_lestkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf,
    py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_en3dvar_lestkf(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__cvt_ens_pdaf,
                                              pdaf_cb.c__cvt_adj_ens_pdaf,
                                              pdaf_cb.c__obs_op_lin_pdaf,
                                              pdaf_cb.c__obs_op_adj_pdaf,
                                              pdaf_cb.c__init_n_domains_p_pdaf,
                                              pdaf_cb.c__init_dim_l_pdaf,
                                              pdaf_cb.c__init_dim_obs_l_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_hyb3dvar(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_hyb3dvar(pdaf_cb.c__init_dim_obs_pdaf,
                                        pdaf_cb.c__obs_op_pdaf,
                                        pdaf_cb.c__cvt_ens_pdaf,
                                        pdaf_cb.c__cvt_adj_ens_pdaf,
                                        pdaf_cb.c__cvt_pdaf,
                                        pdaf_cb.c__cvt_adj_pdaf,
                                        pdaf_cb.c__obs_op_lin_pdaf,
                                        pdaf_cb.c__obs_op_adj_pdaf,
                                        pdaf_cb.c__init_n_domains_p_pdaf,
                                        pdaf_cb.c__init_dim_l_pdaf,
                                        pdaf_cb.c__init_dim_obs_l_pdaf,
                                        pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_hyb3dvar_estkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_hyb3dvar_estkf(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__cvt_ens_pdaf,
                                              pdaf_cb.c__cvt_adj_ens_pdaf,
                                              pdaf_cb.c__cvt_pdaf,
                                              pdaf_cb.c__cvt_adj_pdaf,
                                              pdaf_cb.c__obs_op_lin_pdaf,
                                              pdaf_cb.c__obs_op_adj_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_hyb3dvar_lestkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__cvt_ens_pdaf, py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_hyb3dvar_lestkf(
                                               pdaf_cb.c__init_dim_obs_pdaf,
                                               pdaf_cb.c__obs_op_pdaf,
                                               pdaf_cb.c__cvt_ens_pdaf,
                                               pdaf_cb.c__cvt_adj_ens_pdaf,
                                               pdaf_cb.c__cvt_pdaf,
                                               pdaf_cb.c__cvt_adj_pdaf,
                                               pdaf_cb.c__obs_op_lin_pdaf,
                                               pdaf_cb.c__obs_op_adj_pdaf,
                                               pdaf_cb.c__init_n_domains_p_pdaf,
                                               pdaf_cb.c__init_dim_l_pdaf,
                                               pdaf_cb.c__init_dim_obs_l_pdaf,
                                               pdaf_cb.c__prepoststep_pdaf,
                                               &outflag)

    return outflag


def assim_offline_3dvar_nondiagr(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prodrinva_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    r"""Offline 3DVar DA for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline_3dvar_all`
    for diagonal observation error covariance matrix.

    When 3DVar is used, the background error covariance matrix
    has to be modelled for cotrol variable transformation.
    This is a deterministic filtering scheme
    so no ensemble and parallelisation is needed.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. Iterative optimisation:
            1. py__cvt_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_pdaf
            6. core DA algorithm
        5. py__cvt_pdaf
        6. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_3dvar_nondiagr(pdaf_cb.c__init_dim_obs_pdaf,
                                              pdaf_cb.c__obs_op_pdaf,
                                              pdaf_cb.c__prodrinva_pdaf,
                                              pdaf_cb.c__cvt_pdaf,
                                              pdaf_cb.c__cvt_adj_pdaf,
                                              pdaf_cb.c__obs_op_lin_pdaf,
                                              pdaf_cb.c__obs_op_adj_pdaf,
                                              pdaf_cb.c__prepoststep_pdaf,
                                              &outflag)

    return outflag


def assim_offline_en3dvar_estkf_nondiagr(py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__prepoststep_pdaf):
    r"""Offline 3DEnVar for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline_3dvar_all`
    for diagonal observation error covariance matirx.

    Here, the background error covariance matrix is
    estimated by an ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An ESTKF is used along with 3DEnVar to generate ensemble perturbations.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core 3DEnVar algorithm
        5. py__cvt_ens_pdaf
        6. ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf (for ensemble mean)
            3. py__obs_op_pdaf (for each ensemble member)
            4. py__prodRinvA_pdaf
            5. core ESTKF algorithm
        7. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_en3dvar_estkf_nondiagr(
                                                      pdaf_cb.c__init_dim_obs_pdaf,
                                                      pdaf_cb.c__obs_op_pdaf,
                                                      pdaf_cb.c__prodrinva_pdaf,
                                                      pdaf_cb.c__cvt_ens_pdaf,
                                                      pdaf_cb.c__cvt_adj_ens_pdaf,
                                                      pdaf_cb.c__obs_op_lin_pdaf,
                                                      pdaf_cb.c__obs_op_adj_pdaf,
                                                      pdaf_cb.c__prepoststep_pdaf,
                                                      &outflag)

    return outflag


def assim_offline_en3dvar_lestkf_nondiagr(py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__obs_op_lin_pdaf, py__obs_op_adj_pdaf,
    py__prodrinva_l_pdaf, py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    r"""Offline 3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF using non-diagonal observation
    error covariance matrix.

    PDAFlocal-OMI modules require fewer user-supplied
    functions and improved efficiency.

    3DEnVar for a single DA step where the ensemble anomaly
    is generated by LESTKF
    using non-diagnoal observation error covariance matrix.
    The background error covariance matrix is estimated by ensemble.
    The 3DEnVar only calculates the analysis of the ensemble mean.
    An LESTKF is used to generate ensemble perturbations.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. Starting the iterative optimisation:
            1. py__cvt_ens_pdaf
            2. py__obs_op_lin_pdaf
            3. py__prodRinvA_pdaf
            4. py__obs_op_adj_pdaf
            5. py__cvt_adj_ens_pdaf
            6. core DA algorithm
        5. py__cvt_ens_pdaf
        6. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__prodRinvA_l_pdaf
                4. core DA algorithm
        7. py__prepoststep_state_pdaf


    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations
    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains
    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain
    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast
    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_en3dvar_lestkf_nondiagr(
                                                       pdaf_cb.c__init_dim_obs_pdaf,
                                                       pdaf_cb.c__obs_op_pdaf,
                                                       pdaf_cb.c__prodrinva_pdaf,
                                                       pdaf_cb.c__cvt_ens_pdaf,
                                                       pdaf_cb.c__cvt_adj_ens_pdaf,
                                                       pdaf_cb.c__obs_op_lin_pdaf,
                                                       pdaf_cb.c__obs_op_adj_pdaf,
                                                       pdaf_cb.c__prodrinva_l_pdaf,
                                                       pdaf_cb.c__init_n_domains_p_pdaf,
                                                       pdaf_cb.c__init_dim_l_pdaf,
                                                       pdaf_cb.c__init_dim_obs_l_pdaf,
                                                       pdaf_cb.c__prepoststep_pdaf,
                                                       &outflag)

    return outflag


def assim_offline_hyb3dvar_estkf_nondiagr(py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prepoststep_pdaf):
    r"""Offline hybrid 3DEnVar for a single DA step
    using non-diagnoal observation error covariance matrix.

    See :func:`pyPDAF.PDAF3.assim_offline_3dvar_all`
    for diagonal observation error covariance matrix.

    Here the background error covariance is hybridised by
    a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    ESTKF in this implementation.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. the iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core 3DEnVar algorithm
        5. py__cvt_pdaf
        6. py__cvt_ens_pdaf
        7. Perform ESTKF:
            1. py__init_dim_obs_pdaf
            2. py__obs_op_pdaf
               (for ensemble mean)
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. py__prodRinvA_pdaf
            5. core ESTKF algorithm
        8. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector
    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector
    py__obs_op_pdaf : Callable
        Full observation operator
    py__prodrinva_pdaf : Callable
        Provide product R^-1 A
    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector
    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix
    py__obs_op_lin_pdaf : Callable
        Linearized observation operator
    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator
    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine
    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_hyb3dvar_estkf_nondiagr(
                                                       pdaf_cb.c__init_dim_obs_pdaf,
                                                       pdaf_cb.c__obs_op_pdaf,
                                                       pdaf_cb.c__prodrinva_pdaf,
                                                       pdaf_cb.c__cvt_ens_pdaf,
                                                       pdaf_cb.c__cvt_adj_ens_pdaf,
                                                       pdaf_cb.c__cvt_pdaf,
                                                       pdaf_cb.c__cvt_adj_pdaf,
                                                       pdaf_cb.c__obs_op_lin_pdaf,
                                                       pdaf_cb.c__obs_op_adj_pdaf,
                                                       pdaf_cb.c__prepoststep_pdaf,
                                                       &outflag)

    return outflag


def assim_offline_hyb3dvar_lestkf_nondiagr(py__init_dim_obs_pdaf,
    py__obs_op_pdaf, py__prodrinva_pdaf, py__cvt_ens_pdaf,
    py__cvt_adj_ens_pdaf, py__cvt_pdaf, py__cvt_adj_pdaf,
    py__obs_op_lin_pdaf, py__obs_op_adj_pdaf, py__prodrinva_l_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__prepoststep_pdaf):
    r"""Offline hybrid 3DEnVar for a single DA step
    using non-diagonal observation error covariance matrix.

    Here, the background error covariance is
    hybridised by a static background error covariance,
    and a flow-dependent background error covariance
    estimated from ensemble.
    The 3DVar generates an ensemble mean and
    the ensemble perturbation is generated by
    LESTKF in this implementation.
    This function should be called at each model time step.

    User-supplied functions are executed in the following sequence:
        1. py__prepoststep_state_pdaf
        2. py__init_dim_obs_pdaf
        3. py__obs_op_pdaf
        4. The iterative optimisation:
            1. py__cvt_pdaf
            2. py__cvt_ens_pdaf
            3. py__obs_op_lin_pdaf
            4. py__prodRinvA_pdaf
            5. py__obs_op_adj_pdaf
            6. py__cvt_adj_pdaf
            7. py__cvt_adj_ens_pdaf
            8. core DA algorithm
        5. py__cvt_pdaf
        6. py__cvt_ens_pdaf
        7. Perform LESTKF:
            1. py__init_n_domains_p_pdaf
            2. py__init_dim_obs_pdaf
            3. py__obs_op_pdaf
               (for each ensemble member)
            4. loop over each local domain:
                1. py__init_dim_l_pdaf
                2. py__init_dim_obs_l_pdaf
                3. py__prodRinvA_l_pdaf
                4. core DA algorithm
        8. py__prepoststep_state_pdaf

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prodrinva_pdaf : Callable
        Provide product R^-1 A

    py__cvt_ens_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_ens_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__cvt_pdaf : Callable
        Apply control vector transform matrix to control vector

    py__cvt_adj_pdaf : Callable
        Apply adjoint control vector transform matrix

    py__obs_op_lin_pdaf : Callable
        Linearized observation operator

    py__obs_op_adj_pdaf : Callable
        Adjoint observation operator

    py__prodrinva_l_pdaf : Callable
        Provide product R^-1 A and apply localizations

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prodrinva_pdaf = <void*>py__prodrinva_pdaf
    pdaf_cb.cvt_ens_pdaf = <void*>py__cvt_ens_pdaf
    pdaf_cb.cvt_adj_ens_pdaf = <void*>py__cvt_adj_ens_pdaf
    pdaf_cb.cvt_pdaf = <void*>py__cvt_pdaf
    pdaf_cb.cvt_adj_pdaf = <void*>py__cvt_adj_pdaf
    pdaf_cb.obs_op_lin_pdaf = <void*>py__obs_op_lin_pdaf
    pdaf_cb.obs_op_adj_pdaf = <void*>py__obs_op_adj_pdaf
    pdaf_cb.prodrinva_l_pdaf = <void*>py__prodrinva_l_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assim_offline_hyb3dvar_lestkf_nondiagr(
                                                        pdaf_cb.c__init_dim_obs_pdaf,
                                                        pdaf_cb.c__obs_op_pdaf,
                                                        pdaf_cb.c__prodrinva_pdaf,
                                                        pdaf_cb.c__cvt_ens_pdaf,
                                                        pdaf_cb.c__cvt_adj_ens_pdaf,
                                                        pdaf_cb.c__cvt_pdaf,
                                                        pdaf_cb.c__cvt_adj_pdaf,
                                                        pdaf_cb.c__obs_op_lin_pdaf,
                                                        pdaf_cb.c__obs_op_adj_pdaf,
                                                        pdaf_cb.c__prodrinva_l_pdaf,
                                                        pdaf_cb.c__init_n_domains_p_pdaf,
                                                        pdaf_cb.c__init_dim_l_pdaf,
                                                        pdaf_cb.c__init_dim_obs_l_pdaf,
                                                        pdaf_cb.c__prepoststep_pdaf,
                                                        &outflag)

    return outflag

def assimilate_local(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__init_n_domains_p_pdaf,
    py__init_dim_l_pdaf, py__init_dim_obs_l_pdaf, py__g2l_state_pdaf,
    py__l2g_state_pdaf, py__prepoststep_pdaf, py__next_observation_pdaf,
    int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__g2l_state_pdaf : Callable
        Get local state from full state

    py__l2g_state_pdaf : Callable
        Init full state from local state

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    with nogil:
        c__pdaf3_assimilate_local(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__init_n_domains_p_pdaf,
                                  pdaf_cb.c__init_dim_l_pdaf,
                                  pdaf_cb.c__init_dim_obs_l_pdaf,
                                  pdaf_cb.c__g2l_state_pdaf,
                                  pdaf_cb.c__l2g_state_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_global(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__prepoststep_pdaf,
    py__next_observation_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_global(pdaf_cb.c__collect_state_pdaf,
                                   pdaf_cb.c__distribute_state_pdaf,
                                   pdaf_cb.c__init_dim_obs_pdaf,
                                   pdaf_cb.c__obs_op_pdaf,
                                   pdaf_cb.c__prepoststep_pdaf,
                                   pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_lenkf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__localize_covar_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    """It is recommended to use
    :func:`pyPDAF.PDAF.omi_assimilate_lenkf`
    or :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`.

    PDAF-OMI modules require fewer user-supplied functions
    and improved efficiency.

    Stochastic EnKF (ensemble Kalman filter)
    with covariance localisation [1]_
    for a single DA step without OMI.

    This is the only scheme for covariance localisation in PDAF.

    This function should be called at each model time step.
    The function is a combination of
    :func:`pyPDAF.PDAF.put_state_lenkf`
    and :func:`pyPDAF.PDAF.get_state`.

    User-supplied functions are executed in the following sequence:
        1. py__collect_state_pdaf
        2. py__prepoststep_state_pdaf
        3. py__init_dim_obs_pdaf
        4. py__obs_op_pdaf (for each ensemble member)
        5. py__localize_pdaf
        6. py__add_obs_err_pdaf
        7. py__init_obs_pdaf
        8. py__init_obscovar_pdaf
        9. py__obs_op_pdaf (repeated to reduce storage)
        10. core DA algorith
        11. py__prepoststep_state_pdaf
        12. py__distribute_state_pdaf
        13. py__next_observation_pdaf

    .. deprecated:: 1.0.0

       This function is replaced by
       :func:`pyPDAF.PDAF.omi_assimilate_lenkf`
       and :func:`pyPDAF.PDAF.omi_assimilate_lenkf_nondiagR`

    References
    ----------
    .. [1] Houtekamer, P. L., and H. L. Mitchell (1998):
           Data Assimilation Using an Ensemble Kalman Filter
           Technique.
           Mon. Wea. Rev., 126, 796–811,
           doi: 10.1175/1520-0493(1998)126<0796:DAUAEK>2.0.CO;2.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__localize_covar_pdaf : Callable
        Apply localization to HP and HPH^T

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_lenkf(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__localize_covar_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag


def assimilate_ensrf(py__collect_state_pdaf, py__distribute_state_pdaf,
    py__init_dim_obs_pdaf, py__obs_op_pdaf, py__localize_covar_serial_pdaf,
    py__prepoststep_pdaf, py__next_observation_pdaf):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Routine to collect a state vector

    py__distribute_state_pdaf : Callable
        Routine to distribute a state vector

    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__localize_covar_serial_pdaf : Callable
        Apply localization to HP and BXY for single observation

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    py__next_observation_pdaf : Callable
        Provide information on next forecast

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    pdaf_cb.next_observation_pdaf = <void*>py__next_observation_pdaf
    cdef int  outflag
    with nogil:
        c__pdaf3_assimilate_ensrf(pdaf_cb.c__collect_state_pdaf,
                                  pdaf_cb.c__distribute_state_pdaf,
                                  pdaf_cb.c__init_dim_obs_pdaf,
                                  pdaf_cb.c__obs_op_pdaf,
                                  pdaf_cb.c__localize_covar_serial_pdaf,
                                  pdaf_cb.c__prepoststep_pdaf,
                                  pdaf_cb.c__next_observation_pdaf, &outflag)

    return outflag





def assim_offline_local(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__init_n_domains_p_pdaf, py__init_dim_l_pdaf,
    py__init_dim_obs_l_pdaf, py__g2l_state_pdaf, py__l2g_state_pdaf,
    py__prepoststep_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__init_n_domains_p_pdaf : Callable
        Provide number of local analysis domains

    py__init_dim_l_pdaf : Callable
        Init state dimension for local ana. domain

    py__init_dim_obs_l_pdaf : Callable
        Initialize local dimimension of obs. vector

    py__g2l_state_pdaf : Callable
        Get local state from full state

    py__l2g_state_pdaf : Callable
        Init full state from local state

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.init_n_domains_p_pdaf = <void*>py__init_n_domains_p_pdaf
    pdaf_cb.init_dim_l_pdaf = <void*>py__init_dim_l_pdaf
    pdaf_cb.init_dim_obs_l_pdaf = <void*>py__init_dim_obs_l_pdaf
    pdaf_cb.g2l_state_pdaf = <void*>py__g2l_state_pdaf
    pdaf_cb.l2g_state_pdaf = <void*>py__l2g_state_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_assim_offline_local(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__init_n_domains_p_pdaf,
                                     pdaf_cb.c__init_dim_l_pdaf,
                                     pdaf_cb.c__init_dim_obs_l_pdaf,
                                     pdaf_cb.c__g2l_state_pdaf,
                                     pdaf_cb.c__l2g_state_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_global(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__prepoststep_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_assim_offline_global(pdaf_cb.c__init_dim_obs_pdaf,
                                      pdaf_cb.c__obs_op_pdaf,
                                      pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_lenkf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__localize_covar_pdaf, py__prepoststep_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__localize_covar_pdaf : Callable
        Apply localization to HP and HPH^T

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_pdaf = <void*>py__localize_covar_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_assim_offline_lenkf(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__localize_covar_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag


def assim_offline_ensrf(py__init_dim_obs_pdaf, py__obs_op_pdaf,
    py__localize_covar_serial_pdaf, py__prepoststep_pdaf, int  outflag):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    py__init_dim_obs_pdaf : Callable
        Initialize dimension of full observation vector

    py__obs_op_pdaf : Callable
        Full observation operator

    py__localize_covar_serial_pdaf : Callable
        Apply localization to HP and BXY for single observation

    py__prepoststep_pdaf : Callable
        User supplied pre/poststep routine

    outflag : int
        Status flag

    Returns
    -------
    outflag : int
        Status flag
    """
    pdaf_cb.init_dim_obs_pdaf = <void*>py__init_dim_obs_pdaf
    pdaf_cb.obs_op_pdaf = <void*>py__obs_op_pdaf
    pdaf_cb.localize_covar_serial_pdaf = <void*>py__localize_covar_serial_pdaf
    pdaf_cb.prepoststep_pdaf = <void*>py__prepoststep_pdaf
    with nogil:
        c__pdaf3_assim_offline_ensrf(pdaf_cb.c__init_dim_obs_pdaf,
                                     pdaf_cb.c__obs_op_pdaf,
                                     pdaf_cb.c__localize_covar_serial_pdaf,
                                     pdaf_cb.c__prepoststep_pdaf, &outflag)

    return outflag
