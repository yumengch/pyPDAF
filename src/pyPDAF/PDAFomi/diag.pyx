import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3


def diag_dimobs():
    """diag_dimobs() -> np.ndarray

    Observation dimension for each observation type.

    Returns
    -------
    dim_obs: np.ndarray
        Observation dimension for each observation type. shape: (n_obs,)
    """
    cdef CFI_cdesc_rank1 dim_obs_ptr_cfi
    with nogil:
        c__pdafomi_diag_dimobs(<CFI_cdesc_t *> &dim_obs_ptr_cfi)

    cdef CFI_index_t dim_obs_ptr_subscripts[1]
    dim_obs_ptr_subscripts[0] = dim_obs_ptr_cfi.dim[0].lower_bound
    cdef int * dim_obs_ptr_ptr_np
    dim_obs_ptr_ptr_np = <int *>CFI_address(<CFI_cdesc_t *> &dim_obs_ptr_cfi, dim_obs_ptr_subscripts)
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="fortran", negative_indices=False, cast=False] dim_obs_ptr_np = np.asarray(<int [:dim_obs_ptr_cfi.dim[0].extent:1]> dim_obs_ptr_ptr_np, order="F")
    return dim_obs_ptr_np


def diag_get_hx(int  id_obs):
    """diag_get_hx(id_obs: int) -> Tuple[int, np.ndarray]

    Observed ensemble for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    hx_p_ptr : ndarray[np.float64, ndim=2]
        Pointer to observed ensemble mean
        Array shape: (dim_obs_p_diag, dim_ens)
    """
    cdef CFI_cdesc_rank2 hx_p_ptr_cfi
    cdef int  dim_obs_diag
    with nogil:
        c__pdafomi_diag_get_hx(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &hx_p_ptr_cfi)

    cdef CFI_index_t hx_p_ptr_subscripts[2]
    hx_p_ptr_subscripts[0] = hx_p_ptr_cfi.dim[0].lower_bound
    hx_p_ptr_subscripts[1] = hx_p_ptr_cfi.dim[1].lower_bound
    cdef double * hx_p_ptr_ptr_np
    hx_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &hx_p_ptr_cfi, hx_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hx_p_ptr_np = np.asarray(<double [:hx_p_ptr_cfi.dim[0].extent:1,:hx_p_ptr_cfi.dim[1].extent]> hx_p_ptr_ptr_np, order="F")
    return dim_obs_diag, hx_p_ptr_np


def diag_get_hxmean(int  id_obs):
    """diag_get_hxmean(id_obs: int) -> Tuple[int, np.ndarray]

    Observed ensemble mean for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    hxmean_p_ptr : ndarray[np.float64, ndim=1]
        Pointer to observed ensemble mean
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 hxmean_p_ptr_cfi
    cdef int  dim_obs_diag
    with nogil:
        c__pdafomi_diag_get_hxmean(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &hxmean_p_ptr_cfi)

    cdef CFI_index_t hxmean_p_ptr_subscripts[1]
    hxmean_p_ptr_subscripts[0] = hxmean_p_ptr_cfi.dim[0].lower_bound
    cdef double * hxmean_p_ptr_ptr_np
    hxmean_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &hxmean_p_ptr_cfi, hxmean_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxmean_p_ptr_np = np.asarray(<double [:hxmean_p_ptr_cfi.dim[0].extent:1]> hxmean_p_ptr_ptr_np, order="F")
    return dim_obs_diag, hxmean_p_ptr_np


def diag_get_ivar(int  id_obs):
    """diag_get_ivar(id_obs: int) -> Tuple[int, np.ndarray]

    Inverse of observation error variance for given observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    ivar_ptr : ndarray[np.float64, ndim=1]
        Pointer to inverse observation error variances
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 ivar_ptr_cfi
    cdef int  dim_obs_diag
    with nogil:
        c__pdafomi_diag_get_ivar(&id_obs, &dim_obs_diag, <CFI_cdesc_t *> &ivar_ptr_cfi)

    cdef CFI_index_t ivar_ptr_subscripts[1]
    ivar_ptr_subscripts[0] = ivar_ptr_cfi.dim[0].lower_bound
    cdef double * ivar_ptr_ptr_np
    ivar_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &ivar_ptr_cfi, ivar_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] ivar_ptr_np = np.asarray(<double [:ivar_ptr_cfi.dim[0].extent:1]> ivar_ptr_ptr_np, order="F")
    return dim_obs_diag, ivar_ptr_np


def diag_get_obs(int  id_obs):
    """diag_get_obs(id_obs: int) -> Tuple[int, int, np.ndarray, np.ndarray]

    Observation vector and corresponding coordinates for specified observation type.

    Parameters
    ----------
    id_obs : int
        Index of observation type to return

    Returns
    -------
    dim_obs_diag : int
        Observation dimension
    ncoord : int
        Number of observation dimensions
    obs_p_ptr : ndarray[np.float64, ndim=1]
        Pointer to observation vector
        Array shape: (:)
    ocoord_p_ptr : ndarray[np.float64, ndim=2]
        Pointer to Coordinate array
        Array shape: (:,:)
    """
    cdef CFI_cdesc_rank1 obs_p_ptr_cfi
    cdef CFI_cdesc_rank2 ocoord_p_ptr_cfi
    cdef int  dim_obs_diag
    cdef int  ncoord
    with nogil:
        c__pdafomi_diag_get_obs(&id_obs, &dim_obs_diag, &ncoord,
                                <CFI_cdesc_t *> &obs_p_ptr_cfi, <CFI_cdesc_t *> &ocoord_p_ptr_cfi)

    cdef CFI_index_t obs_p_ptr_subscripts[1]
    obs_p_ptr_subscripts[0] = obs_p_ptr_cfi.dim[0].lower_bound
    cdef double * obs_p_ptr_ptr_np
    obs_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &obs_p_ptr_cfi, obs_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] obs_p_ptr_np = np.asarray(<double [:obs_p_ptr_cfi.dim[0].extent:1]> obs_p_ptr_ptr_np, order="F")
    cdef CFI_index_t ocoord_p_ptr_subscripts[2]
    ocoord_p_ptr_subscripts[0] = ocoord_p_ptr_cfi.dim[0].lower_bound
    ocoord_p_ptr_subscripts[1] = ocoord_p_ptr_cfi.dim[1].lower_bound
    cdef double * ocoord_p_ptr_ptr_np
    ocoord_p_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &ocoord_p_ptr_cfi, ocoord_p_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] ocoord_p_ptr_np = np.asarray(<double [:ocoord_p_ptr_cfi.dim[0].extent:1,:ocoord_p_ptr_cfi.dim[1].extent]> ocoord_p_ptr_ptr_np, order="F")
    return dim_obs_diag, ncoord, obs_p_ptr_np, ocoord_p_ptr_np


def diag_nobstypes(int  nobs):
    """diag_nobstypes(nobs: int) -> int

    The number of observation types that are active in an assimilation run.

    Parameters
    ----------
    nobs : int
        Number of observation types (input can be an arbitrary number)

    Returns
    -------
    nobs : int
        Number of observation types
    """
    with nogil:
        c__pdafomi_diag_nobstypes(&nobs)

    return nobs


def diag_obs_rmsd(int  nobs, int  verbose):
    """diag_obs_rmsd(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    Root mean squared distance between observation and obseved model state
    for each observation type.

    Parameters
    ----------
    nobs : int
        Number of observation types
    verbose : int
        Verbosity flag

    Returns
    -------
    nobs : int
        Number of observation types
    rmsd_pointer : ndarray[np.float64, ndim=1]
        Vector of RMSD values
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 rmsd_pointer_cfi
    with nogil:
        c__pdafomi_diag_obs_rmsd(&nobs, <CFI_cdesc_t *> &rmsd_pointer_cfi, &verbose)

    cdef CFI_index_t rmsd_pointer_subscripts[1]
    rmsd_pointer_subscripts[0] = rmsd_pointer_cfi.dim[0].lower_bound
    cdef double * rmsd_pointer_ptr_np
    rmsd_pointer_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &rmsd_pointer_cfi, rmsd_pointer_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] rmsd_pointer_np = np.asarray(<double [:rmsd_pointer_cfi.dim[0].extent:1]> rmsd_pointer_ptr_np, order="F")
    return nobs, rmsd_pointer_np


def diag_stats(int  nobs, int  verbose):
    """diag_stats(nobs: int, verbose: int) -> Tuple[int, np.ndarray]

    A selection of 6 statistics comparing the observations and the
    observed ensemble mean for each observation type.

    Parameters
    ----------
    nobs : int
        Number of observation types
    verbose : int
        Verbosity flag

    Returns
    -------
    nobs : int
        Number of observation types
    obsstats_ptr : ndarray[np.float64, ndim=2]
        Array of observation statistics
        Included statistics are:
        - (1,:) correlations between observation and observed ensemble mean
        - (2,:) centered RMS difference between observation and observed ensemble mean
        - (3,:) mean bias (observation minus observed ensemble mean)
        - (4,:) mean absolute difference between observation and observed ensemble mean
        - (5,:) variance of observations
        - (6,:) variance of observed ensemble mean
        Array shape: (:,:)
    """
    cdef CFI_cdesc_rank2 obsstats_ptr_cfi
    with nogil:
        c__pdafomi_diag_stats(&nobs, <CFI_cdesc_t *> &obsstats_ptr_cfi, &verbose)

    cdef CFI_index_t obsstats_ptr_subscripts[2]
    obsstats_ptr_subscripts[0] = obsstats_ptr_cfi.dim[0].lower_bound
    obsstats_ptr_subscripts[1] = obsstats_ptr_cfi.dim[1].lower_bound
    cdef double * obsstats_ptr_ptr_np
    obsstats_ptr_ptr_np = <double *>CFI_address(<CFI_cdesc_t *> &obsstats_ptr_cfi, obsstats_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] obsstats_ptr_np = np.asarray(<double [:obsstats_ptr_cfi.dim[0].extent:1,:obsstats_ptr_cfi.dim[1].extent]> obsstats_ptr_ptr_np, order="F")
    return nobs, obsstats_ptr_np


