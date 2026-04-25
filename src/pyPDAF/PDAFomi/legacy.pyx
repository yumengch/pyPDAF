import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

def localize_covar_iso(int  i_obs, int  dim, int  locweight,
    double  cradius, double  sradius, double [::1,:] coords,
    double [::1,:] hp, double [::1,:] hph):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    locweight : int
        Localization weight type
    cradius : double
        localization radius
    sradius : double
        support radius for weight functions
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)

    Returns
    -------
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[2]
    hp_extent[0] = hp.shape[0]
    hp_extent[1] = hp.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 hph_cfi
    cdef size_t hph_nbytes = hph.nbytes
    cdef CFI_index_t hph_extent[2]
    hph_extent[0] = hph.shape[0]
    hph_extent[1] = hph.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hph_np = np.asarray(hph, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                      CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0,0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 2, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hph_cfi, &hph[0,0], CFI_attribute_other,
                    CFI_type_double , hph_nbytes, 2, hph_extent)

    c__pdafomi_localize_covar_iso(&i_obs, &dim, &locweight, &cradius,
                                    &sradius, <CFI_cdesc_t *> &coords_cfi,
                                    <CFI_cdesc_t *> &hp_cfi, <CFI_cdesc_t *> &hph_cfi)

    return hp_np, hph_np


def localize_covar_noniso_locweights(int  i_obs, int  dim,
    int [::1] locweights, double [::1] cradius, double [::1] sradius,
    double [::1,:] coords, double [::1,:] hp, double [::1,:] hph):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    locweights : ndarray[np.intc, ndim=1]
        Types of localization function
        Array shape: (:)
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)

    Returns
    -------
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[2]
    hp_extent[0] = hp.shape[0]
    hp_extent[1] = hp.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 hph_cfi
    cdef size_t hph_nbytes = hph.nbytes
    cdef CFI_index_t hph_extent[2]
    hph_extent[0] = hph.shape[0]
    hph_extent[1] = hph.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hph_np = np.asarray(hph, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 locweights_cfi
    cdef size_t locweights_nbytes = locweights.nbytes
    cdef CFI_index_t locweights_extent[1]
    locweights_extent[0] = locweights.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &locweights_cfi, &locweights[0], CFI_attribute_other,
                      CFI_type_int , locweights_nbytes, 1, locweights_extent)

    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                    CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                    CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0,0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 2, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hph_cfi, &hph[0,0], CFI_attribute_other,
                    CFI_type_double , hph_nbytes, 2, hph_extent)

    c__pdafomi_localize_covar_noniso_locweights(&i_obs, &dim,
                                                <CFI_cdesc_t *> &locweights_cfi,
                                                <CFI_cdesc_t *> &cradius_cfi,
                                                <CFI_cdesc_t *> &sradius_cfi,
                                                <CFI_cdesc_t *> &coords_cfi,
                                                <CFI_cdesc_t *> &hp_cfi,
                                                <CFI_cdesc_t *> &hph_cfi)

    return hp_np, hph_np


def localize_covar_noniso(int  i_obs, int  dim, int  locweight,
    double [::1] cradius, double [::1] sradius, double [::1,:] coords,
    double [::1,:] hp, double [::1,:] hph):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    dim : int
        State dimension
    locweight : int
        Localization weight type
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)

    Returns
    -------
    hp : ndarray[np.float64, ndim=2]
        Matrix HP, dimension (nobs, dim)
        Array shape: (:, :)
    hph : ndarray[np.float64, ndim=2]
        Matrix HPH, dimension (nobs, nobs)
        Array shape: (:, :)
    """
    cdef CFI_cdesc_rank2 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[2]
    hp_extent[0] = hp.shape[0]
    hp_extent[1] = hp.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 hph_cfi
    cdef size_t hph_nbytes = hph.nbytes
    cdef CFI_index_t hph_extent[2]
    hph_extent[0] = hph.shape[0]
    hph_extent[1] = hph.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] hph_np = np.asarray(hph, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                    CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0,0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 2, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hph_cfi, &hph[0,0], CFI_attribute_other,
                    CFI_type_double , hph_nbytes, 2, hph_extent)

    c__pdafomi_localize_covar_noniso(&i_obs, &dim, &locweight,
                                        <CFI_cdesc_t *> &cradius_cfi,
                                        <CFI_cdesc_t *> &sradius_cfi,
                                        <CFI_cdesc_t *> &coords_cfi,
                                        <CFI_cdesc_t *> &hp_cfi,
                                        <CFI_cdesc_t *> &hph_cfi)

    return hp_np, hph_np


def localize_covar_serial_iso(int  i_obs, int  iobs_all, int  dim,
    int  dim_obs, int  locweight, double  cradius, double  sradius,
    double [::1,:] coords, double [::1] hp, double [::1] hxy):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    iobs_all : int
        Index of current observation
    dim : int
        State dimension
    dim_obs : int
        Overall full observation dimension
    locweight : int
        Localization weight type
    cradius : double
        localization radius
    sradius : double
        support radius for weight functions
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)

    Returns
    -------
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[1]
    hp_extent[0] = hp.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 hxy_cfi
    cdef size_t hxy_nbytes = hxy.nbytes
    cdef CFI_index_t hxy_extent[1]
    hxy_extent[0] = hxy.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxy_np = np.asarray(hxy, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                      CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 1, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hxy_cfi, &hxy[0], CFI_attribute_other,
                    CFI_type_double , hxy_nbytes, 1, hxy_extent)

    c__pdafomi_localize_covar_serial_iso(&i_obs, &iobs_all, &dim,
                                            &dim_obs, &locweight,
                                            &cradius, &sradius,
                                            <CFI_cdesc_t *> &coords_cfi,
                                            <CFI_cdesc_t *> &hp_cfi,
                                            <CFI_cdesc_t *> &hxy_cfi)

    return hp_np, hxy_np


def localize_covar_serial_noniso_locweights(int  i_obs, int  iobs_all,
    int  dim, int  dim_obs, int [::1] locweights, double [::1] cradius,
    double [::1] sradius, double [::1,:] coords, double [::1] hp,
    double [::1] hxy):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    iobs_all : int
        Index of current observation
    dim : int
        State dimension
    dim_obs : int
        Overall full observation dimension
    locweights : ndarray[np.intc, ndim=1]
        Types of localization function
        Array shape: (:)
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)

    Returns
    -------
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[1]
    hp_extent[0] = hp.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 hxy_cfi
    cdef size_t hxy_nbytes = hxy.nbytes
    cdef CFI_index_t hxy_extent[1]
    hxy_extent[0] = hxy.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxy_np = np.asarray(hxy, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 locweights_cfi
    cdef size_t locweights_nbytes = locweights.nbytes
    cdef CFI_index_t locweights_extent[1]
    locweights_extent[0] = locweights.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &locweights_cfi, &locweights[0], CFI_attribute_other,
                      CFI_type_int , locweights_nbytes, 1, locweights_extent)

    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                    CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                    CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 1, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hxy_cfi, &hxy[0], CFI_attribute_other,
                    CFI_type_double , hxy_nbytes, 1, hxy_extent)

    c__pdafomi_localize_covar_serial_noniso_locweights(&i_obs,
                                                        &iobs_all, &dim,
                                                        &dim_obs,
                                                        <CFI_cdesc_t *> &locweights_cfi,
                                                        <CFI_cdesc_t *> &cradius_cfi,
                                                        <CFI_cdesc_t *> &sradius_cfi,
                                                        <CFI_cdesc_t *> &coords_cfi,
                                                        <CFI_cdesc_t *> &hp_cfi,
                                                        <CFI_cdesc_t *> &hxy_cfi)
    return hp_np, hxy_np


def localize_covar_serial_noniso(int  i_obs, int  iobs_all, int  dim,
    int  dim_obs, int  locweight, double [::1] cradius,
    double [::1] sradius, double [::1,:] coords, double [::1] hp,
    double [::1] hxy):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    iobs_all : int
        Index of current observation
    dim : int
        State dimension
    dim_obs : int
        Overall full observation dimension
    locweight : int
        Localization weight type
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    coords : ndarray[np.float64, ndim=2]
        Coordinates of state vector elements
        Array shape: (:,:)
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)

    Returns
    -------
    hp : ndarray[np.float64, ndim=1]
        Vector HP, dimension (dim)
        Array shape: (:)
    hxy : ndarray[np.float64, ndim=1]
        Matrix HXY, dimension (nobs)
        Array shape: (:)
    """
    cdef CFI_cdesc_rank1 hp_cfi
    cdef size_t hp_nbytes = hp.nbytes
    cdef CFI_index_t hp_extent[1]
    hp_extent[0] = hp.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hp_np = np.asarray(hp, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 hxy_cfi
    cdef size_t hxy_nbytes = hxy.nbytes
    cdef CFI_index_t hxy_extent[1]
    hxy_extent[0] = hxy.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] hxy_np = np.asarray(hxy, dtype=np.float64, order="F")
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    cdef CFI_cdesc_rank2 coords_cfi
    cdef size_t coords_nbytes = coords.nbytes
    cdef CFI_index_t coords_extent[2]
    coords_extent[0] = coords.shape[0]
    coords_extent[1] = coords.shape[1]
    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                      CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    CFI_establish(<CFI_cdesc_t *> &coords_cfi, &coords[0,0], CFI_attribute_other,
                    CFI_type_double , coords_nbytes, 2, coords_extent)

    CFI_establish(<CFI_cdesc_t *> &hp_cfi, &hp[0], CFI_attribute_other,
                    CFI_type_double , hp_nbytes, 1, hp_extent)

    CFI_establish(<CFI_cdesc_t *> &hxy_cfi, &hxy[0], CFI_attribute_other,
                    CFI_type_double , hxy_nbytes, 1, hxy_extent)

    c__pdafomi_localize_covar_serial_noniso(&i_obs, &iobs_all, &dim,
                                            &dim_obs, &locweight,
                                            <CFI_cdesc_t *> &cradius_cfi,
                                            <CFI_cdesc_t *> &sradius_cfi,
                                            <CFI_cdesc_t *> &coords_cfi,
                                            <CFI_cdesc_t *> &hp_cfi,
                                            <CFI_cdesc_t *> &hxy_cfi)

    return hp_np, hxy_np


def init_dim_obs_l_iso_old(int  i_obs, double [::1] coords_l,
    int  locweight, double  cradius, double  sradius, int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweight : int
        Type of localization function
    cradius : double
        Localization cut-off radius (single or vector)
    sradius : double
        Support radius of localization function (single or vector)
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coords_l_cfi, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

    c__pdafomi_init_dim_obs_l_iso_old(&i_obs, <CFI_cdesc_t *> &coords_l_cfi, &locweight,
                                        &cradius, &sradius, &cnt_obs_l)

    return cnt_obs_l


def init_dim_obs_l_noniso_old(int  i_obs, double [::1] coords_l,
    int  locweight, double [::1] cradius, double [::1] sradius, int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweight : int
        Type of localization function
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coords_l_cfi, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                    CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    c__pdafomi_init_dim_obs_l_noniso_old(&i_obs, <CFI_cdesc_t *> &coords_l_cfi,
                                            &locweight, <CFI_cdesc_t *> &cradius_cfi,
                                            <CFI_cdesc_t *> &sradius_cfi, &cnt_obs_l)

    return cnt_obs_l


def init_dim_obs_l_noniso_locweights_old(int  i_obs, double [::1] coords_l,
    int [::1] locweights, double [::1] cradius, double [::1] sradius,
    int  cnt_obs_l):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int
        index into observation arrays
    coords_l : ndarray[np.float64, ndim=1]
        Coordinates of current analysis domain
        Array shape: (:)
    locweights : ndarray[np.intc, ndim=1]
        Types of localization function
        Array shape: (:)
    cradius : ndarray[np.float64, ndim=1]
        Vector of localization cut-off radii
        Array shape: (:)
    sradius : ndarray[np.float64, ndim=1]
        Vector of support radii of localization function
        Array shape: (:)
    cnt_obs_l : int
        Local dimension of current observation vector

    Returns
    -------
    cnt_obs_l : int
        Local dimension of current observation vector
    """
    cdef CFI_cdesc_rank1 coords_l_cfi
    cdef size_t coords_l_nbytes = coords_l.nbytes
    cdef CFI_index_t coords_l_extent[1]
    coords_l_extent[0] = coords_l.shape[0]
    cdef CFI_cdesc_rank1 locweights_cfi
    cdef size_t locweights_nbytes = locweights.nbytes
    cdef CFI_index_t locweights_extent[1]
    locweights_extent[0] = locweights.shape[0]
    cdef CFI_cdesc_rank1 cradius_cfi
    cdef size_t cradius_nbytes = cradius.nbytes
    cdef CFI_index_t cradius_extent[1]
    cradius_extent[0] = cradius.shape[0]
    cdef CFI_cdesc_rank1 sradius_cfi
    cdef size_t sradius_nbytes = sradius.nbytes
    cdef CFI_index_t sradius_extent[1]
    sradius_extent[0] = sradius.shape[0]
    CFI_establish(<CFI_cdesc_t *> &coords_l_cfi, &coords_l[0], CFI_attribute_other,
                      CFI_type_double , coords_l_nbytes, 1, coords_l_extent)

    CFI_establish(<CFI_cdesc_t *> &locweights_cfi, &locweights[0], CFI_attribute_other,
                    CFI_type_int , locweights_nbytes, 1, locweights_extent)

    CFI_establish(<CFI_cdesc_t *> &cradius_cfi, &cradius[0], CFI_attribute_other,
                    CFI_type_double , cradius_nbytes, 1, cradius_extent)

    CFI_establish(<CFI_cdesc_t *> &sradius_cfi, &sradius[0], CFI_attribute_other,
                    CFI_type_double , sradius_nbytes, 1, sradius_extent)

    c__pdafomi_init_dim_obs_l_noniso_locweights_old(&i_obs,
                                                    <CFI_cdesc_t *> &coords_l_cfi,
                                                    <CFI_cdesc_t *> &locweights_cfi,
                                                    <CFI_cdesc_t *> &cradius_cfi,
                                                    <CFI_cdesc_t *> &sradius_cfi,
                                                    &cnt_obs_l)

    return cnt_obs_l


def deallocate_obs(int  i_obs):
    r"""Deallocate OMI-internal obsrevation arrays

    This function should not be called by users
    because it is called internally in PDAF.

    Parameters
    ----------
    i_obs : int
        index of observations
    """
    c__pdafomi_deallocate_obs(&i_obs)



