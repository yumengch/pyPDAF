"""Conversions between NEMO fields and compact PDAF state vectors.

NEMO stores variables on rank-local 2D/3D grids with land/halo points. PDAF
works on compact vectors. This module performs the packing/unpacking and
optional variable transformations described by each ``StateField``.
"""
import numpy as np
import xarray as xr

import log
import model

def field2state(field: xr.DataArray, wet_grid: model.WetGrid) -> np.ndarray:
    """Convert a DataArray field to a compact wet-point state vector.

    Supported field shapes are ``(nav_lev, lat, lon)`` and ``(lat, lon)``.
    This implements only the Fortran ``use_wet_state == 2`` branch.
    """
    if wet_grid.nwet2d == 0:
        return np.zeros(0)
    if field.ndim == 2:
        return field.values.ravel()[wet_grid.wet_pts[3]]
    if field.ndim == 3:
        # slice the horizontal coordinates for wet points
        # make vertical coordinate the last dimension so it changes first
        values_at_points = field.values.reshape(field.shape[0], -1)[:, wet_grid.wet_pts[3]].T
        # find the mask for wet points in the vertical dimension
        wet_mask_point_major = wet_grid.tmask[:, wet_grid.wet_pts[5], wet_grid.wet_pts[6]].T
        return values_at_points[wet_mask_point_major]

    raise ValueError(f"NEMO-PDAF: cannot handle {field.ndim} number of dimensions.")


def state2field(
    state: np.ndarray,
    field: xr.DataArray,
    wet_grid: model.WetGrid,
    offset: int,
) -> xr.DataArray:
    """Fill a model field from a compact wet-point state vector.

    ``field`` is expected to include ``time_counter`` as its first dimension,
    so the values are written into the first time record. After removing time,
    supported field shapes are ``(nav_lev, lat, lon)`` and ``(lat, lon)``.
    """
    values = field.values[0]

    if wet_grid.nwet2d == 0:
        return field

    if values.ndim == 2:
        flat = values.ravel()
        flat[wet_grid.wet_pts[3]] = state[offset:offset + wet_grid.nwet2d]
        return field

    if values.ndim == 3:
        field_values = values.reshape(values.shape[0], -1)
        point_values = field_values[:, wet_grid.wet_pts[3]].T
        wet_mask_point_major = wet_grid.tmask[:, wet_grid.wet_pts[5], wet_grid.wet_pts[6]].T
        point_values[wet_mask_point_major] = state[offset:offset + wet_grid.nwet3d]
        field_values[:, wet_grid.wet_pts[3]] = point_values.T
        return field

    raise ValueError(f"NEMO-PDAF: cannot handle {values.ndim} number of dimensions.")


def transform_field_mv(transform_type, ens_p, fields, limits=0):
    """Transform all state-vector fields in place.

    Transformations are configured per field in ``config.ini``. They are useful
    for positive-definite variables such as concentrations, thicknesses, or
    biogeochemical tracers. The transformation is applied before PDAF analysis
    and reversed before writing NEMO ASM files.

    Parameters
    ----------
    transform_type : int
        ``1`` transforms from NEMO values to transformed values.
        ``2`` transforms from transformed values back to NEMO values.
    state : ndarray
        State vector, modified in place.
    sfields : object, dict, or iterable
        Field metadata. Supports the ``Sfields`` container, a dict of fields,
        or an iterable of field objects/dicts.
    limits : int
        ``0`` no limits, ``11`` clip forecast in original units before
        transform, ``12`` clip forecast after transform, ``21`` clip analysis
        in original units after inverse transform, ``22`` clip analysis before
        inverse transform.
    """

    if transform_type == 1 and limits == 11:
        ens_p = var_limits_mv(ens_p, fields)
    elif transform_type == 2 and limits == 22:
        ens_p = var_limits_mv(ens_p, fields)

    log.info("Apply field transformations")

    for vname in fields:
        field = fields[vname]
        trafo = field.transform
        shift = field.trafo_shift
        variable = field.variable
        if trafo == 0:
            continue

        if transform_type == 1:
            if trafo == 1:
                log.info(f"--- apply log-10 transformation to {variable}")
                ens_p[field.off:field.off + field.dim] = np.log10(
                    ens_p[field.off:field.off + field.dim] + shift)
            elif trafo == 2:
                log.info(f"--- apply ln transformation to {variable}")
                positive = ens_p[field.off:field.off + field.dim] > 0.0
                n_nonpositive = int((~positive).sum())
                n_positive = int(positive.sum())
                ens_p[field.off:field.off + field.dim][positive] = np.log(
                    ens_p[field.off:field.off + field.dim][positive] + shift)
                ens_p[field.off:field.off + field.dim][~positive] = -30.0
                if n_nonpositive > 0:
                    log.info(
                        f"--- number of <=0, >0: {variable} "
                        f"{n_nonpositive:9d}{n_positive:9d}"
                    )
            else:
                raise ValueError(f"ERROR: no valid variable transformation selected, type {trafo}")

        elif transform_type == 2:
            if trafo == 1:
                log.info(f"--- revert log-10 transformation of {variable}")
                ens_p[field.off:field.off + field.dim] = (10.0 ** ens_p[field.off:field.off + field.dim]) - shift
            elif trafo == 2:
                log.info(f"--- revert ln transformation of {variable}")
                valid = ens_p[field.off:field.off + field.dim] > -30.0
                n_omitted = int((~valid).sum())
                n_valid = int(valid.sum())
                ens_p[field.off:field.off + field.dim][valid] = np.exp(ens_p[field.off:field.off + field.dim][valid]) - shift
                ens_p[field.off:field.off + field.dim][~valid] = 0.0
                if n_omitted > 0:
                    log.info(
                        f"--- number of <=0, >0: {variable} "
                        f"{n_omitted:9d}{n_valid:9d}"
                    )
            else:
                raise ValueError(f"ERROR: no valid variable transformation selected, type {trafo}")

    if transform_type == 1 and limits == 12:
        ens_p = var_limits_mv(ens_p, fields)
    elif transform_type == 2 and limits == 21:
        ens_p = var_limits_mv(ens_p, fields)

    return ens_p


def var_limits_mv(ens_p:np.ndarray, fields:dict) -> np.ndarray:
    """Apply configured min/max limits to each state-vector field in place.

    Per-field ``limit`` options are: ``0`` no clipping, ``1`` minimum only,
    ``2`` maximum only, and ``3`` both minimum and maximum. Use this for
    physical bounds such as non-negative concentrations or sea-ice fractions.
    """

    log.info("Apply limits to model fields")

    for vname in fields:
        field = fields[vname]

        if field.limit == 1:
            log.info(f"--- apply min. limit {field.min_limit:12.3e} to {field.variable}")
            affected = ens_p[field.off:field.off + field.dim] < field.min_limit
            ens_p[field.off:field.off + field.dim][affected] = field.min_limit
            count = int(affected.sum())
            log.info(f"--- number of affected values {count:8d}")
        elif field.limit == 2:
            log.info(f"--- apply max. limit {field.max_limit:12.3e} to {field.variable}")
            affected = ens_p[field.off:field.off + field.dim] > field.max_limit
            ens_p[field.off:field.off + field.dim][affected] = field.max_limit
            count = int(affected.sum())
            log.info(f"--- number of affected values {count:8d}")

        elif field.limit == 3:
            log.info(
                    f"--- apply min/max limits of {field.min_limit:12.3e}"
                    f"{field.max_limit:12.3e} to {field.variable}"
                )
            below = ens_p[field.off:field.off + field.dim] < field.min_limit
            above = ens_p[field.off:field.off + field.dim] > field.max_limit
            ens_p[field.off:field.off + field.dim][below] = field.min_limit
            ens_p[field.off:field.off + field.dim][above] = field.max_limit
            log.info(
                f"--- number of affected values "
                f"{int(below.sum()):8d}{int(above.sum()):8d}"
            )

        else:
            log.info(f"--- no limit applied to {field.variable}")

    return ens_p
