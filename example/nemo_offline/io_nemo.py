"""Read NEMO input files and write NEMO ASM output files.

The reader expects one restart file per ensemble member and MPI/filter rank.
The writer creates the background-state and increment files consumed by NEMO's
ASM/IAU machinery. For another model, this module is the primary adapter:
replace the path conventions, NetCDF variable names, and output format while
leaving PDAF state-vector callbacks unchanged.
"""
import configparser
import pathlib
import typing

import numpy as np
import xarray as xr

import log
import model
import parallel
import sfields
import transforms

class Reader:
    """Reading domain information and restart state from NEMO files.

    Parameters
    ----------
    config : dict
        Dictionary holding IO settings:
          - path_rst_root: path to the ensemble restart files (typically run directory)
          - ens_prefix: prefix for ensemble member directories,
                        experiment for each member is in:  ens_prefix{i_member}
          - path_rst_suffix: directory between ens_prefix{i_member} and the restart files
                            e.g., ens_prefix{i_member}/path_rst_suffix/restart_files
          - f_basename_rst: basename of restart files
          - path_dom: path to the domain file
          - fname_dom: filename of the domain file
    pe : parallel.Parallel
        Object with ``mype_filter`` and ``npes_filter``.
    """

    def __init__(self, pe: parallel.Parallel, config:configparser.SectionProxy):
        self.pe = pe

        fname_dom = config.get("fname_dom", "")
        path_dom = pathlib.Path(config.get("path_dom", "."))
        self.path_dom = path_dom / fname_dom
        self.f_basename_rst = config.get("f_basename_rst", "")
        self.path_rst_root = pathlib.Path(config.get("path_rst_root", "."))
        self.ens_prefix = config.get("ens_prefix", "ens_")
        self.path_rst_suffix = pathlib.Path(config.get("path_rst_suffix", ""))

    # def print_io_configuration(self):
    #     """Print the IO configuration to the log."""
    #     path_rst = self.path_rst_root / f"{self.ens_prefix}1{self.path_rst_suffix}"
    #     log.info("  [io_nml]:")
    #     info = f"    path to model domain file          {self.path_dom}"
    #     log.info(info)
    #     info = f"    model domain filename              {self.fname_dom}"
    #     log.info(info)
    #     info = f"    path to restart files              {path_rst}"
    #     log.info(info)
    #     info = f"    basename of restart files          {self.f_basename_rst}"
    #     log.info(info)
    #     info = f"    path to write increment files      {self.path_asm_root}"
    #     log.info(info)

    def read_local_domain_coords(self, vnames: list[str], attnames:list[str]
                                 ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """Read process-local domain metadata from the first ensemble restart."""

        # get restart file path
        path_rst = self.path_rst_root / f"{self.ens_prefix}1" / self.path_rst_suffix
        fname = f"{self.f_basename_rst}_{str(self.pe.mype_filter).zfill(4)}.nc"

        log.info("Reading local domain information from file")
        info = f"   Input file: {path_rst / fname}"
        log.info(info)

        # read domain metadata from restart file
        with xr.open_dataset(path_rst / fname) as ds:
            attrs = {attname: ds.attrs[attname] for attname in attnames}
            variables = {vname: ds[vname].values for vname in vnames}

        # if self.pe.npes_model > 1 and self.verbose_io > 0:
            # info = (
            #     "NEMO-PDAF  RANK"
            #     f"{self.pe.npes_model:6d} {local_domain_attrs['i0']:7d}"
            #     f"{local_domain_attrs['i0'] + local_domain_attrs['ni_p'] - 1:7d}"
            #     f"{local_domain_attrs['j0']:7d}"
            #     f"{local_domain_attrs['j0'] + local_domain_attrs['nj_p'] - 1:7d}"
            #     f"{local_domain_attrs['ni_p']:7d}{local_domain_attrs['nj_p']:7d}"
            # )
            # self.logger.info(info)
        # self.pe.comm_model.Barrier()
        return variables, attrs

    def read_global_domain(self, vnames: list[str], indexers: dict[str, slice]
                           ) -> tuple[dict[str, int], dict[str, np.ndarray]]:
        """Read global grid variables and create the local ``tmask``."""
        log.info("Reading grid variables from file")
        info = f"   Input file: {self.path_dom}"
        log.info(info)

        with xr.open_dataset(self.path_dom) as ds:
            variables = {vname: ds[vname].isel(indexers).values for vname in vnames}
            dims = {str(dim): int(ds.sizes[dim]) for dim in ds.dims}

        # if self.pe.mype_model == 0 and self.verbose_io > 0:
        #     jpiglo = global_domain_attrs["jpiglo"]
        #     jpjglo = global_domain_attrs["jpjglo"]
        #     jpk = global_domain_attrs["jpk"]
        #     self.logger.info("NEMO-PDAF     *** NEMO: grid dimensions ***")
        #     info = f"NEMO-PDAF   {jpiglo:12d}{jpjglo:12d}{jpk:12d}"
        #     self.logger.info(info)
        #     info = f"NEMO-PDAF     Dimension of global 3D grid box{jpiglo * jpjglo * jpk:12d}"
        #     self.logger.info(info)
        #     info = f"NEMO-PDAF     Number of global surface points{jpiglo * jpjglo:12d}"
        #     self.logger.info(info)

        # self.pe.comm_model.Barrier()
        return dims, variables

    def read_restart(self, i_ens:int, fields: dict[str, sfields.StateField]
                     ) -> typing.Iterator[tuple[xr.DataArray, int, int, str]]:
        """Read restart fields for one ensemble member into ``state_p``."""

        log.info("*** Ensemble: Reading model restart file")
        for field in fields.values():
            path_rst = self.path_rst_root / f"{self.ens_prefix}{i_ens}" / self.path_rst_suffix
            fname = f"{self.f_basename_rst}_{str(self.pe.mype_filter).zfill(4)}.nc"
            log.info(f"Reading: {path_rst / fname}")
            log.info(f"{field.variable:>12},  offset{field.off:10d}")

            with xr.open_dataset(path_rst / fname) as ds:
                data = ds[field.name_rest_n].squeeze()

            if field.k_name == "cat":
                if field.operation == "sum_over_cat_ice":
                    data = data.sum(dim="cat", skipna=True)
                elif field.operation == "select_cat_ice":
                    data = data.isel(cat=field.ice_cat - 1)
                else:
                    log.warning(
                        "Unknown operation for combining variables in restart file: "
                        f"{field.operation}. Using first name_rest_n only."
                    )

            if data.ndim > field.ndims:
                data = data.isel(nav_lev=0)

            yield data, field.off, field.dim, field.variable


class Writer:
    """Writing increments and background state for NEMO-PDAF.

    They are meant to be used with ASM module in NEMO.

    Parameters
    ----------
    pe : parallel.Parallel
        Object with ``mype_filter`` and ``npes_filter``.
    config : dict
        Dictionary holding IO settings:
            - do_deflate: whether to deflate the output NetCDF files
            - ens_prefix: prefix for ensemble member directories
            - path_asm_root: path to the directory where increment and
                            background state files are written
    """
    def __init__(self, pe: parallel.Parallel, config:configparser.SectionProxy):
        self.pe = pe

        self.do_deflate = config.getboolean("do_deflate", False)
        self.ens_prefix = config.get("ens_prefix", "ens_")
        self.path_asm_root = pathlib.Path(config.get("path_asm_root", "."))

    def write_asminc_mv(
        self,
        ens_member: int,
        increment: np.ndarray,
        fields: dict[str, sfields.StateField],
        nemo_grid: model.NemoDomain,
    ):
        """Write NEMO ASM increment file compatible with IAU."""
        log.info("        --- Write increment file")

        time = nemo_grid.spatial_time.ndastp + nemo_grid.spatial_time.nn_time0 * 0.0001
        time_counter = np.atleast_1d(nemo_grid.spatial_time.time_counter)
        path = (self.path_asm_root / f"{self.ens_prefix}{ens_member}" /
                f"assim_background_increments_{self.pe.mype_filter:04d}.nc"
        )

        variables = {}
        for field in fields.values():
            data = self._empty_output_field(field, nemo_grid, time_counter, field.name_incr)
            transforms.state2field(increment, data, nemo_grid.wet_grid, field.off)
            variables[data.name] = data

        ds = xr.Dataset(variables,
                        attrs = self._domain_attrs(nemo_grid,
                                                   "Increment for NEMO-PDAF data assimilation"))
        ds["time"] = xr.DataArray(time)
        ds["z_inc_dateb"] = xr.DataArray(time)
        ds["z_inc_datef"] = xr.DataArray(time)
        encoding = {}
        if self.do_deflate:
            encoding = {
                name: {"zlib": True, "complevel": 1}
                for name in ds.data_vars
            }
        log.info(f"Create file: {path}")
        ds.to_netcdf(path, encoding=encoding, unlimited_dims=["time_counter"])

    def write_asmdin_mv(
        self,
        ens_member: int,
        state_f: np.ndarray,
        fields: dict[str, sfields.StateField],
        nemo_grid: model.NemoDomain,
    ):
        """Write NEMO ASM background-state file compatible with direct initialization."""
        log.info("        --- Write background state file")
        path = (self.path_asm_root / f"{self.ens_prefix}{ens_member}" /
                f"assim_background_state_DI_{self.pe.mype_filter:04d}.nc"
        )
        time_counter = np.atleast_1d(nemo_grid.spatial_time.time_counter)

        variables = {}
        for field in fields.values():
            data = self._empty_output_field(field, nemo_grid, time_counter,
                field.name_bkg_din,
            )
            transforms.state2field(state_f, data, nemo_grid.wet_grid, field.off)
            variables[data.name] = data

        ds = xr.Dataset(variables,
                        attrs = self._domain_attrs(nemo_grid,
                                                   "background for NEMO-PDAF data assimilation"))
        ds["rdastp"] = xr.DataArray(nemo_grid.spatial_time.ndastp)
        encoding = {}
        if self.do_deflate:
            encoding = {
                name: {"zlib": True, "complevel": 1}
                for name in ds.data_vars
            }
        log.info(f"Create file: {path}")
        ds.to_netcdf(path, encoding=encoding, unlimited_dims=["time_counter"])

    def _domain_attrs(self, nemo_grid, title):
        local_dom = nemo_grid.local_dom
        return {
            "title": title,
            "DOMAIN_number_total": self.pe.npes_filter,
            "DOMAIN_number": self.pe.mype_filter,
            "DOMAIN_dimensions_ids": np.array([1, 2], dtype=np.int32),
            "DOMAIN_size_global": np.array(
                [nemo_grid.spatial_time.jpiglo, nemo_grid.spatial_time.jpjglo],
                dtype=np.int32,
            ),
            "DOMAIN_size_local": np.array([local_dom.ni_p, local_dom.nj_p], dtype=np.int32),
            "DOMAIN_position_first": np.array([local_dom.i0, local_dom.j0], dtype=np.int32),
            "DOMAIN_position_last": np.array(
                [local_dom.i0 + local_dom.ni_p - 1, local_dom.j0 + local_dom.nj_p - 1],
                dtype=np.int32,
            ),
            "DOMAIN_halo_size_start": np.asarray(local_dom.halo0, dtype=np.int32),
            "DOMAIN_halo_size_end": np.asarray(local_dom.halo1, dtype=np.int32),
            "DOMAIN_type": "BOX",
        }

    def _empty_output_field(
        self,
        field: sfields.StateField,
        nemo_grid: model.NemoDomain,
        time_counter: np.ndarray,
        name: str,
    ) -> xr.DataArray:
        """Create an output DataArray; values are filled by ``transforms.state2field``."""
        coords = {
            "y": np.arange(nemo_grid.local_dom.nj_p),
            "x": np.arange(nemo_grid.local_dom.ni_p),
            "time_counter": np.atleast_1d(time_counter),
            "nav_lon": (("y", "x"), nemo_grid.coords.nav_lon),
            "nav_lat": (("y", "x"), nemo_grid.coords.nav_lat),
        }
        attrs = {"long_name": field.variable}
        if field.unit:
            attrs["units"] = field.unit

        dims: tuple[str, ...]
        shape: tuple[int, ...]
        if field.ndims == 2:
            shape = (1, nemo_grid.local_dom.nj_p, nemo_grid.local_dom.ni_p)
            dims = ("time_counter", "y", "x",)
        elif field.ndims == 3:
            shape = (
                1,
                nemo_grid.local_dom.nk_p,
                nemo_grid.local_dom.nj_p,
                nemo_grid.local_dom.ni_p,
            )
            dims = ("time_counter", "nav_lev", "y", "x")
            coords["nav_lev"] = nemo_grid.coords.nav_lev
        else:
            raise ValueError(f"NEMO-PDAF: cannot handle {field.ndims} number of dimensions.")

        return xr.DataArray(
            np.zeros(shape, dtype=float),
            dims=dims,
            coords=coords,
            name=name,
            attrs=attrs,
        )
