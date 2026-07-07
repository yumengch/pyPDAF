"""State-vector field definitions and offsets for pyPDAF.

The state vector is configured from ``config.ini`` rather than hard-coded in
the PDAF callbacks. Each listed field is assigned a local wet-grid dimension
and offset so NEMO variables can be packed into one PDAF vector.
"""
import configparser
import dataclasses
import typing

import mpi4py.MPI

import log
import parallel
import model


@dataclasses.dataclass
class StateField:
    """Metadata for one model variable inside the PDAF state vector.

    Each configured state field maps a named model variable to a contiguous
    block in the PDAF state. The required config keys for a normal NEMO field
    are ``ndims``, ``variable``, ``name_rest_n``, ``name_incr``, and
    ``name_bkg_din``. The remaining keys are optional hooks for transforms,
    clipping, units, and category handling.

    Attributes
    ----------
    ndims
        ``2`` for surface/horizontal wet points, ``3`` for all wet
        water-column points.
    dim
        PE-local length of the model variable assigned after reading the wet mask.
    off
        Zero-based offset of this field in the PE-local state vector.
    variable
        Descriptive variable name used in logs and NetCDF attributes.
    name_incr
        Variable name written to NEMO ASM increment files.
    name_rest_n
        Variable name read from NEMO restart files.
    name_bkg_din
        Variable name written to NEMO ASM background-state files.
    rst_file
        Per-field restart basename metadata. The current reader uses the
        global ``[io].f_basename_rst`` for all fields.
    k_name
        Vertical/category axis meaning. Use ``lev`` for normal vertical
        levels or ``cat`` for category-based fields.
    operation
        Category handling when ``k_name == "cat"``:
        ``sum_over_cat_ice`` sums categories, ``select_cat_ice`` selects one.
    ice_cat
        One-based category index used with ``select_cat_ice``.
    unit
        Unit string copied to output NetCDF attributes.
    transform
        ``0`` none, ``1`` log10, ``2`` natural log.
    trafo_shift
        Shift used by log transforms: ``log(value + trafo_shift)``.
    limit
        ``0`` none, ``1`` minimum only, ``2`` maximum only, ``3`` both.
    max_limit, min_limit
        Bounds used when ``limit`` enables clipping.
    """
    ndims: int = 0
    dim: int = 0
    off: int = 0
    variable: str = ""
    name_incr: str = ""
    name_rest_n: str = ""
    name_bkg_din: str = ""
    rst_file: str = ""
    k_name: str = "lev"
    operation: str = ""
    ice_cat: int = 0
    unit: str = ""
    transform: int = 0
    trafo_shift: float = 0.0
    limit: int = 0
    max_limit: float = 0.0
    min_limit: float = 0.0

    @classmethod
    def from_config(cls, section: configparser.SectionProxy) -> "StateField":
        """Build a StateField from a configparser section, using defaults for absent keys."""
        _getters: dict[type, typing.Callable] = {int: section.getint, float:
                                                 section.getfloat, str: section.get}
        kwargs = {}
        for f in dataclasses.fields(cls):
            if f.name not in section:
                continue
            getter = _getters[typing.cast(type, f.type)]
            value = getter(f.name)
            if f.type is str:
                value = value.strip("'\"")
            kwargs[f.name] = value

        return cls(**kwargs)


# @dataclass
# class LocalStateField:
#     dim: int = 0
#     off: int = 0


class Sfields:
    """Define model fields and offsets for the PDAF state vector.

    ``[statevector].svnames`` controls which sections are loaded. Each listed
    section is converted to a :class:`StateField`, assigned a wet-grid length,
    and placed into a contiguous PE-local state vector. This is the main place
    to extend the assimilation state for salinity, sea-ice, biogeochemistry, or
    variables from another model.
    """

    def __init__(self, config: configparser.ConfigParser, wet_grid: model.WetGrid,
                 pe:parallel.Parallel) -> None:
        varnames = [v.strip() for v in config.get("statevector", "svnames").split(",")]
        self.fields = self.setup_statevector(config, wet_grid, varnames)
        self.dim_state_p = sum(self.fields[vname].dim for vname in self.fields)
        self.dim_state = pe.comm_filter.reduce(self.dim_state_p, op=mpi4py.MPI.SUM, root=0)
        self._log_setup(pe.mype_filter)

    def setup_statevector(self, config: configparser.ConfigParser, wet_grid: model.WetGrid, varnames: list):
        """Validate fields, assign dimensions/offsets, and compute total size."""
        fields = {}
        offset = 0
        for vname in varnames:
            fields[vname] = StateField.from_config(config[vname])
            if fields[vname].ndims == 2:
                fields[vname].dim = int(wet_grid.sdim2d)
            elif fields[vname].ndims == 3:
                fields[vname].dim = int(wet_grid.sdim3d)
            else:
                raise ValueError(f"NEMO-PDAF: cannot handle {fields[vname].ndims} number of dimensions.")
            fields[vname].off = offset
            offset += fields[vname].dim
        return fields

    def _log_setup(self, rank):
        log.info("*** Setup of state vector ***")
        log.info(f" --- Number of fields in state vector: {len(self.fields):5d}")
        log.info("NEMO-PDAF    pe   ID   variable    ndims       dim     offset")
        if rank == 0:
            screen = 0
        else:
            screen = 3
        for index, vname in enumerate(self.fields, start=1):
            field = self.fields[vname]
            log.info(
                f"{rank:4d}{index:5d}   {field.variable:10s}"
                f"{field.ndims:5d}{field.dim:10d}{field.off:10d}",
            screen=screen, ranks=rank)
        if rank == 0:
            screen = 0
        else:
            screen = 2
        log.info(
            f"NEMO-PDAF  PE {rank:4d}  PE-local full state dimension: {self.dim_state_p:10d}",
        screen=screen, ranks=rank)
        if rank == 0:
            log.info(f"NEMO-PDAF  Global state dimension: {self.dim_state:10d}")
