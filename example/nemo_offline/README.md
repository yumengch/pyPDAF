# pyPDAF Offline NEMO Adapter

This directory contains an offline pyPDAF workflow for assimilating
observations into a NEMO ensemble. The current example reads NEMO restart files,
packs configured wet-grid variables into a PDAF state vector, creates synthetic
SST-like observations, runs a PDAF analysis, and writes NEMO ASM background and
increment files for each ensemble member.

The code is intentionally organised as an adapter. NEMO-specific file formats
and grid conventions live in `io_nemo.py` and `model.py`; PDAF callback wiring
lives in `da_sys.py`, `collector.py`, `obs_factory.py`, `local.py`, and
`prepost.py`. To use the same pattern with another complex model, replace the
domain/IO/observation adapters while keeping the PDAF callback signatures.

## Requirements

- Python environment with `pyPDAF`, `mpi4py`, `numpy`, and `xarray`.
- A working MPI runtime compatible with `mpi4py` and `pyPDAF`.
- NEMO restart files split by ensemble member and rank.
- A NEMO domain configuration file containing grid coordinates and wet-mask
  inputs such as `top_level` and `bottom_level`.

If running with MPI, invoke the same environment through your MPI launcher, for
example:

```bash
mpiexec -n N_PROC python -u main.py
```

Here, `N_PROC` should be the same number as the domains/processors used by NEMO
because the program reads the same number of restart files as `N_PRC`.

## Directory Layout

- `main.py` reads `config.ini`, initialises MPI/PDAF, runs assimilation, and
  finalises PDAF.
- `config.ini` is the user-facing run configuration.
- `parallel.py` initialises PDAF/MPI communicators.
- `da_sys.py` wires together the model domain, state vector, observations,
  localisation, filter options, and PDAF callbacks.
- `io_nemo.py` reads NEMO restart/domain files and writes NEMO ASM output files.
- `model.py` stores NEMO local-domain metadata, coordinates, time values, and
  wet-grid indexing.
- `sfields.py` defines which model variables are included in the PDAF state.
- `transforms.py` converts between model fields and compact wet-point state
  vectors, and applies optional transforms/limits.
- `collector.py` loads the forecast/background ensemble from restart files.
- `obs_factory.py` registers observation modules with PDAF-OMI.
- `obs_sst.py` implements the current synthetic SST observation example.
- `local.py` supplies local-domain callbacks for local PDAF filters.
- `prepost.py` stores the background ensemble and writes analysis output.
- `log.py` provides rank-gated logging.

## Configure a NEMO Run

All normal user settings are in `config.ini`.

### PDAF and Filter

`[pdaf].dim_ens` must match the number of ensemble directories. With
`dim_ens = 2` and `ens_prefix = ens_`, the code expects `ens_1` and `ens_2`.

`[filter_options].filtertype` selects the PDAF algorithm. The sample uses
`7`, LESTKF. Common values are documented directly in `config.ini`. Local
filters call the local-domain functions in `local.py`.

`forget` and `type_forget` control PDAF forgetting/inflation. Use `forget = 1.0`
for no fixed forgetting, or tune the value for ensemble spread.

### NEMO IO

Restart input paths are built as:

```text
path_rst_root/ens_prefix[1-N]/path_rst_suffix/f_basename_rst_000x.nc
```

`000x` is the four-digit PDAF/filter rank. For example, ensemble member 1 on
rank 0 with the sample config is:

```text
./ens_1/ORCA2_00000072_restart_0000.nc
```

Output files are written as:

```text
path_asm_root/ens_prefix[1-N]/assim_background_increments_xxxx.nc
path_asm_root/ens_prefix[1-N]/assim_background_state_DI_xxxx.nc
```

These files are intended for NEMO ASM/IAU workflows.

### State Vector

`[statevector].svnames` is a comma-separated list of field sections. Each field
section defines one variable in the PDAF state vector.

For each field:

- `ndims = 2` stores wet surface points only.
- `ndims = 3` stores every wet point in the water column.
- `name_rest_n` is the variable read from NEMO restart files.
- `name_incr` is the variable written to increment files.
- `name_bkg_din` is the variable written to direct-initialisation background
  files.
- `transform` can be `0` none, `1` log10, or `2` natural log.
- `limit` can be `0` none, `1` minimum, `2` maximum, or `3` both.

To add salinity, for example, add a `[salinity]` section with the correct NEMO
variable names and set:

```ini
[statevector]
svnames = temperature, salinity
```

### Observations

`[observations].obsnames` lists observation sections. The current code
implements `sst`.

The sample SST module creates synthetic observations from the background
ensemble mean plus Gaussian noise. It is a working template for real
observation readers. To assimilate a real SST product, replace the synthetic
observation construction in `obs_sst.py:init_dim` with file input,
quality-control, interpolation indices, and observation errors.

Important SST options:

- `doassim`: `1` to assimilate, `0` to skip.
- `field_name`: state-vector field observed by SST.
- `noise_amp`: synthetic observation standard deviation.
- `cradius`: local observation cut-off radius.
- `sradius`: localisation support or e-folding radius.
- `loc_weight`: `0` constant, `1` exponential, `2` fifth-order polynomial,
  `3` or `4` PDAF regulated localisation variants.

## Adapting to Another Model

The PDAF-facing workflow needs five model-independent concepts:

1. A local domain with dimensions, rank-local coordinates, and global position.
2. An active grid cell mask equivalent to NEMO's wet grid.
3. A way to pack model fields into a compact state vector.
4. A way to unpack analysis vectors back into model restart/increment files.
5. Observation modules that provide values, errors, coordinates, and an
   observation operator.

For another model, start by replacing:

- `model.NemoDomain` with a domain object for your grid.
- `io_nemo.Reader` with input readers for your model state and grid files.
- `io_nemo.Writer` with output writers your model can consume.
- `transforms.field2state` and `transforms.state2field` if the active-cell
  indexing is not NEMO-like.
- `obs_sst.py` or add new observation modules for your observing system.

Keep the callback methods used by `pyPDAF.init`, `pyPDAF.assim_offline`, and
`pyPDAF.PDAFomi` unless you intentionally change the PDAF workflow.

## Typical Workflow

1. Prepare ensemble restart directories, one per member.
2. Confirm `config.ini` paths and `dim_ens`.
3. Choose state fields in `[statevector]`.
4. Choose a PDAF filter in `[filter_options]`.
5. Configure observations in `[observations]` and their sections.
6. Run `main.py` in the Python environment that contains pyPDAF.
7. Feed the generated ASM background/increment files back into NEMO.

## Current Limitations

- The SST observation module currently generates synthetic observations rather
  than reading an external product.
- The adapter assumes NEMO-style rank-split NetCDF restart files and wet masks.

