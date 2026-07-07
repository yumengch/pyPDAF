# pyPDAF

pyPDAF is a Python interface to the
[Parallel Data Assimilation Framework (PDAF)](https://pdaf.awi.de/trac/wiki),
an established Fortran library for data assimilation. PDAF is used with
atmosphere, ocean, hydrology, land-surface, and sea-ice models. pyPDAF makes
the same algorithms available from Python while keeping the performance-critical
assimilation routines in compiled PDAF code.

Data assimilation combines model forecasts with observations to produce an
improved estimate of the model state, called an analysis. With pyPDAF, users
write Python callback functions that describe the model state, observations,
observation operators, and file or memory exchange. PDAF then runs the selected
filter, variational method, or hybrid method.

pyPDAF is useful when you want to:

- build an online data assimilation system around a Python model;
- add offline data assimilation to an existing model workflow through restart
  files;
- test data assimilation ideas without rewriting PDAF algorithms;
- use Python tools for observation processing, diagnostics, or machine-learning
  models while still relying on PDAF for the core algorithms.

## Getting started

The recommended installation route is `conda`:

```bash
conda create -n pypdaf -c conda-forge yumengch::pypdaf
conda activate pypdaf
```

Source builds are also supported with `pip` and `meson-python`. See the
[installation guide](https://yumengch.github.io/pyPDAF/install.html) for
compiler, MPI, BLAS, and platform-specific notes.

## Building a DA system with pyPDAF

A minimal online workflow usually follows this pattern:

1. Configure MPI communicators with `pyPDAF.set_parallel` or
   `pyPDAF.init_parallel`.
2. Initialise PDAF with `pyPDAF.init`.
3. Initialise PDAFomi observation handling with `pyPDAF.PDAFomi.init`.
4. Start the first forecast with `pyPDAF.init_forecast`.
5. Advance the model ensemble to the next observation time.
6. Run an assimilation driver such as `pyPDAF.assimilate`.
7. Repeat the forecast-analysis cycle.
8. Finalise with `pyPDAF.deallocate`.

The user-supplied callback functions provide the model state vector,
observations, observation operator, localisation information, and analysis
distribution. For a conceptual overview, see the
[data assimilation workflow guide](https://yumengch.github.io/pyPDAF/workflow.html).

For users without prior experience with PDAF, we recommend starting with the
tutorial notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yumengch/pyPDAF/).

The [example](example) directory contains online and offline examples. It also
includes an [offline NEMO adapter](example/nemo_offline/README.md) that shows how
to couple pyPDAF to NEMO restart files, build a PDAF state vector from model
variables, assimilate SST-like observations, and write NEMO ASM background and
increment files.

pyPDAF and PDAF use the Message Passing Interface (MPI), so parallel examples
should be started from the command line with `mpiexec`. For example:

```bash
cd example
mpiexec -n 4 python -u online/main.py
```

This runs the online example with four processes.

## Inspecting PDAF options

pyPDAF exposes a few PDAF discovery functions that are useful when choosing a
filter or checking the linked PDAF library:

```python
import pyPDAF

pyPDAF.print_version()
pyPDAF.print_filter_types(verbose=3)
pyPDAF.options_filters(filtertype)
```

`print_version()` reports the PDAF version linked into pyPDAF.
`print_filter_types(verbose=0)` prints the available PDAF filter type numbers.
After choosing a filter type, `options_filters(filtertype)` prints the valid
subtypes and option information for that filter.

The same information is also available from command-line entry points after
installation:

```bash
pypdaf-version
pypdaf-filter-types
pypdaf-filter-options filtertype
```

## Documentation

The current pyPDAF interface supports PDAF-V3.0. Documentation is available at
[yumengch.github.io/pyPDAF](https://yumengch.github.io/pyPDAF/index.html).

Important starting points:

- [Data assimilation workflow](https://yumengch.github.io/pyPDAF/workflow.html)
  explains the forecast-analysis cycle, online/offline modes, and user callback
  responsibilities.
- [API reference](https://yumengch.github.io/pyPDAF/API.html) lists the
  functions exported by `pyPDAF`, `PDAF`, `PDAF3`, `PDAFomi`, `PDAFlocal`, and
  `PDAFlocalomi`.
- [User-supplied functions](https://yumengch.github.io/pyPDAF/user_functions.html)
  describes callback interfaces used by PDAF assimilation drivers.
- [Offline NEMO adapter](example/nemo_offline/README.md) demonstrates the
  pattern for coupling pyPDAF to a rank-split NEMO ensemble workflow.

The main namespaces are:

- `pyPDAF`: high-level PDAF3 workflow functions.
- `pyPDAF.PDAF3`: PDAF3 initialisation and assimilation drivers.
- `pyPDAF.PDAFomi`: observation handling, observation operators, and
  localisation helpers.
- `pyPDAF.PDAF`: lower-level utilities, diagnostics, legacy functions, and
  advanced control.
- `pyPDAF.PDAFlocal`: local analysis helper functions.

Supported PDAF methods include ensemble Kalman filter variants, local ensemble
filters, particle and nonlinear filters, 3DVar, and hybrid ensemble-variational
methods. See the PDAF
[available options](https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF)
page for the underlying algorithm list.


## Having questions?

We welcome issues, pull requests, feature requests, and discussions in the
GitHub issue tracker.

## Contributors

Yumeng Chen, Lars Nerger

pyPDAF is mainly developed and maintained by the National Centre for Earth
Observation and the University of Reading.

<img src="https://github.com/nansencenter/DAPPER/blob/master/docs/images/logos/UoR-logo.png?raw=true" height="50" /> <img src="https://github.com/nansencenter/DAPPER/blob/master/docs/images/logos/nceologo1000.png?raw=true" height="50"/>
