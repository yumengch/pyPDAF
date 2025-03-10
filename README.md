# pyPDAF - A Python interface to the Fortran-written data assimilation library

pyPDAF provides a Python interface to the established
[Parallel Data Assimilation Framework (PDAF)](https://pdaf.awi.de/trac/wiki).
The original framework is used with various regional and global
climate models including atmosphere, ocean, hydrology, land surface
and sea ice models. These models are typically written in Fortran
which can be easily used with PDAF. pyPDAF can become useful
in the following scenarios:
- With an increasing number of Python-coded numerical models,
  especially machine learning models, pyPDAF is a convenient tool
  to implement data assimilation (DA) systems purely in Python.
- Alternatively, pyPDAF can be used to set up offline
  data assimilation system. In such a system, the model fields in
  restart files are replaced by analyses generated by pyPDAF.
  This can be an attractive alternative to the original Fortran
  implementations considering the simplicity of code implementation
  and package management in Python.

The interface inherits the efficiency of the data assimilation
algorithms in Fortran, and the flexibility to be applied to different
models and observations. This means that users of pyPDAF can couple
the DA algorithms with any types of model and observations without
the need to coding the actual DA algorithms. This allows the users
to focus on the specific research problems. The framework includes
various ensemble DA algorithms including many variants of ensemble
Kalman filters, particle filters and other non-linear filters.
It also provides framework for variants of 3DVar. A full list of
supported methods can be found
[here](https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF)

## Getting Started
It is recommended to install pyPDAF via `conda`:
```bash
conda create -n pypdaf -c conda-forge yumengch::pypdaf==1.0.2
```
You can also install locally from the source code using `pip`
by setting up `setup.cfg` and `cmake` configurations
with examples given in [PDAFBuild](PDAFBuild).

## Building a DA system with pyPDAF
For users without prior experience with PDAF, we highly recommend to
start with the tutorial here:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yumengch/pyPDAF/).

To construct a parallel ensemble DA system,
in [example](example) directory, we provide both `online`
and `offline` examples.
pyPDAF and PDAF both utilise `Message Passing Interface (MPI)`
parallelisation. Hence, to run the example, it needs to be executed
from commandline using `mpiexec`. For example,
```bash
mpiexec -n 4 python -u example/online/main.py
```
will run the example with 4 processes.
The example is based on
the [tutorials](http://pdaf.awi.de/trac/wiki/FirstSteps) of the original PDAF.


## Documentation:
The most up-to-date pyPDAF has interface with ```PDAF-V2.3.1```.
A [documentation](https://yumengch.github.io/pyPDAF/index.html) is provided.
The interface follows the naming convention of PDAF.
One major difference is the localisation functions
in the [Observation Module Infrastructure (OMI)](https://pdaf.awi.de/trac/wiki/PDAF_OMI_Overview).
In PDAF, one can simply call `PDAFomi_init_dim_obs_l` or `PDAFomi_localize_covar`.
In pyPDAF, these subroutines are broken into three functions:
`pyPDAF.PDAF.omi_init_dim_obs_l_iso`,
`pyPDAF.PDAF.omi_init_dim_obs_l_noniso`,
`pyPDAF.PDAF.omi_init_dim_obs_l_noniso_locweights` for isotropic,
non-isotropic and horizontal and vertically separated
non-isotropic localisation. The suffix is applied similarly
to `pyPDAF.PDAF.omi_localize_covar`.
Details of the application of these localisation for
[`PDAFomi_init_dim_obs_l`](https://pdaf.awi.de/trac/wiki/OMI_observation_modules#init_dim_obs_l_OBSTYPE)
and
[`PDAFomi_localize_covar`](https://pdaf.awi.de/trac/wiki/OMI_observation_modules#localize_covar_OBSTYPE)
can be found in PDAF documentation.


## Having questions?
We welcome issues, pull requests, feature requests and any other discussions in the issues section.

## Contributors:
Yumeng Chen, Lars Nerger

pyPDAF is mainly developed and maintained by National Centre for Earth Observation and University of Reading.

<img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/UoR-logo.png?raw=true" height="50" /> <img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/nceologo1000.png?raw=true" height="50">
