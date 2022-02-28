# pyPDAF
A (incomplete) Python interface to the Fortran-written data assimilation library - [PDAF](http://pdaf.awi.de/trac/wiki)

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/yumengch/pyPDAF/test_build)


## Prerequisite:
- `PDAF-V1.16`
- `Fortran compiler: e.g.:gfortran/intel fortran`
- `a message passing interface (MPI) implementation: e.g. openMPI/MPICH`
- `Python>=3.8`


## Installation:
- Currently, Fortran-written PDAF is compiled together with pyPDAF. Hence, the Fortran compiler options including needs to be specified in [`setup.cfg`](setup.cfg).
- Install Python package: ```pip install -e .```

## Run example:
```bash
mpiexec -n 8 python -u example/main.py```
Will run the model with 4 ensemble members where each member uses 2 processes. 

## Note:
Currently, it only interfaces with limited subroutines of ```PDAF-V1.16``` with an example for online coupling with PDAF using a simple model based on the [tutorial](http://pdaf.awi.de/trac/wiki/FirstSteps) from PDAF. Some interface in Python changes slightly due to different ways to handling return values in Python from Fortran. It is possible to check Python interface [here](https://yumengch.github.io/pyPDAF/index.html). More subroutines will be supported in future release. 

## Development
The interface between Fortran and Python relies on [Cython](https://cython.readthedocs.io/en/stable/index.html) and the feature of [iso_c_binding](https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html) in Fortran. 


## Contributors:
Yumeng Chen, Lars Nerger

pyPDAF is mainly developed and maintainde by National Centre for Earth Observation and University of Reading.

<img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/UoR-logo.png?raw=true" height="50" /> <img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/nceologo1000.png?raw=true" height="50">
