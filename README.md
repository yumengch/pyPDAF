# pyPDAF
A Python interface to the Fortran-written data assimilation library - [PDAF](http://pdaf.awi.de/trac/wiki)

![GitHub Workflow Status](https://github.com/yumengch/pyPDAF/actions/workflows/test_build.yaml/badge.svg)


## Prerequisite:
- `PDAF-V2.1`
- `Fortran compiler: e.g.:gfortran/intel fortran`
- `a message passing interface (MPI) implementation: e.g. openMPI/MPICH`
- `Python>=3.8`


## Installation:
- Currently, Fortran-written PDAF is compiled together with pyPDAF. Hence, the Fortran compiler options need to be specified in the PDAF section of [`setup.cfg`](setup.cfg).
- Install Python package: ```pip install .```

## Run example:
```bash
mpiexec -n 4 python -u example/main.py
```
Will run the model with 4 ensemble members where each member uses 1 process. 

## Documentation:
Currently, it interfaces with subroutines of ```PDAF-V2.1``` with an example for online coupling with PDAF using a simple model based on the [tutorial](http://pdaf.awi.de/trac/wiki/FirstSteps) from PDAF. Some interface in Python changes slightly due to different ways to handling return values in Python from Fortran. 

[This documentation](https://yumengch.github.io/pyPDAF/index.html) contains more details on installation, Python interface, and implementation guide on pyPDAF. The simplified interfaces in PDAF are not supported. 

## Contributors:
Yumeng Chen, Lars Nerger

pyPDAF is mainly developed and maintainde by National Centre for Earth Observation and University of Reading.

<img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/UoR-logo.png?raw=true" height="50" /> <img src="https://github.com/nansencenter/DAPPER/blob/master/docs/imgs/nceologo1000.png?raw=true" height="50">
