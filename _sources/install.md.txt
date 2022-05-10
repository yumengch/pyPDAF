# Installation

## Prerequisite:
- `PDAF-V2.0`
- `Fortran compiler: e.g.:gfortran/intel fortran`
- `a message passing interface (MPI) implementation: e.g. openMPI/MPICH`
- `Python>=3.8`

---
**NOTE**
- As a parallel framwork, PDAF makes use of MPI. The MPI is used by pyPDAF as well. To avoid unexpected errors, we recommend using the same MPI implementation (OpenMP/MPICH) for both the PDAF and pyPDAF.
- pyPDAF uses [MPI4py](https://mpi4py.readthedocs.io/en/stable/). To specify the MPI implementation for MPI4py, the following method can be used:
```bash
export CC=/path/to/mpicc python
env MPICC=/path/to/mpicc python -m pip install mpi4py
```
---

## Install pyPDAF:
- ```pip install -e .```
