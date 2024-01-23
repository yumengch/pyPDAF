# Installation

## Prerequisite:
- `Fortran compiler: e.g.:gfortran/intel fortran`
- `a message passing interface (MPI) implementation: e.g. openMPI/MPICH/MS-MPI`
- `BLAS and LAPACK installation or Intel MKL library compatiable with the Fortran compiler`
- `Python>=3.8`

---
**NOTE**
- pyPDAF uses [MPI4py](https://mpi4py.readthedocs.io/en/stable/). The MPI4py and the compile-time MPI should use the same MPI implementations to avoid any issues. To specify the MPI implementation for MPI4py, the following method can be used:
```bash
export CC=/path/to/mpicc python
env MPICC=/path/to/mpicc python -m pip install mpi4py
```
---

## Install pyPDAF:
- First, provide path to the compiler and libraries in `setup.cfg`
- ```pip install -e .```
