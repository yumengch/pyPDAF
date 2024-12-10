
# Installation

There are two ways of installing pyPDAF.

## Conda
The easiest approach is using `conda`. Currently, `pyPDAF` is available from `conda` for `Windows`, `Linux` and `MacOS`. The installation can be obtained via:
```bash
conda create -n pypdaf -c conda-forge yumengch::pypdaf==1.0.1
```
After installation, `pyPDAF` can be used by activating the conda environment `conda activate pypdaf`.

## Source code
In some cases, it might be desirable to compile pyPDAF using your preferred compilers and environment. In this case, pyPDAF can be installed from source
```bash
git clone https://github.com/yumengch/pyPDAF.git
cd pyPDAF
git submodule update --init --recursive
pip install -v .
```
The installation congifuration is specified in `setup.cfg`. Example setup configurations are provided in [PDAFBuild directory](https://github.com/yumengch/pyPDAF/tree/main/PDAFBuild).

Each entry of the configuration file is listed here:
- `pwd`: This is the path to the current directory of the repository, e.g. `/home/users/xxxx/pyPDAF`
- `PDAF_dir`: This is the absolute, or relative (to `pwd`) path where PDAF directory is located. One do not need to change this option.
- `cmake_config_path` is the path to the `.cmake` file for PDAF compilation. This file contains CMake configurations for PDAF. There are example configurations in the [`PDAFBuild` directory](https://github.com/yumengch/pyPDAF/tree/main/PDAFBuild)
- `condaBuild` is a switch for building a conda package. If the purpose of compiling is not building a conda package, this option should always be `False`.
- `c_compiler` and `fortran_compiler` specifies the compiler being used. `c_compiler` can be: `gcc`, `msvc`, `icc`, `clang` and `fortran_compiler` can be: `gfortran` and `ifort`. The `fortran_compiler` should be consistent with the compiler used in cmake configuration file.
- `useMKL` decides if you use Intel's Math Kernel Library (MKL). If `True` is given, `MKLROOT` must be specified which is the absolute path to the static MKL library
- `LAPACK_PATH` and `LAPACK_Flag` is the path to the BLAS and LAPACK directory and the linking flag respectively. They can be delimited by `,`. For example, we can have `LAPACK_Flag=blas,lapack`. Do not give `-lblas` as `setuptools` deal with the format to the linker.
- `MPI_LIB_PATH` is the path to the MPI library. In windows, this path is usually  `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64`.