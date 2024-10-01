
# Installation

There are two ways of installing pyPDAF. 
- The easiest approach is using `conda`. Currently, `pyPDAF` is available from `conda` for `Windows`, `Linux` and `MacOS`. The installation can be obtained via:
```bash
conda create -n pyPDAF -c yumengch -c conda-forge pyPDAF
```
After installation, `pyPDAF` can be used by activating the conda environment `conda activate pyPDAF`.

- In HPC or cluster environment, it might not be desirable to use compilers and MPI implementation provided by conda. In this case, pyPDAF can be installed from source
```bash
git clone https://github.com/yumengch/pyPDAF.git
cd pyPDAF
git submodule update --init --recursive
pip install -v .
```
The `pip` command compiles both `PDAF` and its C interface. To customise the compiler options with the local machine, it is necessary to specify the compiler, compiler options, and path to the dependent libraries. There are a few components for compilation:
1. The PDAF and PDAFc compilation:
   In pyPDAF, PDAF is compiled with `CMake`. It utilises the `.cmake` files in `PDAFBuild` directory in which examples of `.cmake` files are included for linux, mac and windows systems. It might be necessary to adapt the `CMAKE_Fortran_COMPILER`, `MPI_Fortran_INCLUDE_PATH` and `MPI_Fortran_MODULE_DIR`entries based on your local system. If the  `CMAKE_Fortran_COMPILER` is a MPI wrapper such as `mpif90` or `ftn` in cray, the `MPI_Fortran_xxx` entries can be left empty.
2. The linking of pyPDAF and PDAFc:
    The linking is handled by `setuptools` of Python and Cython extension modules. The `setuptools` is configured by `setup.cfg`.  Here is an entry by entry explanation:
   - `pwd` is the root path of the repository directory. This should be an absolute path starting from `/`, e.g.`/home/user/pyPDAF`
   - `PDAF_dir` is the relative path to the PDAF directory. By default it is `PDAF_V2.2.1`
   - `cmake_config_path` is the path to the `.cmake` file for PDAF compilation, e.g. `/home/user/pyPDAF/PDAFBuild/linux_gfortran_openmpi_pypdaf.cmake`
   - `c_compiler` and `fortran_compiler` specifies the compiler being used. `c_compiler` can be: `gcc`, `msvc`, `icc`, `clang` and `fortran_compiler` can be: `gfortran` and `ifort`. The `fortran_compiler` should be consistent with the compiler used in cmake configuration file.
   - `condaBuild` -- a switch for building a conda package. Usually this can be ignored
   - `useMKL` decides if you use Intel's Math Kernel Library (MKL). If `True` is given, `MKLROOT` must be specified which is the absolute path to the static MKL library
   - `LAPACK_PATH` and `LAPACK_Flag` is the path to the BLAS and LAPACK directory and the linking flag respectively. They can be delimited by `,`. For example, we can have `LAPACK_Flag=blas,lapack`. Do not give `-lblas` as `setuptools` deal with the format to the linker.
   - `MPI_LIB_PATH` is the path to the MPI library. In windows, this path is usually  `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64`.