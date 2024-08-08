# Installation

There are two ways of installing pyPDAF. 
- The easiest approach is using `conda`. Currently, `pyPDAF` is available from `conda` for `Windows`, `Linux` and `MacOS`. The installation can be obtained via:
```bash
conda create -n pyPDAF -c yumengch -c conda-forge pyPDAF
```
You can start to use `pyPDAF` by `conda activate pyPDAF`.
- In HPC or cluster environment, it might not be desirable to use compilers and MPI implementation provided by conda. In this case, pyPDAF can be installed from source
```bash
git clone https://github.com/yumengch/pyPDAF.git
cd pyPDAF
git submodule update --init --recursive
pip install -v .
```
The `pip` command compiles both `PDAF V2.1` and its C interface. To customise the compiler options with the local machine, it is necessary to specify the compiler, compiler options, path to the dependent libraries. In our case, the dependent library is `BLAS`, `LAPACK`, and `MPI` implementation. 
   - The installation requires `Cython`, `mpi4py`, and `numpy` package.
   - The Fortran compiler options need to be specified in the PDAF section of [`setup.cfg`](setup.cfg). Note that the `-fPIC` compiler option is required to create a Python package. Note that these are only relevant on non-Windows machines. For Windows machines, `MSVC` and `Intel Fortran compilers` are used by default and adaptations for other compilers will need changes in `CMakeLists.txt` in [PDAFBuild/CMakeLists.txt](PDAFBuild/CMakeLists.txt) and [pyPDAF/fortran/CMakeLists.txt](pyPDAF/fortran/CMakeLists.txt).
   - Options in pyPDAF section of `setup.cfg` requires the following options:
      - `pwd` is the absolute path to the pyPDAF repository directory
      - `CC` is the C compiler used by Cython, e.g. `CC=mpicc` for GNU compiler or `CC=mpiicc` for Intel compiler. This option is not usable in Windows as only `MSVC` is supported.
      - `condaBuild` -- ignore this option as is only relevant for `conda build` scenario
      - `useMKL` decides if you use Intel's Math Kernel Library (MKL). If `True` is given, `MKLROOT` must be specified which is the absolute path to the static MKL library
      - `LAPACK_PATH` and `LAPACK_LIBRARY` is the path to the BLAS and LAPACK directory and the linking flag respectively. They can be delimited by `,`. For example, we can have `LAPACK_LIBRARY=blas,lapack`. Do not give `-lblas` as `setuptools` deal with the format to the linker.
      - `MPI_INC_PATH`, `MPI_MOD_PATH`, and `MPI_LIB_PATH` are only relevant in Windows, which is the path to `.h` file, `.f90` file, and `.lib` file respectively. These paths are usually `C:\Program Files (x86)\Microsoft SDKs\MPI\Include\x64`, `C:\Program Files (x86)\Microsoft SDKs\MPI\Include`, and `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64` respectively.
