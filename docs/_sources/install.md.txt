
# Installation

There are two ways of installing pyPDAF.

## Conda
The easiest approach is using `conda`. Currently, `pyPDAF` is available from
`conda` for `Windows`, `Linux` and `MacOS`. The installation can be obtained via:
```bash
conda create -n pypdaf -c conda-forge yumengch::pypdaf
```
After installation, `pyPDAF` can be used by activating the conda environment
`conda activate pypdaf`.

## Source code
In some cases, it is desirable to compile pyPDAF so that one can use their
favourite compiler as well as MPI and BLAS implementation.

In this case, pyPDAF source code can be obtained from source
```bash
git clone --recurse-submodules https://github.com/yumengch/pyPDAF.git
cd pyPDAF
```

The package can be installed with:
```bash
python -m pip install . -v \
    --config-settings=setup-args="-Dblas_lib=[LIBS]" \
    --config-settings=setup-args="-Dincdirs=[INCDIRS]" \
    --config-settings=setup-args="-Dlibdirs=[LIBDIRS]" \
    --config-settings=setup-args="-Dmpi_mod=MPIF90"^
    --config-settings=setup-args="-Dbuildtype=release"
```
Here, `LIBS`, `INCDIRS`, and `LIBDIRS` are elements of a list, separated by
`,` similar to a Python list.
   - `LIBS` are all the required library names for BLAS libraries.
   - `LIBDIRS` are directories of these libraries
   - `INCDIRS` are include directories of libraries
   - `MPIF90` is the path to `mpi.f90`. This is useful for the case where
     `mpi.mod` is not directly provided but requires compiling by the user. This
     is optional.

One can adjust the compiler, compiler and linker flags by changing environment
variables such as `CC`, `FC`, `CFLAGS` and `FCFLAGS`. See [flags](https://mesonbuild.com/Reference-tables.html#compiler-and-linker-flag-environment-variables) and [compiler](https://mesonbuild.com/Reference-tables.html#compiler-and-linker-selection-variables) table for references. One could also directly
modify [`meson.build`](meson.build) which might require more knowledge of meson.

An example of installing pyPDAF in Linux or Mac:
```bash
CC=mpicc FC=mpifort python -m pip install . -v -Cbuild-dir=build \
    --config-settings=setup-args="-Dblas_lib=['openblas']" \
    --config-settings=setup-args="-Dincdirs=['/usr/lib/']" \
    --config-settings=setup-args="-Dlibdirs=['/usr/include']" \
    --config-settings=setup-args="-Dbuildtype=release"
```
In Windows, one can use
```cmd
set CXX=clang-cl
set CC=clang-cl
set FC=flang-new
set MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include
set MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64

python -m pip install . -v ^
    -Cbuild-dir=build --config-settings=setup-args="-Dblas_lib=openblas"^
    --config-settings=setup-args="-Dincdirs="C:\Program Files (x86)\blas\include"^
    --config-settings=setup-args="-Dlibdirs="C:\Program Files (x86)\blas\lib"^
    --config-settings=setup-args="-Dmpi_mod=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.f90"^
    --config-settings=setup-args="-Dbuildtype=release"
```
where the variable `MSMPI_INC` and `MSMPI_LIB64` are required environment
variable for using `MSMPI`.

Please [raise an issue](https://github.com/yumengch/pyPDAF/issues/new) if you
have any questions or problems with this.
