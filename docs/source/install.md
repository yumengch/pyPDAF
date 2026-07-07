
# Installation

There are two common ways to install pyPDAF. Most users should start with the
conda package. Build from source when you need a specific compiler, MPI
implementation, BLAS library, or development checkout.

## Conda

The easiest approach is using `conda`. Currently, `pyPDAF` is available for
Windows, Linux, and macOS. Create an environment with:

```bash
conda create -n pypdaf -c conda-forge yumengch::pypdaf
conda activate pypdaf
```

This installs pyPDAF together with its Python dependencies and the compiled PDAF
interface libraries.

## Source code

Build from source when you need control over the compiler, MPI, or BLAS setup.
The source tree includes PDAF as a submodule, so clone recursively:

```bash
git clone --recurse-submodules https://github.com/yumengch/pyPDAF.git
cd pyPDAF
```

pyPDAF uses `meson-python` as its build backend. The general installation
command is:

```bash
python -m pip install . -v \
    --config-settings=setup-args="-Dblas_lib=[LIBS]" \
    --config-settings=setup-args="-Dincdirs=[INCDIRS]" \
    --config-settings=setup-args="-Dlibdirs=[LIBDIRS]" \
    --config-settings=setup-args="-Dmpi_mod=MPIF90" \
    --config-settings=setup-args="-Dbuildtype=release"
```

The values passed to Meson are:

- `LIBS`: BLAS or LAPACK library names, for example `openblas`.
- `LIBDIRS`: directories containing these libraries.
- `INCDIRS`: include directories required by the libraries.
- `MPIF90`: optional path to `mpi.f90`, useful when the MPI installation does
  not provide a ready-to-use `mpi.mod`.

Compilers and compiler flags are selected through standard Meson environment
variables such as `CC`, `CXX`, `FC`, `CFLAGS`, `CXXFLAGS`, and `FCFLAGS`. See
the Meson reference tables for
[compiler and linker flags](https://mesonbuild.com/Reference-tables.html#compiler-and-linker-flag-environment-variables)
and
[compiler selection](https://mesonbuild.com/Reference-tables.html#compiler-and-linker-selection-variables).

### Linux and macOS example

```bash
CC=mpicc FC=mpifort python -m pip install . -v -Cbuild-dir=build \
    --config-settings=setup-args="-Dblas_lib=['openblas']" \
    --config-settings=setup-args="-Dincdirs=['/usr/include']" \
    --config-settings=setup-args="-Dlibdirs=['/usr/lib']" \
    --config-settings=setup-args="-Dbuildtype=release"
```

### Windows example

On Windows, one possible setup uses Clang/Flang and Microsoft MPI:

```console
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

The `MSMPI_INC` and `MSMPI_LIB64` environment variables are required by the
Microsoft MPI toolchain.

Please [raise an issue](https://github.com/yumengch/pyPDAF/issues/new) if you
have any questions or problems with this.
