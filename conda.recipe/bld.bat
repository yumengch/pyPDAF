set CXX=clang-cl
set CC=clang-cl
set MSMPI_INC=%LIBRARY_INC%
set MSMPI_LIB64=%LIBRARY_LIB%

"%PYTHON%" -m pip install . --no-build-isolation -v^
 -Cbuild-dir=build --config-settings=setup-args="-Dblas_lib=openblas"^
 --config-settings=setup-args="-Dincdirs=%LIBRARY_INC%"^
 --config-settings=setup-args="-Dlibdirs=%LIBRARY_LIB%"^
 --config-settings=setup-args="-Dmpi_mod=%LIBRARY_INC%\mpi.f90"^
 --config-settings=setup-args="-Dbuildtype=release"