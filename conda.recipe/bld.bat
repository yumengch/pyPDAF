@echo on

set CXX=clang-cl
set CC=clang-cl
set FC=flang-new
set MSMPI_INC=%LIBRARY_INC%
set MSMPI_LIB64=%LIBRARY_LIB%
set OPENMP_ARGS=
if /I "%openmp%"=="true" set OPENMP_ARGS=--config-settings=setup-args="-Dopenmp=true"
if /I "%OPENMP%"=="true" set OPENMP_ARGS=--config-settings=setup-args="-Dopenmp=true"

"%PYTHON%" -m pip install . -v --no-build-isolation ^
    -Cbuild-dir=build ^
    --config-settings=setup-args="-Dlink_args_blas=-l"%LIBRARY_LIB%"\openblas.lib,-l"%LIBRARY_LIB%"\FortranRuntime.lib,-l"%LIBRARY_LIB%"\FortranDecimal.lib" ^
    --config-settings=setup-args="-Dincdir_blas="%LIBRARY_INC%","%BUILD_PREFIX%"\Library\include" ^
    --config-settings=setup-args="-Dmpi_mod=%LIBRARY_INC%\mpi.f90" ^
    %OPENMP_ARGS% ^
    --config-settings=setup-args="-Dbuildtype=release"
