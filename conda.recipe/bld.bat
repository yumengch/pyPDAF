set CXX=clang-cl
set CC=clang-cl
set MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include
set MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64

"%PYTHON%" -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv^
 -Cbuild-dir=build --config-settings=setup-args="-Dblas_lib=openblas"^
 --config-settings=setup-args="-Dincdirs=C:\Users\ia923171\AppData\Local\anaconda3\envs\devel\Library\include"^
 --config-settings=setup-args="-Dlibdirs=C:\Users\ia923171\AppData\Local\anaconda3\envs\devel\Library\lib"^
 --config-settings=setup-args="-Dmpi_mod=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.f90"^
 --config-settings=setup-args="-Dbuildtype=release"