# Set the PDAF library name, it can be pdaf-var or pdaf-d
set(PDAF_NAME "pdaf-var")

# set MPI library information
set(MPI_Fortran_INCLUDE_PATH
    "C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
    CACHE STRING "path to the MPI Fortran mpi.h include directory")
set(MPI_Fortran_MODULE_INCLUDE_PATH
    "C:/Program Files (x86)/Microsoft SDKs/MPI/Include/x64"
    CACHE STRING "path to the include directory to compile MPI Fortran module")
set(MPI_Fortran_MODULE_DIR
    "C:/Users/cymji/Desktop/pyPDAF/PDAF_V2.2.1/mpi_include" 
    CACHE STRING "path to the installed module directory of MPI_Fortran")
# We have to compile MPI-MPI mod file ourselves in Windows
set(MPI_Fortran_MODULE_SRC_FILE "C:/Program Files (x86)/Microsoft SDKs/MPI/Include/mpi.f90"
    CACHE STRING "path to the MPI Fortran module source file" )
# This overrides any language flags and must be given before add_library
set_source_files_properties(${MPI_Fortran_MODULE_SRC_FILE} PROPERTIES COMPILE_FLAGS "/O2 /Qdiag-disable:10448 /Qdiag-disable:10423")

# ideally, we should provide compiler options based on target or source files, 
# but this will force us to use multiple config files or encode it in the CMakelist.txt file
# set(PDAF_FLAGS_RELEASE -O3 -ffree-line-length-none -fdefault-real-8 -fPIC -fopenmp)
set(PDAF_FLAGS_RELEASE /O2 /4R8 /Qdiag-disable:10448 /Qdiag-disable:10423)
set(LEGACY_FLAGS_RELEASE /O2 /Qdiag-disable:10448 /Qdiag-disable:10423)
# Set compiler flags for Debug configurations
set(PDAF_FLAGS_DEBUG /Od /debug:full /traceback /check:all /Qtrapuv /4R8 /Qdiag-disable:10448 /Qdiag-disable:10423)
set(LEGACY_FLAGS_DEBUG /Od /debug:full /traceback /check:all /Qtrapuv /Qdiag-disable:10448 /Qdiag-disable:10423)
