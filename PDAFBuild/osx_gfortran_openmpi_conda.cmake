# Set the PDAF library name, it can be pdaf-var or pdaf-d
set(PDAF_NAME "pdaf-var")

# set compiler executable
set(CMAKE_Fortran_COMPILER "mpif90")
# Set compiler flags for Release/Production configurations
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -ffree-line-length-none -fdefault-real-8 -fPIC  -mmacosx-version-min=11.0")
# Set compiler flags for Debug configurations
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Wall -Wextra -g -pedantic -fcheck=all -fbacktrace -ffree-line-length-none -fdefault-real-8 -fPIC")

# set MPI library information
set(MPI_Fortran_INCLUDE_PATH
    ""
    CACHE STRING "path to the include directory of MPI_Fortran")
set(MPI_Fortran_MODULE_DIR
    "" 
    CACHE STRING "path to the module directory of MPI_Fortran")

# Check if the build type is not explicitly set and default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()