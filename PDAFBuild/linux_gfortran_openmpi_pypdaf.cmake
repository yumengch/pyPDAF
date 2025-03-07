# Set the PDAF library name, it can be pdaf-var or pdaf-d
set(PDAF_NAME "pdaf-var")

# Set preprocessor and linker for current system
# AR is used to generate static libraries
set(CMAKE_AR "ar" CACHE FILEPATH "Archiver")
# RANKLIB is used to deal with static libraries
set(CMAKE_RANLIB "ranlib" CACHE FILEPATH "Ranlib")
# Here a linker is specified it can be the compiler used
set(CMAKE_LINKER "mpif90" CACHE FILEPATH "Linker")
# CPP is the C preprocessor
set(CMAKE_CPP "cpp" CACHE FILEPATH "C Preprocessor")

# set compiler executable
set(CMAKE_Fortran_COMPILER "mpif90")
# Set compiler flags for Release/Production configurations
# set(PDAF_FLAGS_RELEASE -O3 -ffree-line-length-none -fdefault-real-8 -fPIC -fopenmp)
set(PDAF_FLAGS_RELEASE -O3 -ffree-line-length-none -fdefault-real-8 -fPIC)
set(LEGACY_FLAGS_RELEASE -O3 -fPIC)
# Set compiler flags for Debug configurations
set(PDAF_FLAGS_DEBUG -O0 -Wall -Wextra -g -pedantic -fcheck=all -fbacktrace -ffree-line-length-none -fdefault-real-8 -fPIC)
set(LEGACY_FLAGS_DEBUG -O0 -Wall -Wextra -g -fPIC)

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