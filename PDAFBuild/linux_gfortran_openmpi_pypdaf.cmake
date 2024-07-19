# Set the PDAF library name, it can be pdaf-var or pdaf-d
set(PDAF_NAME "pdaf-var")

# Set preprocessor and linker for current system
# AR is used to generate static libraries
set(CMAKE_AR "ar" CACHE FILEPATH "Archiver")
# RANKLIB is used to deal with static libraries
set(CMAKE_RANLIB "ranlib" CACHE FILEPATH "Ranlib")
# Here a linker is specified it can be the compiler used
set(CMAKE_LINKER "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openmpi-5.0.3-ummuvvrxsncu6txyqbvehgqkpijndmql/bin/mpif90" CACHE FILEPATH "Linker")
# CPP is the C preprocessor
set(CMAKE_CPP "cpp" CACHE FILEPATH "C Preprocessor")

# set compiler executable
set(CMAKE_Fortran_COMPILER "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openmpi-5.0.3-ummuvvrxsncu6txyqbvehgqkpijndmql/bin/mpif90")
# Set compiler flags for Release/Production configurations
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffree-line-length-none -fdefault-real-8 -fPIC")
# Set compiler flags for Debug configurations
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Wall -Wextra -g -pedantic -fcheck=all -fbacktrace -ffree-line-length-none -fdefault-real-8 -fPIC")

# set MPI library information
set(MPI_Fortran_INCLUDE_PATH
    "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openmpi-5.0.3-ummuvvrxsncu6txyqbvehgqkpijndmql/include/"
    CACHE STRING "path to the include directory of MPI_Fortran")
set(MPI_Fortran_MODULE_DIR
    "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openmpi-5.0.3-ummuvvrxsncu6txyqbvehgqkpijndmql/lib/" 
    CACHE STRING "path to the module directory of MPI_Fortran")

# set BLAS information
set(BLAS_NAME "openblas")
set(BLAS_INCLUDE_PATH "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openblas-0.3.26-fhoiyp4jbfxtxry3d6skn5he4hqh7ur7/include")
set(BLAS_LIB_PATH "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/openblas-0.3.26-fhoiyp4jbfxtxry3d6skn5he4hqh7ur7/lib")

# set LAPACK information
set(LAPACK_NAME "scalapack")
set(LAPACK_INCLUDE_PATH "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/netlib-scalapack-2.2.0-izrouv2qyqsuchli4s3jybeksc6m7zdo/include")
set(LAPACK_LIB_PATH "/home/username/software/spack-0.22.1/opt/spack/linux-ubuntu22.04-skylake/gcc-13.2.0/netlib-scalapack-2.2.0-izrouv2qyqsuchli4s3jybeksc6m7zdo/lib")

# set PDAF information
set(PDAF_INCLUDE_PATH "/home/username/Documents/pyPDAF/PDAF_V2.2.1/include")