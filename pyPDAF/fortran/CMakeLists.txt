cmake_minimum_required(VERSION 3.12)
project(PDAF VERSION 2.1.0 LANGUAGES Fortran)

# In windows, MPI.mod must be compiled by ourselves
if(WIN32)
  add_library(mpimod STATIC ${MPI_Fortran_MODULE_DIR}/mpi.f90)
  set_source_files_properties(${MPI_Fortran_MODULE_DIR}/mpi.f90 PROPERTIES COMPILE_OPTIONS "/Qdiag-disable:10448")
  target_include_directories(mpimod PUBLIC ${MPI_Fortran_INCLUDE_PATH})
endif()

# List all source files
set(SRC_FILES
    U_PDAF_interface_c_binding.F90
    PDAF_c_binding.F90
    PDAFomi_obs_c_binding.F90
)

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3;-ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3;-N 1023")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   if(WIN32)
      set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "/O3;/Qdiag-disable:10448")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
   else()
      # set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
      set_source_files_properties (${SRC_FILES} PROPERTIES COMPILE_OPTIONS "-O3")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -mkl")
      #KB For newer versions of CMake - this is the way to go
      enable_language(C)
      set(BLA_VENDOR Intel10_64lp_seq)
      #KB https://cmake.org/cmake/help/latest/module/FindBLAS.html
  endif()
endif()

# Add the library target
add_library(pdafc STATIC ${SRC_FILES})

# Set library properties, such as include directories
target_include_directories(pdafc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pdafc PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(pdafc PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_include_directories(pdafc PUBLIC ${MPI_Fortran_MODULE_DIR})


