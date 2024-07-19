cmake_minimum_required(VERSION 3.12)
project(PDAF VERSION 2.2.1 LANGUAGES Fortran)

if(NOT PDAF_PATH)
    message(FATAL_ERROR "PDAF_PATH is not given. Providing a PDAF_PATH in commandline by specifying -DPDAF_PATH=PATH_TO_PDAF_DIRECTORY")
endif()

add_subdirectory(${PDAF_PATH}/src ${PDAF_PATH}/src/build)

# List all source files
set(SRC_FILES
    U_PDAF_interface_c_binding.F90
    PDAF_c_binding.F90
    PDAFomi_obs_c_binding.F90
)

# check if Config Path is provided
message(STATUS ${Config_PATH})
if(NOT Config_PATH)
    message(FATAL_ERROR "Config_PATH is not given. Providing a Config_PATH in commandline by specifying -DConfig_PATH=PATH_TO_CONFIG_FILE")
endif()

# import configuration file
include(${Config_PATH})

# Check if the build type is not explicitly set and default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add the library target
add_library(PDAFc STATIC ${SRC_FILES})

# Add include directory for PDAFc
target_link_libraries(PDAFc PRIVATE ${PDAF_NAME})
target_include_directories(PDAFc PUBLIC ${PDAF_PATH}/include)
target_link_directories(PDAFc PUBLIC ${PDAF_PATH}/lib)

# Set MPI include directory
target_include_directories(PDAFc PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_include_directories(PDAFc PUBLIC ${MPI_Fortran_MODULE_DIR})
# Set PDAF include directory
target_include_directories(PDAFc PUBLIC ${PDAF_PATH}/include)

install(TARGETS PDAFc
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../lib
    ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../lib
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin  # For Windows DLLs
)