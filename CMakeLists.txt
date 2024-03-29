cmake_minimum_required(VERSION 3.11)
project(mgmc CXX)

option(MGMC_USE_OMP "Compile MGMC with OpenMP for shared memory parallelism" ON)
option(MGMC_USE_MPI "Compile MGMC with MPI for distributed memory parallelism" OFF)

# Get FetchContent for downloading dependencies
include(FetchContent)

#===============================================================================
# Get YAML-CPP version 0.7.0
message(STATUS "Downloading yaml-cpp 0.7.0")
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG        yaml-cpp-0.7.0
)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")
FetchContent_MakeAvailable(yaml-cpp)

#===============================================================================
# Get PapillonNDL latest development version
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/vendor/papillon-ndl")

#===============================================================================
# Get docopt
message(STATUS "Downloading docopt.cpp v0.6.3")
FetchContent_Declare(docopt.cpp
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
  GIT_TAG        v0.6.3
)
if(NOT docopt.cpp_POPULATED)
  FetchContent_Populate(docopt.cpp)
  add_subdirectory(${docopt.cpp_SOURCE_DIR} ${docopt.cpp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#===============================================================================
# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

#===============================================================================
# Get PCG-C++
message(STATUS "Downloading PCG-C++")
FetchContent_Declare(PCG_CXX
  GIT_REPOSITORY https://github.com/HunterBelanger/pcg-cpp.git
  GIT_TAG        feature/cmake
)
set(PCG_CPP_INSTALL OFF CACHE BOOL "Install PCG-CPP")
FetchContent_MakeAvailable(PCG_CXX)

#===============================================================================
# Get SOBOL
add_subdirectory("${CMAKE_SOURCE_DIR}/vendor/sobol")

#===============================================================================
include(sourcelist.cmake) # For MGMC_SOURCE_FILES list

# Get the Git hash for define
execute_process(COMMAND git describe --dirty --always
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                OUTPUT_VARIABLE MGMC_GIT_HASH_RAW
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )

add_executable(mgmc src/main.cpp ${MGMC_SOURCE_FILES})
target_include_directories(mgmc PRIVATE include)
target_compile_features(mgmc PRIVATE cxx_std_17)
target_compile_options(mgmc PRIVATE -W -Wall -Wextra -Wpedantic)
target_link_libraries(mgmc PUBLIC docopt_s yaml-cpp PapillonNDL::PapillonNDL NDArray::NDArray PCG::PCG_CXX sobol)
target_compile_definitions(mgmc PUBLIC MGMC_GIT_HASH=\"${MGMC_GIT_HASH_RAW}\")

# Find OpenMP if desired
if(MGMC_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(mgmc PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

# Find MPI if desired
if(MGMC_USE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_CXX_FOUND)
    target_link_libraries(mgmc PUBLIC MPI::MPI_CXX)
    target_compile_definitions(mgmc PUBLIC MGMC_USE_MPI)
  endif()
endif()

#===============================================================================
# Copy input files to build directory so that the user can try them right
# after the build.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input_files/c5g7.yaml"
               "${CMAKE_CURRENT_BINARY_DIR}/c5g7.yaml")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input_files/ref_sqr_c5g7.yaml"
               "${CMAKE_CURRENT_BINARY_DIR}/ref_sqr_c5g7.yaml")

# Allow for install of executable
install(TARGETS mgmc
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
