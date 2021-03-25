cmake_minimum_required(VERSION 3.11)
project(mgmc CXX)

option(MGMC_USE_OMP "Compile MGMC with OpenMP for shared memory parallelism" ON)
option(MGMC_GO_FAST "Compile MGMC with flags which are not for the faint of heart" OFF)
option(MGMC_TESTS "Compile MGMC unit tests" OFF)

find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND)
  message(STATUS "Could not find a local install of yaml-cpp")
  message(STATUS "Will download and build yaml-cpp instead")

  # YAML-CPP has horrible default cmake options. Need to change them
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
  set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")

  include(FetchContent)

  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
    GIT_TAG        yaml-cpp-0.6.3
  )

  FetchContent_MakeAvailable(yaml-cpp)
else()
  message(STATUS "Using local install of yaml-cpp")
endif()

include(sourcelist.cmake) # For MGMC_SOURCE_FILES list

add_executable(mgmc src/main.cpp vendor/docopt/docopt.cpp ${MGMC_SOURCE_FILES})
target_include_directories(mgmc PRIVATE include vendor vendor/ndarray/include)
target_compile_features(mgmc PRIVATE cxx_std_17)
target_compile_options(mgmc PRIVATE -W -Wall -Wextra -Wpedantic)
target_compile_options(mgmc PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(mgmc PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_compile_options(mgmc PRIVATE $<$<BOOL:MGMC_GO_FAST>:-Ofast -ffast-math>)
target_link_libraries(mgmc PUBLIC yaml-cpp)

# Find OpenMP if desired
if(MGMC_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(mgmc PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()
