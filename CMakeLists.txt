cmake_minimum_required(VERSION 3.9)
project(mgmc CXX)

option(MGMC_USE_OMP "Compile MGMC with OpenMP for shared memory parallelism" ON)
option(MGMC_GO_FAST "Compile MGMC with flags which are not for the faint of heart" OFF)
option(MGMC_TESTS "Compile MGMC unit tests" OFF)

find_package(yaml-cpp REQUIRED)

include(sourcelist.cmake) # For SOURCE_FILES list

add_executable(mgmc src/main.cpp vendor/docopt/docopt.cpp ${MGMC_SOURCE_FILES})
target_include_directories(mgmc PRIVATE include vendor vendor/ndarray/include ${YAML_CPP_INCLUDE_DIR})
target_compile_features(mgmc PRIVATE cxx_std_17)
target_compile_options(mgmc PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
target_compile_options(mgmc PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(mgmc PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_compile_options(mgmc PRIVATE $<$<BOOL:MGMC_GO_FAST>:-Ofast -ffast-math>)
target_link_libraries(mgmc PUBLIC ${YAML_CPP_LIBRARIES})

# Find OpenMP if desired
if(MGMC_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(mgmc PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()
