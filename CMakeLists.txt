cmake_minimum_required(VERSION 3.15.0)

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

if (NOT DEFINED CEREALIZABLE)
	option(CEREALIZABLE "Enable serialisation of various classes" ON)
endif()
if (CEREALIZABLE)
	add_definitions(-DCERALIZABLE)
	add_compile_definitions(CEREALIZABLE)
endif()

# start project
if (NOT DEFINED VERSION_NR)
project(pylimer_tools VERSION "0.2.4" LANGUAGES CXX)
else()
project(pylimer_tools VERSION ${VERSION_NR} LANGUAGES CXX)
endif()
add_definitions( -DOVERALL_PROJECT_VERSION="${CMAKE_PROJECT_VERSION}")

include(${CMAKE_CURRENT_LIST_DIR}/vendor/myConfigureFile.cmake)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_CURRENT_LIST_DIR}/vendor/CMakeLists.txt)

# prepare python bindings
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/vendor/pybind11)

message("Using toolchain file ${CMAKE_TOOLCHAIN_FILE}")

# add sources & targets
add_subdirectory(src/pylimer_tools_cpp)
add_dependencies(pylimer_tools igraph)
add_dependencies(pylimer_tools nlopt)
add_dependencies(pylimer_tools cerealLib)
add_dependencies(pylimer_tools Spectra)

pybind11_add_module(pylimer_tools_cpp MODULE 
  src/pylimer_tools_cpp/pybind11/readers.cpp
  src/pylimer_tools_cpp/pybind11/entities.cpp
  src/pylimer_tools_cpp/pybind11/writers.cpp
  src/pylimer_tools_cpp/pybind11/generators.cpp
  src/pylimer_tools_cpp/pybind11/topology.cpp
  src/pylimer_tools_cpp/pybind11/calculations.cpp
  src/pylimer_tools_cpp/pybind11/simulations.cpp
  src/pylimer_tools_cpp/pybind11/pylimer_tools_cpp.cpp
)

# add other dependencies
include_directories(${igraph_INCLUDE_DIRS})
include_directories(${nlopt_INCLUDE_DIRS})
include_directories(${cereal_INCLUDE_DIRS})
include_directories(${Spectra_INCLUDE_DIRS})

message("pylimer_tools_cpp includes igraph from ${igraph_INCLUDE_DIRS}, nlopt from ${nlopt_INCLUDE_DIRS}, Spectra from ${Spectra_INCLUDE_DIRS} and cereal from ${cereal_INCLUDE_DIRS}")

# link
target_link_libraries(pylimer_tools_cpp 
  PUBLIC 
    pylimer_tools
)

if(CMAKE_COMPILER_IS_GNUCXX)
target_link_libraries(pylimer_tools_cpp 
  PUBLIC 
    stdc++fs
)
endif()

# compile & install
target_compile_definitions(pylimer_tools_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS pylimer_tools_cpp LIBRARY DESTINATION .)
