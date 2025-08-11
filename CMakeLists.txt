# ==============================================================================
# Pylimer Tools - CMake Configuration
# ==============================================================================
# This project provides C++ tools for polymer simulations with Python bindings.
# Key features:
# - High-performance C++ core library with Python bindings via pybind11
# - Dependencies: igraph, nlopt, Eigen3, cereal, Spectra
# - Supports OpenMP for parallel computations
# - Optional features: code coverage, memory leak analysis
# ==============================================================================

cmake_minimum_required(VERSION 3.15.0)

# Modern CMake policies for better behavior
if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)  # Fix CMAKE_CXX_STANDARD behavior on macOS
endif ()
if (POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)  # Only interpret if() arguments as variables or keywords when unquoted
endif ()
if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)  # find_package() uses PackageName_ROOT variables
endif ()

# Project definition
if (NOT DEFINED VERSION_NR)
    project(pylimer_tools VERSION "0.3.2" LANGUAGES CXX)
else ()
    project(pylimer_tools VERSION ${VERSION_NR} LANGUAGES CXX)
endif ()

# Global project configuration
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Disable compiler-specific extensions
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # For clang-tidy and other tools

# Set default build type if not specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# Improve build performance
if (NOT DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
    include(ProcessorCount)
    ProcessorCount(N)
    if (NOT N EQUAL 0)
        set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
        message(STATUS "Building with ${N} parallel jobs")
    endif ()
endif ()

# Project version and metadata
add_compile_definitions(OVERALL_PROJECT_VERSION="${CMAKE_PROJECT_VERSION}")

# Include build configuration and dependencies
include(${CMAKE_CURRENT_LIST_DIR}/vendor/myConfigureFile.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/vendor/CMakeLists.txt)

# ==============================================================================
# Python bindings configuration
# ==============================================================================

# Prepare pybind11 for Python bindings
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/vendor/pybind11)

# Add main library target
add_subdirectory(src/pylimer_tools_cpp)

# Add tests (only if testing is enabled)
if (BUILD_TESTING OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_subdirectory(tests)
endif ()

# ==============================================================================
# Python module definition
# ==============================================================================

# Create the Python extension module
pybind11_add_module(pylimer_tools_cpp MODULE
    src/pylimer_tools_cpp/pybind11/readers.cpp
    src/pylimer_tools_cpp/pybind11/entities.cpp
    src/pylimer_tools_cpp/pybind11/writers.cpp
    src/pylimer_tools_cpp/pybind11/topology.cpp
    src/pylimer_tools_cpp/pybind11/calculations.cpp
    src/pylimer_tools_cpp/pybind11/simulations.cpp
    src/pylimer_tools_cpp/pybind11/generators.cpp
    src/pylimer_tools_cpp/pybind11/pylimer_tools_cpp.cpp
)

# Link against the main library
target_link_libraries(pylimer_tools_cpp PRIVATE pylimer_tools)

# Ensure the Python module gets all the dependency include directories
target_include_directories(pylimer_tools_cpp
    SYSTEM PRIVATE
        ${igraph_INCLUDE_DIRS}
        ${nlopt_INCLUDE_DIRS}
        ${cereal_INCLUDE_DIRS}
        ${Spectra_INCLUDE_DIRS}
)

# Add dependencies to ensure proper build order
add_dependencies(pylimer_tools_cpp pylimer_tools)
if (TARGET igraph::igraph)
    add_dependencies(pylimer_tools_cpp igraph::igraph)
endif ()
if (TARGET nlopt)
    add_dependencies(pylimer_tools_cpp nlopt)
endif ()
if (TARGET cerealLib)
    add_dependencies(pylimer_tools_cpp cerealLib)
endif ()
if (TARGET Spectra)
    add_dependencies(pylimer_tools_cpp Spectra)
endif ()

# Configure compilation and installation
target_compile_definitions(pylimer_tools_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

# Apply compiler flags and optimizations
include(${CMAKE_CURRENT_LIST_DIR}/vendor/myCompilerDefaultFlags.cmake)
if (COMMAND OUTPUT_FLAGS)
    OUTPUT_FLAGS(pylimer_tools_cpp)
endif ()

# Platform-specific linking
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Intel")
    target_link_libraries(pylimer_tools_cpp PRIVATE stdc++fs)
endif ()

# Installation
install(TARGETS pylimer_tools_cpp LIBRARY DESTINATION .)

message(STATUS "CMake configuration completed successfully")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if (CMAKE_TOOLCHAIN_FILE)
    message(STATUS "Using toolchain file: ${CMAKE_TOOLCHAIN_FILE}")
endif ()
