# ==============================================================================
# Pylimer Tools Core Library
# ==============================================================================
# This file defines the main C++ library with all core functionality.
# The library is linked against Python bindings and used in tests.
# ==============================================================================

cmake_minimum_required(VERSION 3.15)

# Core library definition
project(pylimer_tools_core VERSION ${CMAKE_PROJECT_VERSION} LANGUAGES CXX)

# Create the main library target
add_library(pylimer_tools 
    # Implementation files
    pylimer_tools.cpp
    entities/AtomGraphParent.cpp
    entities/Molecule.cpp
    entities/Universe.cpp
    entities/UniverseSequence.cpp
    entities/NeighbourList.cpp
    entities/EigenNeighbourList.cpp
    io/AveFileReader.cpp
    io/DataFileParser.cpp
    io/DumpFileParser.cpp
    calc/Correlator.cpp
    calc/NormalModeAnalyzer.cpp
    topo/EntanglementDetector.cpp
    sim/OutputSupportingSimulation.cpp
    sim/MEHPForceEvaluator.cpp
    sim/MEHPForceRelaxation.cpp
    sim/MEHPForceBalance.cpp
    sim/MEHPForceBalance2.cpp
    sim/DPDSimulator.cpp
    utils/MCUniverseGenerator.cpp
)

# Set target properties
set_target_properties(pylimer_tools PROPERTIES
    VERSION ${PROJECT_VERSION}
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

# ==============================================================================
# Dependencies
# ==============================================================================

# Include dependency configuration
include(${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/CMakeLists.txt)

# Link dependencies using modern CMake target-based approach
target_link_libraries(pylimer_tools 
    PUBLIC 
        Eigen3::Eigen
    PRIVATE
        ${STD_FS_LIB}
)

# Add explicit dependencies to ensure proper build order
if (TARGET igraph::igraph)
    add_dependencies(pylimer_tools igraph::igraph)
    target_link_libraries(pylimer_tools PUBLIC igraph::igraph)
else ()
    # Fallback to variable-based linking if target doesn't exist
    target_link_libraries(pylimer_tools PUBLIC ${igraph_LIBRARIES})
endif ()
if (TARGET nlopt)
    add_dependencies(pylimer_tools nlopt)
    target_link_libraries(pylimer_tools PUBLIC nlopt)
else ()
    target_link_libraries(pylimer_tools PUBLIC ${nlopt_LIBRARIES})
endif ()
if (TARGET cerealLib)
    add_dependencies(pylimer_tools cerealLib)
endif ()
if (TARGET Spectra)
    add_dependencies(pylimer_tools Spectra)
endif ()

# Include directories for dependencies
target_include_directories(pylimer_tools 
    SYSTEM PUBLIC
        ${igraph_INCLUDE_DIRS}
        ${nlopt_INCLUDE_DIRS}
        ${cereal_INCLUDE_DIRS}
        ${Spectra_INCLUDE_DIRS}
    PUBLIC
        # Make the source directory available to targets linking against this library
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# ==============================================================================
# Filesystem library detection (for C++17/20 compatibility)
# ==============================================================================

# Modern approach: Check if std::filesystem is available without extra linking
include(CheckCXXSourceCompiles)

set(TEST_FILESYSTEM_CODE "
#include <filesystem>
int main() {
    std::filesystem::path p(\"/tmp\");
    return p.string().length();
}")

# First try without any additional libraries
check_cxx_source_compiles("${TEST_FILESYSTEM_CODE}" STD_FS_NO_LIB_NEEDED)

if (NOT STD_FS_NO_LIB_NEEDED AND NOT MSVC)
    # Try with stdc++fs
    set(CMAKE_REQUIRED_LIBRARIES stdc++fs)
    check_cxx_source_compiles("${TEST_FILESYSTEM_CODE}" STD_FS_NEEDS_STDCXXFS)
    
    if (NOT STD_FS_NEEDS_STDCXXFS)
        # Try with c++fs
        set(CMAKE_REQUIRED_LIBRARIES c++fs)
        check_cxx_source_compiles("${TEST_FILESYSTEM_CODE}" STD_FS_NEEDS_CXXFS)
    endif ()
    
    # Reset CMAKE_REQUIRED_LIBRARIES
    unset(CMAKE_REQUIRED_LIBRARIES)
endif ()

# Set the appropriate filesystem library
if (STD_FS_NEEDS_STDCXXFS)
    set(STD_FS_LIB stdc++fs)
    message(STATUS "Using stdc++fs for std::filesystem")
elseif (STD_FS_NEEDS_CXXFS)
    set(STD_FS_LIB c++fs)
    message(STATUS "Using c++fs for std::filesystem")
else ()
    set(STD_FS_LIB "")
    message(STATUS "std::filesystem available without additional linking")
endif ()

# ==============================================================================
# Optional features and dependencies
# ==============================================================================

# Code coverage option
if (NOT DEFINED CODE_COVERAGE)
    option(CODE_COVERAGE "Collect coverage from test library" OFF)
endif ()

# LAPACK support
if (LAPACKE_FOUND)
    target_include_directories(pylimer_tools SYSTEM PUBLIC ${LAPACKE_LIBRARY_DIRS})
    target_link_libraries(pylimer_tools PUBLIC ${LAPACKE_LIBRARIES})
    message(STATUS "LAPACKE found and linked")
endif ()

# OpenMP support
if (WIN32)
    message(STATUS "Windows detected: OpenMP disabled")
elseif (CODE_COVERAGE)
    message(STATUS "Code coverage enabled: OpenMP disabled")
else ()
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(pylimer_tools PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(pylimer_tools PUBLIC OPENMP_FOUND=1)
        message(STATUS "OpenMP found and enabled")
    else ()
        message(STATUS "OpenMP not found")
    endif ()
endif ()

# ==============================================================================
# Build configuration
# ==============================================================================

# Development and debug features
if (CODE_COVERAGE)
    set(ENABLE_COVERAGE ON)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        message(WARNING "Code coverage enabled: switching from Release to Debug build")
        set(CMAKE_BUILD_TYPE Debug)
    endif ()
    
    find_package(codecov)
    if (codecov_FOUND)
        add_coverage(pylimer_tools)
        message(STATUS "Code coverage enabled")
    else ()
        message(STATUS "Code coverage requested but codecov not found")
    endif ()
endif ()

# Memory leak analysis
if (LEAK_ANALYSIS AND NOT WIN32)
    target_compile_options(pylimer_tools PRIVATE -fno-omit-frame-pointer -fsanitize=address)
    target_link_options(pylimer_tools PRIVATE -fno-omit-frame-pointer -fsanitize=address)
    message(STATUS "Memory leak analysis enabled")
endif ()

# Apply compiler flags
if (COMMAND OUTPUT_FLAGS)
    OUTPUT_FLAGS(pylimer_tools)
endif ()

message(STATUS "Pylimer tools core library configured")
message(STATUS "Dependencies: igraph, nlopt, Eigen3, cereal, Spectra")
if (OpenMP_CXX_FOUND AND NOT WIN32 AND NOT CODE_COVERAGE)
    message(STATUS "OpenMP: enabled")
else ()
    message(STATUS "OpenMP: disabled")
endif ()
