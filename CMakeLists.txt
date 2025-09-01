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

# Function to get version from our Python script
function(get_version_from_script VERSION_VAR)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if (Python3_Interpreter_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/bin/get_version.py --ensure-file
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            OUTPUT_VARIABLE SCRIPT_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if (SCRIPT_VERSION)
            set(${VERSION_VAR} ${SCRIPT_VERSION} PARENT_SCOPE)
        else()
            set(${VERSION_VAR} "0.0.0" PARENT_SCOPE)  # Fallback
        endif()
    else()
        # Fallback: try to read from version file directly
        if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/_version.txt)
            file(READ ${CMAKE_CURRENT_LIST_DIR}/_version.txt FILE_VERSION)
            string(STRIP "${FILE_VERSION}" FILE_VERSION)
            if (FILE_VERSION)
                set(${VERSION_VAR} ${FILE_VERSION} PARENT_SCOPE)
            else()
                set(${VERSION_VAR} "0.0.0" PARENT_SCOPE)
            endif()
        else()
            set(${VERSION_VAR} "0.0.0" PARENT_SCOPE)  # Final fallback
        endif()
    endif()
endfunction()

# Project definition
if (NOT DEFINED VERSION_NR)
    get_version_from_script(AUTO_VERSION)
    project(pylimer_tools VERSION ${AUTO_VERSION} LANGUAGES CXX)
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

# Detect if building with Pyodide
if (CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR 
    CMAKE_C_COMPILER MATCHES "pyodide" OR 
    CMAKE_CXX_COMPILER MATCHES "pyodide" OR
    DEFINED ENV{PYODIDE_BUILD} OR
    PYODIDE_BUILD)
    set(BUILDING_WITH_PYODIDE ON)
    message(STATUS "Detected Pyodide build environment")
    
    # NOTE: Using pybind11's FindPython (PYBIND11_FINDPYTHON=ON) triggers a
    # call to CMake's FindPython with the Development.Module component.
    # In the Pyodide cross-compilation environment this currently fails
    # (cannot satisfy Development.Module / Python_INCLUDE_DIRS even when we
    # manually set Python_INCLUDE_DIRS) leading to configuration abort.
    # We therefore disable it and fall back to the legacy include discovery
    # mechanism, while manually providing the include directory below.
    set(PYBIND11_FINDPYTHON OFF CACHE BOOL "Disable FindPython for Pyodide cross build" FORCE)
    # Hint pybind11 that we really are cross-compiling so it avoids
    # interrogating the target interpreter.
    set(PYBIND11_USE_CROSSCOMPILING ON CACHE BOOL "Enable pybind11 cross-compiling mode" FORCE)
    # clang-scan-deps currently struggles with the emscripten wrapper driver;
    # disable CMake's C++20 module/dependency scanning to prevent fatal
    # 'stdlib.h' / 'iostream' not found errors during scanning.
    set(CMAKE_CXX_SCAN_FOR_MODULES OFF CACHE BOOL "Disable module scanning for Pyodide" FORCE)
    
    # Get Python include directory from pyodide config
    execute_process(
        COMMAND pyodide config get python_include_dir
        OUTPUT_VARIABLE PYODIDE_PYTHON_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    
    if (PYODIDE_PYTHON_INCLUDE_DIR)
        set(Python_INCLUDE_DIRS ${PYODIDE_PYTHON_INCLUDE_DIR} CACHE PATH "Python include directory from pyodide config" FORCE)
        message(STATUS "Set Python_INCLUDE_DIRS to: ${Python_INCLUDE_DIRS}")
        
        set(PYTHON_INCLUDE_DIRS ${PYODIDE_PYTHON_INCLUDE_DIR} CACHE PATH "Python include directory for pybind11" FORCE)
        
    else()
        message(WARNING "Could not get python_include_dir from pyodide config, falling back to default behavior")
    endif()
else()
    set(BUILDING_WITH_PYODIDE OFF)
    message(STATUS "Not building with Pyodide, using standard configuration")
endif()

# Prepare pybind11 for Python bindings
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/vendor/pybind11)

# Add main library target
add_subdirectory(src/pylimer_tools_cpp)

# Add tests unless disabled for Pyodide cross build (tests depend on full libc++)
if (NOT BUILDING_WITH_PYODIDE AND (BUILD_TESTING OR CMAKE_BUILD_TYPE STREQUAL "Debug"))
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

# Pyodide/Emscripten specific link flags
if (BUILDING_WITH_PYODIDE)
    set_target_properties(pylimer_tools_cpp PROPERTIES
        LINK_FLAGS
        "-flto -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=8GB -s MODULARIZE=1"
    )
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
