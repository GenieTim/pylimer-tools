# ==============================================================================
# Compiler Optimization and Warning Configuration
# ==============================================================================
# This file sets up modern, performance-oriented compiler flags for C++20
# with appropriate warnings and optimization levels per build type.
# ==============================================================================

# Ensure C++20 standard (redundant but explicit)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set Release as default build type if not specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified")
endif ()

# ==============================================================================
# Base compiler flags - applied to all build types
# ==============================================================================

# Always include debug symbols (even in Release) for better profiling/debugging
if (NOT MSVC)
    add_compile_options(-g)
endif ()

# ==============================================================================
# Build-type specific optimizations
# ==============================================================================

# Release optimizations
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

# Debug flags with comprehensive warnings
if (MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /Wall /D_DEBUG")
else ()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wextra -Wpedantic -D_DEBUG")
endif ()

# RelWithDebInfo - optimized but with debug info (good for profiling)
if (MSVC)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi /DNDEBUG")
else ()
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
endif ()

# ==============================================================================
# High-performance optimizations (optional)
# ==============================================================================

if (NOT DEFINED HIGH_PERFORMANCE)
    option(HIGH_PERFORMANCE "Enable aggressive compiler optimizations for maximum performance" OFF)
endif ()

if (HIGH_PERFORMANCE)
    if (MSVC)
        add_compile_options(/O2 /fp:fast)
    else ()
        # Enable aggressive optimizations, but be careful with -march=native in CI/packaging
        add_compile_options(-O3 -ffast-math)

        # Only use -march=native if not cross-compiling and not in CI
        if (NOT CMAKE_CROSSCOMPILING AND NOT DEFINED ENV{CI})
            add_compile_options(-march=native)
            message(STATUS "HIGH_PERFORMANCE: Enabled -march=native (detected native compilation)")
        else ()
            message(STATUS "HIGH_PERFORMANCE: Skipping -march=native (cross-compilation or CI detected)")
        endif ()
    endif ()
    message(STATUS "HIGH_PERFORMANCE: Aggressive optimizations enabled")
endif ()

# ==============================================================================
# Compiler-specific optimizations and warnings
# ==============================================================================

if (CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR
        CMAKE_C_COMPILER MATCHES "emcc" OR
        CMAKE_CXX_COMPILER MATCHES "em\\+\\+" OR
        DEFINED ENV{PYODIDE_BUILD} OR
        PYODIDE_BUILD OR
        BUILDING_WITH_PYODIDE)
    # Emscripten/WebAssembly/Pyodide-specific optimizations
    message(STATUS "Configuring Emscripten/WebAssembly optimizations")

    # Enable WebAssembly SIMD instructions for vector operations
    add_compile_options(-msimd128)
    add_compile_options(-mrelaxed-simd)

    # Aggressive optimizations for numerical computations
    add_compile_options(-ffast-math)  # Fast floating-point math

    # Conditionally disable RTTI to reduce binary size, but only if Cereal is disabled
    # since Cereal's polymorphic serialization requires RTTI (uses typeid())
    if (NOT DEFINED CEREALIZABLE)
        # Default: check if CEREALIZABLE option will be enabled (it defaults to ON)
        option(CEREALIZABLE "Enable serialisation of various classes" ON)
    endif ()

    if (NOT CEREALIZABLE)
        add_compile_options(-fno-rtti)
        message(STATUS "Emscripten: Disabled RTTI for smaller binary size (Cereal serialization disabled)")
    else ()
        message(STATUS "Emscripten: Keeping RTTI enabled for Cereal serialization")
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        # Maximum optimization for release builds
        add_compile_options(-O3)

        # Link-time optimization
        add_compile_options(-flto)
        add_link_options(-flto)

        # Disable assertions and other runtime checks for maximum performance
        add_compile_options(-DNDEBUG)
        add_link_options(-s ASSERTIONS=0)

        message(STATUS "Emscripten Release: Enabled maximum optimizations with SIMD")
    else ()
        # Debug builds: moderate optimization with debug info
        add_compile_options(-O1 -g)
        add_link_options(
            "SHELL:-s ASSERTIONS=1"
            "SHELL:-s SAFE_HEAP=1"
        )

        message(STATUS "Emscripten Debug: Enabled debug optimizations")
    endif ()

    # Essential Emscripten link options for pybind11 modules
#         "SHELL:-s DISABLE_EXCEPTION_CATCHING=0"
    add_link_options(
        "SHELL:-s ALLOW_MEMORY_GROWTH=1"
        "SHELL:-s MAXIMUM_MEMORY=8GB"
        "SHELL:-s MODULARIZE=1"
        "SHELL:-s EXPORT_ES6=1"
    )
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # GCC-specific optimizations and warnings
    if (DEFINED ENV{CI})
        # Reduced warning set for CI to avoid hiding actual errors
        add_compile_options(
                -Wall -Wextra
                -Wformat=2 -Winit-self
                -Wmissing-declarations -Woverloaded-virtual
                -Wshadow -Wswitch-default -Wundef
                -Wno-unused-parameter  # Common in template-heavy code
        )
        message(STATUS "CI environment detected: Using reduced warning set for GCC")
    else ()
        # Full warning set for development
        add_compile_options(
                -Wall -Wextra -Wpedantic
                -Wcast-align -Wcast-qual -Wctor-dtor-privacy
                -Wdisabled-optimization -Wformat=2 -Winit-self
                -Wmissing-declarations -Wmissing-include-dirs
                -Wold-style-cast -Woverloaded-virtual -Wredundant-decls
                -Wshadow -Wsign-conversion -Wsign-promo
                -Wstrict-overflow=5 -Wswitch-default -Wundef
                -Wno-unused-parameter  # Common in template-heavy code
        )
    endif ()


    # Enable more aggressive optimization for Release builds
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(-flto)  # Link-time optimization
        set(CMAKE_AR gcc-ar)
        set(CMAKE_RANLIB gcc-ranlib)
    endif ()

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Clang-specific optimizations and warnings
    if (DEFINED ENV{CI})
        # Reduced warning set for CI to avoid hiding actual errors
        add_compile_options(
                -Wall -Wextra
                -Wformat=2 -Winit-self
                -Wmissing-declarations -Woverloaded-virtual
                -Wshadow -Wswitch-default -Wundef
                -Wno-unused-parameter
        )
        message(STATUS "CI environment detected: Using reduced warning set for Clang")
    else ()
        # Full warning set for development
        add_compile_options(
                -Wall -Wextra -Wpedantic
                -Wcast-align -Wcast-qual -Wctor-dtor-privacy
                -Wdisabled-optimization -Wformat=2 -Winit-self
                -Wmissing-declarations -Wold-style-cast -Woverloaded-virtual
                -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo
                -Wstrict-overflow=5 -Wswitch-default -Wundef
                -Wno-unused-parameter
        )
    endif ()


    # Enable LTO for Release builds
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(-flto)
    endif ()

elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # MSVC-specific optimizations
    if (DEFINED ENV{CI})
        # Reduced warning set for CI to avoid hiding actual errors
        add_compile_options(
                /W3                # Standard warning level instead of W4
                /permissive-       # Strict conformance
                /Zc:__cplusplus   # Correct __cplusplus macro
        )
        message(STATUS "CI environment detected: Using reduced warning set for MSVC")
    else ()
        # Full warning set for development
        add_compile_options(
                /W4                # High warning level
                /permissive-       # Strict conformance
                /Zc:__cplusplus   # Correct __cplusplus macro
        )
    endif ()


    # Suppress some noisy MSVC warnings
    add_compile_options(
            /wd4251  # 'identifier' : class 'type' needs to have dll-interface
            /wd4275  # non dll-interface class used as base for dll-interface class
    )

    # Enable whole program optimization for Release
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(/GL)  # Whole program optimization
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
        set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG")
    endif ()

endif ()

# Enable position-independent code for shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ==============================================================================
# Utility function to print compilation flags for debugging
# ==============================================================================

function(OUTPUT_FLAGS target_name)
    if (NOT TARGET ${target_name})
        message(WARNING "OUTPUT_FLAGS called with non-existent target: ${target_name}")
        return()
    endif ()

    get_target_property(COMPILE_OPTIONS ${target_name} COMPILE_OPTIONS)
    get_target_property(COMPILE_DEFS ${target_name} COMPILE_DEFINITIONS)
    get_target_property(INCLUDE_DIRS ${target_name} INCLUDE_DIRECTORIES)

    message(STATUS "=== ${target_name} Configuration ===")
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

    if (COMPILE_OPTIONS)
        message(STATUS "Compile options: ${COMPILE_OPTIONS}")
    endif ()

    if (COMPILE_DEFS)
        message(STATUS "Compile definitions: ${COMPILE_DEFS}")
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
    elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
        message(STATUS "Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        message(STATUS "RelWithDebInfo flags: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    endif ()

    message(STATUS "==================================")
endfunction()

message(STATUS "Compiler configuration completed for ${CMAKE_BUILD_TYPE} build")
