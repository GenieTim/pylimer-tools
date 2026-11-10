# ==============================================================================
# OpenMP Configuration
# ==============================================================================
# Finds and configures OpenMP for multi-threaded parallel computation.
# On macOS with AppleClang, OpenMP is not bundled with Xcode Command Line Tools.
# Homebrew provides keg-only libomp (e.g. /opt/homebrew/opt/libomp on Apple Silicon
# or /usr/local/opt/libomp on Intel). Furthermore, system scientific libraries
# like openblas (linked via LAPACKE) are dynamically linked against Homebrew's
# libomp.dylib. If pylimer-tools links against a different libomp (such as one
# from a Homebrew LLVM installation picked up via PATH or LDFLAGS), dyld loads both
# runtime dylibs into memory, triggering fatal error:
#   "OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized."
#
# This module locates Homebrew's libomp, purges any stale LLVM libomp cache entries,
# and enforces linking against the matching Homebrew libomp runtime across all targets.
# ==============================================================================

if (NOT DEFINED openmp_LOADED)
    if (POLICY CMP0074)
        cmake_policy(SET CMP0074 NEW)  # find_package() uses PackageName_ROOT
    endif ()

    if (WIN32)
        message(STATUS "Windows detected: OpenMP disabled")
    elseif (CODE_COVERAGE)
        message(STATUS "Code coverage enabled: OpenMP disabled")
    else ()
        # On macOS with AppleClang or Clang, locate and enforce Homebrew's libomp
        if (APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
            set(_brew_libomp "")

            # 1. Check brew --prefix libomp
            find_program(BREW_EXECUTABLE brew)
            if (BREW_EXECUTABLE)
                execute_process(
                    COMMAND ${BREW_EXECUTABLE} --prefix libomp
                    OUTPUT_VARIABLE _brew_libomp_prefix
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                )
                if (_brew_libomp_prefix AND EXISTS "${_brew_libomp_prefix}/include/omp.h" AND EXISTS "${_brew_libomp_prefix}/lib/libomp.dylib")
                    set(_brew_libomp "${_brew_libomp_prefix}")
                endif ()
            endif ()

            # 2. Check standard Homebrew installation prefixes
            if (NOT _brew_libomp)
                foreach (_candidate IN ITEMS "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
                    if (EXISTS "${_candidate}/include/omp.h" AND EXISTS "${_candidate}/lib/libomp.dylib")
                        set(_brew_libomp "${_candidate}")
                        break ()
                    endif ()
                endforeach ()
            endif ()

            # 3. Check OpenMP_ROOT or ENV{OpenMP_ROOT} if explicitly set and not pointing to LLVM
            if (NOT _brew_libomp)
                foreach (_hint IN ITEMS "${OpenMP_ROOT}" "$ENV{OpenMP_ROOT}")
                    if (_hint AND EXISTS "${_hint}/include/omp.h" AND EXISTS "${_hint}/lib/libomp.dylib")
                        if (NOT _hint MATCHES "llvm")
                            set(_brew_libomp "${_hint}")
                            break ()
                        endif ()
                    endif ()
                endforeach ()
            endif ()

            if (_brew_libomp)
                # Purge stale cache entries if pointing to LLVM to prevent duplicate runtime errors
                if (DEFINED OpenMP_libomp_LIBRARY AND OpenMP_libomp_LIBRARY MATCHES "llvm")
                    message(STATUS "Overriding stale LLVM libomp (${OpenMP_libomp_LIBRARY}) with Homebrew libomp (${_brew_libomp})")
                    unset(OpenMP_libomp_LIBRARY CACHE)
                    unset(OpenMP_CXX_FLAGS CACHE)
                    unset(OpenMP_CXX_LIB_NAMES CACHE)
                    unset(OpenMP_C_FLAGS CACHE)
                    unset(OpenMP_C_LIB_NAMES CACHE)
                    unset(FIND_PACKAGE_MESSAGE_DETAILS_OpenMP CACHE)
                    unset(FIND_PACKAGE_MESSAGE_DETAILS_OpenMP_CXX CACHE)
                    unset(FIND_PACKAGE_MESSAGE_DETAILS_OpenMP_C CACHE)
                endif ()

                # Enforce Homebrew libomp in CMake cache
                set(OpenMP_ROOT "${_brew_libomp}" CACHE PATH "OpenMP root directory" FORCE)
                set(OpenMP_CXX_INCLUDE_DIR "${_brew_libomp}/include" CACHE PATH "OpenMP CXX include directory" FORCE)
                set(OpenMP_C_INCLUDE_DIR "${_brew_libomp}/include" CACHE PATH "OpenMP C include directory" FORCE)
                set(OpenMP_libomp_LIBRARY "${_brew_libomp}/lib/libomp.dylib" CACHE FILEPATH "OpenMP libomp library" FORCE)
                list(PREPEND CMAKE_PREFIX_PATH "${_brew_libomp}")
                message(STATUS "Using macOS Homebrew libomp at: ${_brew_libomp}")
            endif ()
        endif ()

        find_package(OpenMP)
        if (OpenMP_CXX_FOUND)
            # On macOS, ensure imported targets explicitly link against the intended libomp dylib
            # so that linker search flags (-L/opt/homebrew/opt/llvm/lib) cannot redirect the link.
            if (APPLE AND _brew_libomp)
                if (TARGET OpenMP::OpenMP_CXX)
                    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY
                        INTERFACE_LINK_LIBRARIES "${_brew_libomp}/lib/libomp.dylib"
                    )
                    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY
                        INTERFACE_INCLUDE_DIRECTORIES "${_brew_libomp}/include"
                    )
                endif ()
                if (TARGET OpenMP::OpenMP_C)
                    set_property(TARGET OpenMP::OpenMP_C PROPERTY
                        INTERFACE_LINK_LIBRARIES "${_brew_libomp}/lib/libomp.dylib"
                    )
                    set_property(TARGET OpenMP::OpenMP_C PROPERTY
                        INTERFACE_INCLUDE_DIRECTORIES "${_brew_libomp}/include"
                    )
                endif ()
            endif ()
            message(STATUS "OpenMP found and enabled (version ${OpenMP_CXX_VERSION})")
        else ()
            message(STATUS "OpenMP not found")
        endif ()
    endif ()

    set(openmp_LOADED ON)
endif ()
