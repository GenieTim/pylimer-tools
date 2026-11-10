# ==============================================================================
# OpenMP Configuration
# ==============================================================================
# Finds and configures OpenMP for multi-threaded parallel computation.
# On macOS with AppleClang, OpenMP is not bundled with Xcode Command Line Tools.
# Homebrew provides keg-only libomp (e.g. /opt/homebrew/opt/libomp on Apple Silicon
# or /usr/local/opt/libomp on Intel). This module locates Homebrew's libomp and
# configures OpenMP_ROOT, OpenMP_CXX_INCLUDE_DIR, and OpenMP_libomp_LIBRARY so
# CMake finds the matching runtime and prevents conflicts with other libomp builds.
# ==============================================================================

if (NOT DEFINED openmp_LOADED)
    if (WIN32)
        message(STATUS "Windows detected: OpenMP disabled")
    elseif (CODE_COVERAGE)
        message(STATUS "Code coverage enabled: OpenMP disabled")
    else ()
        # On macOS with AppleClang or Clang, point to Homebrew's libomp
        if (APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
            set(_libomp_hints "")
            if (DEFINED OpenMP_ROOT AND EXISTS "${OpenMP_ROOT}")
                list(APPEND _libomp_hints "${OpenMP_ROOT}")
            endif ()
            if (DEFINED ENV{OpenMP_ROOT} AND EXISTS "$ENV{OpenMP_ROOT}")
                list(APPEND _libomp_hints "$ENV{OpenMP_ROOT}")
            endif ()

            # Query brew --prefix libomp if brew is available
            find_program(BREW_EXECUTABLE brew)
            if (BREW_EXECUTABLE)
                execute_process(
                    COMMAND ${BREW_EXECUTABLE} --prefix libomp
                    OUTPUT_VARIABLE _brew_libomp_prefix
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                )
                if (_brew_libomp_prefix AND EXISTS "${_brew_libomp_prefix}")
                    list(APPEND _libomp_hints "${_brew_libomp_prefix}")
                endif ()
            endif ()

            # Standard Homebrew fallback paths
            list(APPEND _libomp_hints "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")

            foreach (_hint IN LISTS _libomp_hints)
                if (EXISTS "${_hint}/include/omp.h" AND EXISTS "${_hint}/lib/libomp.dylib")
                    if (NOT DEFINED OpenMP_ROOT)
                        set(OpenMP_ROOT "${_hint}" CACHE PATH "OpenMP root directory" FORCE)
                    endif ()
                    if (NOT DEFINED OpenMP_CXX_INCLUDE_DIR)
                        set(OpenMP_CXX_INCLUDE_DIR "${_hint}/include" CACHE PATH "OpenMP CXX include directory" FORCE)
                    endif ()
                    if (NOT DEFINED OpenMP_C_INCLUDE_DIR)
                        set(OpenMP_C_INCLUDE_DIR "${_hint}/include" CACHE PATH "OpenMP C include directory" FORCE)
                    endif ()
                    if (NOT DEFINED OpenMP_libomp_LIBRARY)
                        set(OpenMP_libomp_LIBRARY "${_hint}/lib/libomp.dylib" CACHE FILEPATH "OpenMP libomp library" FORCE)
                    endif ()
                    message(STATUS "Found macOS Homebrew libomp at: ${_hint}")
                    break ()
                endif ()
            endforeach ()
        endif ()

        find_package(OpenMP)
        if (OpenMP_CXX_FOUND)
            message(STATUS "OpenMP found and enabled (version ${OpenMP_CXX_VERSION})")
        else ()
            message(STATUS "OpenMP not found")
        endif ()
    endif ()

    set(openmp_LOADED ON)
endif ()
