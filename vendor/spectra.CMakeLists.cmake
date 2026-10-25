# load the Spectra library
include(FetchContent)
message(STATUS "Using FetchContent to load Spectra")

if (NOT DEFINED FETCHCONTENT_SOURCE_DIR_SPECTRA AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/../tests/build/_deps/spectra-src")
    set(FETCHCONTENT_SOURCE_DIR_SPECTRA "${CMAKE_CURRENT_LIST_DIR}/../tests/build/_deps/spectra-src")
endif ()

FetchContent_Declare(
        Spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG v1.0.1
)
# FetchContent_MakeAvailable(Spectra)

# Get the Spectra source without building it yet so we can configure it
FetchContent_GetProperties(Spectra)
if(NOT spectra_POPULATED)
    # set policy CMP0169 to OLD in order to allow FetchContent_Populate without FetchContent_MakeAvailable
    cmake_policy(SET CMP0169 OLD)

    FetchContent_Populate(Spectra)
    
    # Temporarily override the install command to prevent Spectra from installing
    # This avoids the export conflict with FetchContent Eigen
    function(install)
        # Silently ignore all install commands from Spectra
    endfunction()
    
    # Silence Spectra tests/examples and exclude from all
    set(SPECTRA_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SPECTRA_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    
    # Add subdirectory excluding targets from 'all'
    add_subdirectory(${spectra_SOURCE_DIR} ${spectra_BINARY_DIR} EXCLUDE_FROM_ALL)
    
    # Restore the original install function
    unset(install)
endif()

SET(Spectra_INCLUDE_DIRS "${spectra_SOURCE_DIR}/include")
