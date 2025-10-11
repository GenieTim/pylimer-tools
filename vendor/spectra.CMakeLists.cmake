# load the Spectra library
include(FetchContent)
message(STATUS "Using FetchContent to load Spectra")

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
    
    # Now add the subdirectory with our overridden install function
    add_subdirectory(${spectra_SOURCE_DIR} ${spectra_BINARY_DIR})
    
    # Restore the original install function
    unset(install)
endif()

SET(Spectra_INCLUDE_DIRS "${spectra_SOURCE_DIR}/include")
