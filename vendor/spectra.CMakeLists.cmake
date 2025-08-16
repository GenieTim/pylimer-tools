# load the Spectra library
include(FetchContent)
message(STATUS "Using FetchContent to load Spectra")

FetchContent_Declare(
        Spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG v1.0.1
)
FetchContent_MakeAvailable(Spectra)

SET(Spectra_INCLUDE_DIRS "${Spectra_SOURCE_DIR}/include")
