include(FetchContent)

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif()

if (NOT DEFINED least_squares_cpp_LOADED)
FetchContent_Declare(least_squares_cpp
        GIT_REPOSITORY git@github.com:Rookfighter/least-squares-cpp.git
        GIT_TAG <the git commit>
        SOURCE_SUBDIR "MADE-UP-DIRECTORY"
)

FetchContent_MakeAvailable(the-project)

# I recommend using SYSTEM when dealing with 3rd party code.
# Avoids the hassle of warnings from a library you don't own.
target_include_directories(foo SYSTEM PRIVATE
        "${FETCHCONTENT_BASE_DIR}/the-project-src/include"
)
