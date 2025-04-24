include(FetchContent)

if (NOT DEFINED cpp_optim_LOADED)
    if (NOT TARGET cpp_optim)

        if (WIN32)
            set(LIBRARY_PREFIX "")
            set(LIBRARY_SUFFIX ".lib")
        else ()
            set(LIBRARY_PREFIX "lib")
            set(LIBRARY_SUFFIX ".a")
        endif ()

        FetchContent_Declare(
                cpp_optim
                GIT_REPOSITORY https://github.com/PatWie/CppNumericalSolvers
                GIT_TAG 9e21e01736237a1e10e57bffca1b00fe6655d9f6 # origin/v2 #
        )

        FetchContent_MakeAvailable(cpp_optim)

        set(optim_INCLUDE_DIRS ${cpp_optim_SOURCE_DIR}/include)
        message("Hoping optim_INCLUDE_DIRS will be found in: ${optim_INCLUDE_DIRS}")

        set(cpp_optim_LOADED ON)
    endif ()
endif ()
