# ==============================================================================
# Common CMake arguments forwarded to ExternalProject vendor dependencies
# ==============================================================================

if (NOT DEFINED VENDOR_FORWARD_CMAKE_ARGS)
    set(VENDOR_FORWARD_CMAKE_ARGS
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
    )
    if (CMAKE_TOOLCHAIN_FILE)
        list(APPEND VENDOR_FORWARD_CMAKE_ARGS "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    endif ()
    if (CMAKE_OSX_ARCHITECTURES)
        list(APPEND VENDOR_FORWARD_CMAKE_ARGS "-DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}")
    endif ()
    if (CMAKE_OSX_SYSROOT)
        list(APPEND VENDOR_FORWARD_CMAKE_ARGS "-DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}")
    endif ()
    if (CMAKE_C_COMPILER)
        list(APPEND VENDOR_FORWARD_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
    endif ()
    if (CMAKE_CXX_COMPILER)
        list(APPEND VENDOR_FORWARD_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
    endif ()
endif ()
