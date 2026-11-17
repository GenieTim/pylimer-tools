# ==============================================================================
# Common CMake arguments forwarded to ExternalProject vendor dependencies
# ==============================================================================

if (NOT DEFINED VENDOR_FORWARD_CMAKE_ARGS)
    set(VENDOR_FORWARD_CMAKE_ARGS
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
        "-DCMAKE_C_VISIBILITY_PRESET=hidden"
        "-DCMAKE_CXX_VISIBILITY_PRESET=hidden"
        "-DCMAKE_VISIBILITY_INLINES_HIDDEN=ON"
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

    # Forward the WebAssembly flags defined in myCompilerDefaultFlags.cmake so
    # that vendor libraries share the ABI (wasm exceptions), visibility and
    # optimization settings of the main module.
    if (DEFINED PYLIMER_WASM_COMPILE_FLAGS)
        list(JOIN PYLIMER_WASM_COMPILE_FLAGS " " _vendor_wasm_common_flags)
        list(JOIN PYLIMER_WASM_LINK_FLAGS " " _vendor_wasm_link_flags)

        set(_vendor_c_base "${CMAKE_C_FLAGS}")
        if (NOT _vendor_c_base AND DEFINED ENV{CFLAGS})
            set(_vendor_c_base "$ENV{CFLAGS}")
        endif ()
        set(_vendor_cxx_base "${CMAKE_CXX_FLAGS}")
        if (NOT _vendor_cxx_base AND DEFINED ENV{CXXFLAGS})
            set(_vendor_cxx_base "$ENV{CXXFLAGS}")
        endif ()
        set(_vendor_exe_ld_base "${CMAKE_EXE_LINKER_FLAGS}")
        if (NOT _vendor_exe_ld_base AND DEFINED ENV{LDFLAGS})
            set(_vendor_exe_ld_base "$ENV{LDFLAGS}")
        endif ()

        list(APPEND VENDOR_FORWARD_CMAKE_ARGS
            "-DCMAKE_C_FLAGS=${_vendor_c_base} ${_vendor_wasm_common_flags}"
            "-DCMAKE_CXX_FLAGS=${_vendor_cxx_base} ${_vendor_wasm_common_flags} -fvisibility-inlines-hidden"
            "-DCMAKE_EXE_LINKER_FLAGS=${_vendor_exe_ld_base} ${_vendor_wasm_link_flags}"
        )
    endif ()
endif ()
