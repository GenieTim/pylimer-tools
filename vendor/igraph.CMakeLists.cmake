include(ExternalProject)
# include(FetchContent)

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif ()

# download, compile & install igraph
if (NOT DEFINED igraph_LOADED)
    find_package(igraph "0.10.15...<1.1.0")

    # if (${igraph_FOUND} AND TARGET igraph::igraph)
    # 	message("Found igraph library")
    # 	set(igraph_LIBRARIES igraph::igraph)
    # 	get_target_property(l_igraph_INCLUDE_DIRS igraph::igraph INTERFACE_INCLUDE_DIRECTORIES)
    # 	set(igraph_INCLUDE_DIRS "${l_igraph_INCLUDE_DIRS}")
    # else()
    if (NOT TARGET igraphLib)

        if (WIN32)
            set(LIBRARY_PREFIX "")
            set(LIBRARY_SUFFIX ".lib")
        else ()
            set(LIBRARY_PREFIX "lib")
            set(LIBRARY_SUFFIX ".a")
        endif ()

        set(igraph_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/igraph${vendor_suffix}")
        set(igraph_EXTRA_CMAKE_ARGS "")
        if (DEFINED ENV{BISON_EXECUTABLE})
            set(igraph_EXTRA_CMAKE_ARGS "-DBISON_EXECUTABLE=$ENV{BISON_EXECUTABLE}")
        endif ()

        # For WebAssembly / Pyodide builds we cannot link against system GMP.
        # igraph optionally uses GMP for big integer support; disable it to
        # avoid pulling in a non-wasm libgmp.a from the host toolchain.
        if (BUILDING_WITH_PYODIDE)
            # Aggressively disable GMP: igraph changed option names across versions; provide all.
            list(APPEND igraph_EXTRA_CMAKE_ARGS
                -DIGRAPH_USE_INTERNAL_BLAS=ON
                -DIGRAPH_USE_INTERNAL_LAPACK=ON
                -DIGRAPH_USE_INTERNAL_ARPACK=ON
                -DIGRAPH_USE_INTERNAL_GLPK=ON
                -DIGRAPH_USE_INTERNAL_GMP=ON
                -DIGRAPH_USE_INTERNAL_PLFIT=ON
                -DHAVE_GMP=0
                #
                -DIGRAPH_GRAPHML_SUPPORT=OFF
                -DIGRAPH_WARNINGS_AS_ERRORS=OFF
                -DIGRAPH_ENABLE_LTO=AUTO
            )
        endif ()
        
        # Add Windows-specific configuration
        if (WIN32)
            list(APPEND igraph_EXTRA_CMAKE_ARGS 
#                 "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
#                 "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                "-DCMAKE_POSITION_INDEPENDENT_CODE=OFF"
            )
        else ()
            list(APPEND igraph_EXTRA_CMAKE_ARGS "-DCMAKE_POSITION_INDEPENDENT_CODE=ON")
        endif ()

        set(igraph_BUILD_TYPE "Release")
        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(igraph_BUILD_TYPE "Debug")
        elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
            set(igraph_BUILD_TYPE "RelWithDebInfo")
        endif ()

        # Set expected build products for different platforms and build types
        if (WIN32)
            set(EXPECTED_IGRAPH_LIB "${igraph_PREFIX_PATH}/igraphLib-install/lib/igraph.lib")
        else ()
            set(EXPECTED_IGRAPH_LIB "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}")
        endif ()

        ExternalProject_Add(
                igraphLib
                GIT_REPOSITORY https://github.com/GenieTim/igraph.git # https://github.com/igraph/igraph.git
                GIT_TAG ca4def86c89f1a3b9a15003b14abde60bc864931 # 0.10.15
                PREFIX ${igraph_PREFIX_PATH}
                PATCH_COMMAND git apply --check ${CMAKE_CURRENT_LIST_DIR}/patches/igraph.patch && git apply ${CMAKE_CURRENT_LIST_DIR}/patches/igraph.patch || true
                INSTALL_DIR ${igraph_PREFIX_PATH}/igraphLib-install
                CMAKE_ARGS ${igraph_EXTRA_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${igraph_PREFIX_PATH}/igraphLib-install -DCMAKE_BUILD_TYPE=${igraph_BUILD_TYPE} -DCMAKE_INSTALL_LIBDIR=${igraph_PREFIX_PATH}/igraphLib-install/lib -DIGRAPH_GRAPHML_SUPPORT=OFF
                BUILD_COMMAND ${CMAKE_COMMAND} --build ${igraph_PREFIX_PATH}/src/igraphLib-build --config ${igraph_BUILD_TYPE}
                BUILD_BYPRODUCTS ${EXPECTED_IGRAPH_LIB}
        )
        # FetchContent_MakeAvailable(igraphLib)
        
        add_library(igraph::igraph STATIC IMPORTED)
        add_dependencies(igraph::igraph igraphLib)

        if (MSVC)
            set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include" "${igraph_PREFIX_PATH}/src/igraphLib/msvc/include")
        else ()
            set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include")
        endif ()
        
        # Find the actual library files with more robust patterns
        if (WIN32)
            # On Windows, try multiple possible names and locations
            file(GLOB igraph_LIBRARIES 
                "${igraph_PREFIX_PATH}/igraphLib-install/lib/igraph.lib"
                "${igraph_PREFIX_PATH}/igraphLib-install/lib/libigraph.lib"
                "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}"
                "${igraph_PREFIX_PATH}/igraphLib-install/lib/Release/igraph.lib"
                "${igraph_PREFIX_PATH}/igraphLib-install/lib/Debug/igraph.lib"
            )
        else ()
            file(GLOB igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph.*")
        endif ()
        
        if (NOT igraph_LIBRARIES)
            # message("WARNING: igraph_LIBRARIES empty")
            set(igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}")
            # file(GLOB_RECURSE igraph_LIBRARIES "${igraph_PREFIX_PATH}/*.a")
        endif ()
        message("Hoping igraph_LIBRARIES will be compiled to: ${igraph_LIBRARIES}")
        
        # Use the first found library if multiple exist
        list(GET igraph_LIBRARIES 0 igraph_LIBRARY_MAIN)
        if (EXISTS "${igraph_LIBRARY_MAIN}")
            set(igraph_LIBRARIES "${igraph_LIBRARY_MAIN}")
            message(STATUS "Using igraph library: ${igraph_LIBRARIES}")
        else ()
            message(WARNING "igraph library not found at expected location: ${igraph_LIBRARY_MAIN}")
        endif ()
        
        set_target_properties(igraph::igraph PROPERTIES 
            IMPORTED_LOCATION "${igraph_LIBRARIES}"
#             INTERFACE_INCLUDE_DIRECTORIES "${igraph_INCLUDE_DIRS}"
        )
        
        set(igraph_LOADED ON)
    endif ()
    # endif()
endif ()

# find_package(igraph REQUIRED)
# if(igraph_FOUND)
#   include_directories(${igraph_INCLUDE_DIRS})
#   target_link_libraries(pylimer_tools igraph)
# 	message("Found igraph for pylimer_tools_cpp")
# else()
# 	message(WARNING "DID NOT FIND igraph")
# endif()
