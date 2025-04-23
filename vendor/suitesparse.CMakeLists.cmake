include(FetchContent)

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif ()

if (NOT DEFINED suitesparse_LOADED)

    set(suitesparse_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/suitesparse${vendor_suffix}")

    include(ExternalProject)
    ExternalProject_Add(suitesparseLib
            GIT_REPOSITORY git@github.com:DrTimothyAldenDavis/SuiteSparse.git
            GIT_SUBMODULES_RECURSE ON
            GIT_TAG stable
            GIT_SHALLOW ON
            CMAKE_ARGS -DJUST_INSTALL_suitesparse=ON -DSKIP_PORTABILITY_TEST=ON -DBUILD_TESTS=OFF
            PREFIX ${suitesparse_PREFIX_PATH}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${suitesparse_PREFIX_PATH}/suitesparseLib-install -DINSTALL_LIBDIR=${CMAKE_CURRENT_LIST_DIR}/suitesparse/suitesparseLib-install/lib -DCMAKE_INSTALL_LIBDIR=${CMAKE_CURRENT_LIST_DIR}/suitesparse/suitesparseLib-install/lib -Dsuitesparse_GUILE=OFF -Dsuitesparse_OCTAVE=OFF -Dsuitesparse_MATLAB=OFF -Dsuitesparse_SWIG=OFF -Dsuitesparse_PYTHON=OFF -DBUILD_SHARED_LIBS=OFF
            INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/suitesparse${vendor_suffix}/suitesparseLib-install
    )

    # Create an INTERFACE library instead of STATIC IMPORTED
    add_library(suitesparse INTERFACE)
    add_dependencies(suitesparse suitesparseLib)

    # find the compiled library files
    file(GLOB suitesparse_LIBRARIES
            "${suitesparse_PREFIX_PATH}/suitesparseLib-install/lib/*.a"
            "${suitesparse_PREFIX_PATH}/suitesparseLib-install/lib/*.lib"
            "${suitesparse_PREFIX_PATH}/suitesparseLib-install/lib/*.so"
            "${suitesparse_PREFIX_PATH}/suitesparseLib-install/lib/*.dylib")

    set(suitesparse_INCLUDE_DIRS "${suitesparse_PREFIX_PATH}/suitesparseLib-install/include/suitesparse")
    # Set the include directories and link libraries as interface properties
    set_target_properties(suitesparse PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${suitesparse_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${suitesparse_LIBRARIES}")

    set(suitesparse_LOADED ON)
endif ()