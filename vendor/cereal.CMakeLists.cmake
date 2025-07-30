include(FetchContent)

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif ()

if (NOT DEFINED cereal_LOADED)
    # include(FetchContent)
    # FetchContent_Declare(
    # 	cereal
    # 	GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    # 	GIT_TAG        v1.3.2
    # )
    # FetchContent_MakeAvailable(cereal)
    set(cereal_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/cereal${vendor_suffix}")

    include(ExternalProject)
    ExternalProject_Add(cerealLib
            GIT_REPOSITORY https://github.com/GenieTim/cereal.git # https://github.com/USCiLab/cereal.git
            GIT_SUBMODULES_RECURSE ON
            GIT_SHALLOW ON
            CMAKE_ARGS -DJUST_INSTALL_CEREAL=ON -DSKIP_PORTABILITY_TEST=ON -DBUILD_TESTS=OFF
            PREFIX ${cereal_PREFIX_PATH}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${cereal_PREFIX_PATH}/cerealLib-install -DINSTALL_LIBDIR=${CMAKE_CURRENT_LIST_DIR}/cereal/cerealLib-install/lib -DCMAKE_INSTALL_LIBDIR=${CMAKE_CURRENT_LIST_DIR}/cereal/cerealLib-install/lib -Dcereal_GUILE=OFF -Dcereal_OCTAVE=OFF -Dcereal_MATLAB=OFF -Dcereal_SWIG=OFF -Dcereal_PYTHON=OFF -DBUILD_SHARED_LIBS=OFF
            INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/cereal${vendor_suffix}/cerealLib-install
    )

    # add_library(cereal STATIC IMPORTED)
    # add_dependencies(cereal cerealLib)
    # set(cereal_INCLUDE_DIRS "${cereal_PREFIX_PATH}/cerealLib-install/include")
    # file(GLOB cereal_LIBRARIES "${cereal_PREFIX_PATH}/cerealLib-install/lib/${LIBRARY_PREFIX}cereal*")
    # set_target_properties(cereal PROPERTIES IMPORTED_LOCATION ${cereal_LIBRARIES})
    set(cereal_INCLUDE_DIRS "${cereal_PREFIX_PATH}/src/cerealLib/include")

    set(cereal_LOADED ON)
endif ()


if (NOT DEFINED CEREALIZABLE)
    option(CEREALIZABLE "Enable serialisation of various classes" ON)
endif ()
if (CEREALIZABLE)
    add_definitions(-DCEREALIZABLE)
    add_compile_definitions(CEREALIZABLE)
    message(STATUS "Enabling serialization for classes in the project. Use the CEREALIZABLE CMake option to enable/disable this feature.")
else ()
    message(STATUS "Disabling serialization for classes in the project. Use the CEREALIZABLE CMake option to enable/disable this feature.")
endif ()
