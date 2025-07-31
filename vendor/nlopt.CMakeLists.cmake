include(ExternalProject)
# include(FetchContent)

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif ()

# download, compile & install nlopt
# TODO: enable possibility of including externally installed nlopt library
if (NOT DEFINED nlopt_LOADED)
    if (NOT TARGET nloptLib)

        if (WIN32)
            set(LIBRARY_PREFIX "")
            set(LIBRARY_SUFFIX ".lib")
        else ()
            set(LIBRARY_PREFIX "lib")
            set(LIBRARY_SUFFIX ".a")
        endif ()

        set(nlopt_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/nlopt${vendor_suffix}")
        
        set(nlopt_BUILD_TYPE "Release")
        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(nlopt_BUILD_TYPE "Debug")
        endif ()

        ExternalProject_Add(
                nloptLib
                GIT_REPOSITORY https://github.com/stevengj/nlopt
                GIT_TAG 11cff2c773b4b98821915a72179f4667c307ce6d # 2.9.1
                PREFIX ${nlopt_PREFIX_PATH}
                INSTALL_DIR ${nlopt_PREFIX_PATH}/nloptLib-install
                CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${nlopt_PREFIX_PATH}/nloptLib-install -DINSTALL_LIBDIR=${nlopt_PREFIX_PATH}/nloptLib-install/lib -DCMAKE_INSTALL_LIBDIR=${nlopt_PREFIX_PATH}/nloptLib-install/lib -DNLOPT_GUILE=OFF -DNLOPT_OCTAVE=OFF -DNLOPT_MATLAB=OFF -DNLOPT_SWIG=OFF -DNLOPT_PYTHON=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=${nlopt_BUILD_TYPE}
                BUILD_COMMAND ${CMAKE_COMMAND} --build ${nlopt_PREFIX_PATH}/src/nloptLib-build --config ${nlopt_BUILD_TYPE}
                BUILD_BYPRODUCTS ${nlopt_PREFIX_PATH}/nloptLib-install/lib/${LIBRARY_PREFIX}nlopt${LIBRARY_SUFFIX}
        )
        # FetchContent_MakeAvailable(nloptLib)
        add_library(nlopt STATIC IMPORTED)
        add_dependencies(nlopt nloptLib)
        if (MSVC)
            set(nlopt_INCLUDE_DIRS "${nlopt_PREFIX_PATH}/nloptLib-install/include" "${nlopt_PREFIX_PATH}/src/nloptLib/msvc/include")
        else ()
            set(nlopt_INCLUDE_DIRS "${nlopt_PREFIX_PATH}/nloptLib-install/include")
        endif ()
        file(GLOB nlopt_LIBRARIES "${nlopt_PREFIX_PATH}/nloptLib-install/lib/${LIBRARY_PREFIX}nlopt*")
        if (NOT nlopt_LIBRARIES)
            # message("WARNING: nlopt_LIBRARIES empty")
            # TODO: this is somewhat unreliable
            set(nlopt_LIBRARIES "${nlopt_PREFIX_PATH}/nloptLib-install/lib/${LIBRARY_PREFIX}nlopt${LIBRARY_SUFFIX}")
            # file(GLOB_RECURSE nlopt_LIBRARIES "${nlopt_PREFIX_PATH}/*.a")
        endif ()
        message("Hoping nlopt_LIBRARIES will be compiled to: ${nlopt_LIBRARIES}")
        set_target_properties(nlopt PROPERTIES 
            IMPORTED_LOCATION ${nlopt_LIBRARIES}
#             INTERFACE_INCLUDE_DIRECTORIES "${nlopt_INCLUDE_DIRS}"
        )
        set(nlopt_LOADED ON)
    endif ()
endif ()

# find_package(nlopt REQUIRED)
# if(nlopt_FOUND)
#   include_directories(${nlopt_INCLUDE_DIRS})
#   target_link_libraries(pylimer_tools nlopt)
# 	message("Found nlopt for pylimer_tools_cpp")
# else()
# 	message(WARNING "DID NOT FIND nlopt")
# endif()
