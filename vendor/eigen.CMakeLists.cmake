# load the Eigen library
if (NOT DEFINED eigen_LOADED)
    find_package(Eigen3 3.4 NO_MODULE) # 3.4
    if (TARGET Eigen3::Eigen)  #(${Eigen3_FOUND}) # AND (${Eigen3_VERSION} VERSION_GREATER_EQUAL 3.4)
        message(STATUS "Found Eigen3 Version: ${Eigen3_VERSION} Path: ${Eigen3_DIR}")
    else ()
        message(STATUS "Using FetchContent to load Eigen3")
        include(FetchContent)
        FetchContent_Declare(
                Eigen
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen
                GIT_TAG 3.4.0
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
        )
        set(EIGEN_BUILD_DOC OFF)
        set(EIGEN_BUILD_PKGCONFIG OFF)
        set(EIGEN_BUILD_TESTING OFF)
        set(EIGEN_MPL2_ONLY ON)

        FetchContent_MakeAvailable(Eigen)

        # Create alias target if it doesn't exist
        if(NOT TARGET Eigen3::Eigen)
            add_library(Eigen3::Eigen ALIAS eigen)
        endif()

        # Synchronize Eigen3_* variables from Eigen_* variables if needed
        get_cmake_property(_allVariables VARIABLES)
        foreach(_var ${_allVariables})
            if(_var MATCHES "^Eigen_(.+)$")
                set(_suffix ${CMAKE_MATCH_1})
                set(_eigen3_var "Eigen3_${_suffix}")

                # Check if Eigen3_* variable needs to be set
                if(NOT DEFINED ${_eigen3_var} OR
                   "${${_eigen3_var}}" STREQUAL "" OR
                   "${${_eigen3_var}}" STREQUAL "OFF" OR
                   "${${_eigen3_var}}" MATCHES ".*-NOTFOUND$")
                    set(${_eigen3_var} ${${_var}})
                    message(STATUS "Set ${_eigen3_var} to ${${_var}}")
                endif()
            endif()
        endforeach()

        if (NOT EXISTS ${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/eigen3-config.cmake AND
                NOT EXISTS ${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/Eigen3Config.cmake)
            file(WRITE ${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/eigen3-config.cmake
                    [=[
# Redirect to the FetchContent Eigen target
if(NOT TARGET Eigen3::Eigen AND TARGET eigen)
    add_library(Eigen3::Eigen ALIAS eigen)
endif()
set(Eigen3_FOUND TRUE)
set(EIGEN3_FOUND TRUE)
            ]=])
        endif ()
    endif ()

    set(eigen_LOADED ON)
    message(STATUS "Eigen include directories: ${eigen_INCLUDE_DIRS}, libraries ${eigen_LIBRARIES}")
endif ()

# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKEXT.cmake)
find_package(LAPACKE)
# include(FindLAPACKLibs)

# if (LAPACKLIBS_FOUND)
if (LAPACKE_FOUND AND LAPACKE_INCLUDE_DIRS)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_LIST_DIR}/lapacke-extra)
    add_definitions(-DEIGEN_USE_LAPACKE)
    add_definitions(-DHAVE_LAPACK_CONFIG_H)
    add_definitions(-DLAPACK_COMPLEX_CPP)
    # add_definitions(-DLAPACK_ILP64) # not supported by Eigen, unfortunately
    add_definitions(-DLAPACK_DISABLE_NAN_CHECK)
    message(STATUS "Found lapacke, include dirs ${LAPACKE_INCLUDE_DIRS}, libraries ${LAPACKE_LIBRARIES}")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(BLA_VENDOR Intel10_64lp)
    add_definitions(-DEIGEN_USE_MKL_ALL)
endif ()
