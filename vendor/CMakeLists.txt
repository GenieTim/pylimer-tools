# Install (incl. download) dependencies where applicable

if (NOT DEFINED vendor_suffix)
    set(vendor_suffix "")
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

include(${CMAKE_CURRENT_LIST_DIR}/myCompilerDefaultFlags.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/igraph.CMakeLists.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/cppNumericalSolvers.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/nlopt.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/eigen.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cereal.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/spectra.CMakeLists.cmake)

if (TARGET Eigen3)
    add_dependencies(Spectra Eigen3)
endif()

add_dependencies(Spectra Eigen3::Eigen)
