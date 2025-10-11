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

if (TARGET Eigen3::Eigen)
    add_dependencies(Spectra Eigen3::Eigen)
endif()

# Ensure the dependencies are set up correctly for export
if (TARGET eigen AND TARGET Spectra)
    # Make sure Spectra can find the eigen target for export
    set_target_properties(Spectra PROPERTIES
        INTERFACE_LINK_LIBRARIES "eigen"
    )
endif()
