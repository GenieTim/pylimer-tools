# Install (incl. download) dependencies where applicable

if (NOT DEFINED vendor_suffix)
	set(vendor_suffix "")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/igraph.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cppNumericalSolvers.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/nlopt.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/eigen.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cereal.CMakeLists.cmake)
