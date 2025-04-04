# Install (incl. download) dependencies where applicable

if (NOT DEFINED vendor_suffix)
	set(vendor_suffix "")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

include(${CMAKE_CURRENT_LIST_DIR}/myCompilerDefaultFlags.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/CPM.cmake)

# CPMAddPackage(
# 	NAME igraph
# 	# GITHUB_REPOSITORY GenieTim/igraph
# 	# VERSION 0.10.15
# 	# GIT_TAG 6f8c9394cbe7f4986861befd607850b91492ace1
# 	URL https://github.com/GenieTim/igraph/archive/6f8c9394cbe7f4986861befd607850b91492ace1.tar.gz
# )

CPMAddPackage(
	NAME igraph
	# GITHUB_REPOSITORY igraph/igraph
	# VERSION 0.10.15
	# GIT_TAG 635b432eff0a89580ac9bb98068d2fbc8ef374f2
	# GIT_SHALLOW OFF
	URL https://github.com/igraph/igraph/releases/download/0.10.15/igraph-0.10.15.tar.gz
)

if (NOT DEFINED igraph_LIBRARIES)
	set(igraph_LIBRARIES igraph)
	set(igraph_INCLUDE_DIRS ${igraph_SOURCE_DIR}/include)
endif()

CPMAddPackage("gh:GenieTim/cereal#8a8ef6bd5e5416de839f1fb164b734292d2641e4")

if (NOT DEFINED cereal_LIBRARIES)
	set(cereal_LIBRARIES cereal)
	set(cereal_INCLUDE_DIRS ${cereal_SOURCE_DIR}/include)
endif()

CPMAddPackage("gh:stevengj/nlopt#11cff2c773b4b98821915a72179f4667c307ce6d")

if (NOT DEFINED nlopt_LIBRARIES)
	set(nlopt_LIBRARIES nlopt)
	set(nlopt_INCLUDE_DIRS ${nlopt_BINARY_DIR})
endif()

CPMAddPackage("gh:yixuan/spectra#v1.0.1")

if (NOT DEFINED Spectra_LIBRARIES)
	set(Spectra_LIBRARIES Spectra)
	set(Spectra_INCLUDE_DIRS ${spectra_SOURCE_DIR}/include)
endif()

set(USE_CCACHE ON)

CPMAddPackage(
  NAME Ccache.cmake
  GITHUB_REPOSITORY TheLartians/Ccache.cmake
  VERSION 1.2.5
)

# include(${CMAKE_CURRENT_LIST_DIR}/igraph.CMakeLists.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/cppNumericalSolvers.CMakeLists.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/nlopt.CMakeLists.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/eigen.CMakeLists.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/cereal.CMakeLists.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/spectra.CMakeLists.cmake)
