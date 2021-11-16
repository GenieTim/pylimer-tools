# Install (incl. download) dependencies where applicable
include(ExternalProject)

# download, compile & install igraph
if (NOT DEFINED igraph_LOADED)
	if (NOT TARGET igraphLib)
	ExternalProject_Add(
			igraphLib
			GIT_REPOSITORY https://github.com/igraph/igraph.git
			GIT_TAG 0.9.5
			PREFIX ${CMAKE_CURRENT_LIST_DIR}/igraph
			INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install
			CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install
	)
	endif()
	# add_library(igraphLib)
	set(igraph_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/igraph")
	if (MSVC)
	set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include;${igraph_PREFIX_PATH}/src/igraphLib/msvc/include")
	else()
	set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include")
	endif()
	# set(igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/libigraph.a")
	file(GLOB igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/libigraph.*")
	set(igraph_LOADED ON)
endif()

# find_package(igraph REQUIRED)
# if(igraph_FOUND)
#   include_directories(${igraph_INCLUDE_DIRS})
#   target_link_libraries(pylimer_tools igraph)
# 	message("Found igraph for pylimer_tools_cpp")
# else()
# 	message(WARNING "DID NOT FIND igraph")
# endif()
