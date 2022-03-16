# Install (incl. download) dependencies where applicable
include(ExternalProject)
# include(FetchContent)

# download, compile & install igraph
# TODO: enable possibility of including externally installed igraph library
if (NOT DEFINED igraph_LOADED)
	if (NOT TARGET igraphLib)

		if(WIN32)
			set(LIBRARY_PREFIX "")
			set(LIBRARY_SUFFIX ".lib")
		else()
			set(LIBRARY_PREFIX "lib")
			set(LIBRARY_SUFFIX ".a")
		endif()

		ExternalProject_Add(
				igraphLib
				GIT_REPOSITORY https://github.com/igraph/igraph.git
				GIT_TAG c20705b81cfa0a84cb110931a3958cd71110b216 # 0.9.7
				PREFIX ${CMAKE_CURRENT_LIST_DIR}/igraph
				INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install
				CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install -DCMAKE_INSTALL_LIBDIR=${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install/lib -DCMAKE_C_FLAGS="-fPIC"
				BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_LIST_DIR}/igraph/src/igraphLib-build --config Release
				BUILD_BYPRODUCTS ${CMAKE_CURRENT_LIST_DIR}/igraph/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}
		)
		# FetchContent_MakeAvailable(igraphLib)
		add_library(igraph SHARED IMPORTED)
		add_dependencies(igraph igraphLib)
		set(igraph_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/igraph")
		if (MSVC)
			set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include;${igraph_PREFIX_PATH}/src/igraphLib/msvc/include")
		else()
			set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include")
		endif()
		file(GLOB igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph.*")
		if (NOT igraph_LIBRARIES)
			# message("WARNING: igraph_LIBRARIES empty")
			set(igraph_LIBRARIES "${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}")
			# file(GLOB_RECURSE igraph_LIBRARIES "${igraph_PREFIX_PATH}/*.a")
		endif()
		message("Hoping igraph_LIBRARIES will be compiled to: ${igraph_LIBRARIES}")
		set_target_properties(igraph PROPERTIES IMPORTED_LOCATION ${igraph_LIBRARIES})
		set(igraph_LOADED ON)
	endif()
endif()

# find_package(igraph REQUIRED)
# if(igraph_FOUND)
#   include_directories(${igraph_INCLUDE_DIRS})
#   target_link_libraries(pylimer_tools igraph)
# 	message("Found igraph for pylimer_tools_cpp")
# else()
# 	message(WARNING "DID NOT FIND igraph")
# endif()
