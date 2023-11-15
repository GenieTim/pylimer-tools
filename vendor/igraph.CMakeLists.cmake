include(ExternalProject)
# include(FetchContent)

if (NOT DEFINED vendor_suffix)
	set(vendor_suffix "")
endif()

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

		set(igraph_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/igraph${vendor_suffix}")

		ExternalProject_Add(
				igraphLib
				GIT_REPOSITORY https://github.com/igraph/igraph.git
				GIT_TAG 6559f7e92f64d6f71a61063132dace9ce72cf680 # 0.10.7
				PREFIX ${igraph_PREFIX_PATH}
				INSTALL_DIR ${igraph_PREFIX_PATH}/igraphLib-install
				CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${igraph_PREFIX_PATH}/igraphLib-install -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_LIBDIR=${igraph_PREFIX_PATH}/igraphLib-install/lib -DIGRAPH_GRAPHML_SUPPORT=ON -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true
				BUILD_COMMAND ${CMAKE_COMMAND} --build ${igraph_PREFIX_PATH}/src/igraphLib-build --config Release
				BUILD_BYPRODUCTS ${igraph_PREFIX_PATH}/igraphLib-install/lib/${LIBRARY_PREFIX}igraph${LIBRARY_SUFFIX}
				UPDATE_DISCONNECTED ON
		)
		# FetchContent_MakeAvailable(igraphLib)
		add_library(igraph SHARED IMPORTED)
		add_dependencies(igraph igraphLib)
		if (MSVC)
			set(igraph_INCLUDE_DIRS "${igraph_PREFIX_PATH}/igraphLib-install/include" "${igraph_PREFIX_PATH}/src/igraphLib/msvc/include")
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
