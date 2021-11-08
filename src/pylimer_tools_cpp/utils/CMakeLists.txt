add_library(pylimer_utils INTERFACE)

set(CMAKE_CXX_STANDARD 20)

# add libraries
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
  target_link_libraries(pylimer_utils ${igraph_LIBRARIES})
	message("Found igraph")
else()
	message(WARNING "DID NOT FIND igraph")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pylimer_utils ${Boost_LIBRARIES})
endif()

target_sources(pylimer_utils
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/DataFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/DumpFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/GraphUtils.h
		${CMAKE_CURRENT_LIST_DIR}/VectorUtils.h
)

set_target_properties(pylimer_utils PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(pylimer_utils
	INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}
)
