add_library(pylimer_tools_bindings)

set_target_properties(pylimer_tools_bindings PROPERTIES LINKER_LANGUAGE CXX)

# add dependencies
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
  target_link_libraries(pylimer_tools_bindings ${igraph_LIBRARIES})
	message("Found igraph")
else()
	message(WARNING "DID NOT FIND igraph")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pylimer_tools_bindings ${Boost_LIBRARIES})
endif()

target_include_directories(pylimer_tools_bindings
	PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
