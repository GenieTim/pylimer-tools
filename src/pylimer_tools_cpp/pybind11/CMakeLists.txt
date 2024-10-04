add_library(pylimer_tools_bindings)

set_target_properties(pylimer_tools_bindings PROPERTIES LINKER_LANGUAGE CXX)

# add dependencies

include(${CMAKE_CURRENT_LIST_DIR}/../../../vendor/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/../../../vendor/myConfigureFile.cmake)

add_dependencies(pylimer_tools_bindings igraph::igraph)
include_directories(${igraph_INCLUDE_DIRS})
target_link_libraries(pylimer_tools_bindings PUBLIC ${igraph_LIBRARIES})

message("pybind11 includes igraph from ${igraph_INCLUDE_DIRS}")

target_include_directories(pylimer_tools_bindings
	PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
