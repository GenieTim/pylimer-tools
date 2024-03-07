add_library(pylimer_tools_bindings)

set_target_properties(pylimer_tools_bindings PROPERTIES LINKER_LANGUAGE CXX)

# add dependencies

include(${CMAKE_CURRENT_SOURCE_DIR}/../../../vendor/CMakeLists.txt)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../version_config.h.in ${CMAKE_BINARY_DIR}/generated/version_config.h )
include_directories( ${CMAKE_BINARY_DIR}/generated/ )

add_dependencies(pylimer_tools_bindings igraphLib)
include_directories(${igraph_INCLUDE_DIRS})
target_link_libraries(pylimer_tools_bindings PUBLIC ${igraph_LIBRARIES})

message("pybind11 includes igraph from ${igraph_INCLUDE_DIRS}")

target_include_directories(pylimer_tools_bindings
	PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
