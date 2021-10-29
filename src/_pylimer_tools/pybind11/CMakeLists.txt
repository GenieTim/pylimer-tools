add_library(pylimer_tools_bindings)

set_target_properties(pylimer_tools_bindings PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(pylimer_tools_bindings
	PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
