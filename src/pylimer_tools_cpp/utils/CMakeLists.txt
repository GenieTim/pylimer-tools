add_library(pylimer_utils INTERFACE)

set(CMAKE_CXX_STANDARD 20)

# add libraries
include_directories(/usr/local/Cellar/igraph/0.9.4/include/)
link_directories(/usr/local/Cellar/igraph/0.9.4/lib)
link_directories(/usr/local/opt/igraph/lib/)

include_directories(/usr/local/Cellar/boost/1.76.0/include)

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
