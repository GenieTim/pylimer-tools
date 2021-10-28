add_library(pylimer_entities)
# add libraries
include_directories(/usr/local/Cellar/igraph/0.9.4/include/)
link_directories(/usr/local/Cellar/igraph/0.9.4/lib)
link_directories(/usr/local/opt/igraph/lib/)

include_directories(/usr/local/Cellar/boost/1.76.0/include)

target_sources(pylimer_entities
	PRIVATE
		UniverseSequence.cpp
		Molecule.cpp
		Universe.cpp
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/Box.h
		${CMAKE_CURRENT_LIST_DIR}/Atom.h
		${CMAKE_CURRENT_LIST_DIR}/Universe.h
		${CMAKE_CURRENT_LIST_DIR}/Molecule.h
		${CMAKE_CURRENT_LIST_DIR}/UniverseSequence.h
)

target_include_directories(pylimer_entities
	PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
