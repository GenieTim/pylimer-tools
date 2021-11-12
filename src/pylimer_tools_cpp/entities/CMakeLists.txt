add_library(pylimer_entities)

# add dependencies
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
  target_link_libraries(pylimer_entities ${igraph_LIBRARIES})
  message("Igraph found in entities. ${CMAKE_TOOLCHAIN_FILE}")
else()
	message(WARNING "DID NOT FIND igraph")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pylimer_entities ${Boost_LIBRARIES})
endif()

target_sources(pylimer_entities
	PRIVATE
		Molecule.cpp
		Universe.cpp
		UniverseSequence.cpp
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/Box.h
		${CMAKE_CURRENT_LIST_DIR}/Atom.h
		${CMAKE_CURRENT_LIST_DIR}/Universe.h
		${CMAKE_CURRENT_LIST_DIR}/Molecule.h
		${CMAKE_CURRENT_LIST_DIR}/UniverseSequence.h
)
