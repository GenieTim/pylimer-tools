add_library(pylimer_entities)

# add dependencies
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
  target_link_libraries(pylimer_entities ${igraph_LIBRARIES})
	message("Found igraph")
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

# optional additional dependencies
find_package(OpenMP)
if (OpenMP_FOUND)
	add_definitions(-DOPENMP_FOUND=true)
  include_directories(${OpenMP_INCLUDE_DIRS})
	target_link_libraries(pylimer_entities PRIVATE ${OpenMP_LIBRARIES})
	message(INFO "Found OpenMP")
else()
	message("Did not find OpenMP")
endif()
