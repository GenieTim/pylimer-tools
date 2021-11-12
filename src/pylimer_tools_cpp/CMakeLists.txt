cmake_minimum_required(VERSION 3.10)

# set the project name
project(pylimer_tools VERSION 1.0)
add_library(pylimer_tools pylimer_tools.cpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add libraries
if (DEFINED ADDITIONAL_INCLUDES)
	include_directories(ADDITIONAL_INCLUDES)
endif()

if (DEFINED ADDITIONAL_LIBS)
	link_directories(ADDITIONAL_LIBS)
endif()

# add dependencies
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
  target_link_libraries(pylimer_tools ${igraph_LIBRARIES})
	message("Found igraph")
else()
	message(WARNING "DID NOT FIND igraph")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pylimer_tools ${Boost_LIBRARIES})
endif()

# optional additional dependencies
find_package(OpenMP)
if (OpenMP_FOUND)
	add_definitions(-DOPENMP_FOUND=true)
  include_directories(${OpenMP_INCLUDE_DIRS})
	target_link_libraries(pylimer_tools PRIVATE ${OpenMP_LIBRARIES})
	message(INFO "Found OpenMP")
else()
	message("Did not find OpenMP")
endif()

target_sources(pylimer_tools
	PRIVATE
		${CMAKE_CURRENT_LIST_DIR}/entities/Molecule.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/Universe.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/UniverseSequence.cpp

		${CMAKE_CURRENT_LIST_DIR}/utils/DataFileParser.cpp
		${CMAKE_CURRENT_LIST_DIR}/utils/DumpFileParser.cpp
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/entities/Box.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Atom.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Universe.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Molecule.h
		${CMAKE_CURRENT_LIST_DIR}/entities/UniverseSequence.h

		${CMAKE_CURRENT_LIST_DIR}/utils/DataFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/utils/DumpFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/utils/GraphUtils.h
		${CMAKE_CURRENT_LIST_DIR}/utils/VectorUtils.h

		${CMAKE_CURRENT_LIST_DIR}/calc/MEHPanalysis.h
		${CMAKE_CURRENT_LIST_DIR}/calc/MMTanalysis.h
)

set_target_properties(pylimer_tools PROPERTIES VERSION ${PROJECT_VERSION})

# add other needed packages
# target_link_libraries(pylimer_tools 
# 	PRIVATE
# )
