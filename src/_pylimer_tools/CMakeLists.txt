cmake_minimum_required(VERSION 3.10)

# set the project name
project(pylimer_tools VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add subdirectories
add_subdirectory(utils)
add_subdirectory(entities)

# add libraries
if (DEFINED ADDITIONAL_INCLUDES)
	include_directories(ADDITIONAL_INCLUDES)
endif()

if (DEFINED ADDITIONAL_LIBS)
	link_directories(ADDITIONAL_LIBS)
endif()

# add the executable
add_library(pylimer_tools pylimer_tools.cpp)
set_target_properties(pylimer_tools PROPERTIES VERSION ${PROJECT_VERSION})


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

# add other needed packages
target_link_libraries(pylimer_tools 
	PRIVATE
		pylimer_entities
)
