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
include_directories(/usr/local/Cellar/igraph/0.9.4/include/)
link_directories(/usr/local/Cellar/igraph/0.9.4/lib)
link_directories(/usr/local/opt/igraph/lib/)

# add the executable
add_library(pylimer_tools pylimer_tools.cpp)
set_target_properties(pylimer_tools PROPERTIES VERSION ${PROJECT_VERSION})

# add other needed packages
find_library(igraph SHARED IMPORTED)
target_link_libraries(pylimer_tools 
	PRIVATE
		pylimer_entities
	PUBLIC 
		igraph
)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pylimer_tools ${Boost_LIBRARIES})
endif()
