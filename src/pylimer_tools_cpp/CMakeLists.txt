cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0077 NEW)

# set the project name
project(pylimer_tools VERSION 1.0)
add_library(pylimer_tools pylimer_tools.cpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# TODO: remove
# set(CMAKE_CXX_FLAGS "-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -g")

# add libraries
if (DEFINED ADDITIONAL_INCLUDES)
	include_directories(ADDITIONAL_INCLUDES)
endif()

if (DEFINED ADDITIONAL_LIBS)
	link_directories(ADDITIONAL_LIBS)
endif()

# add dependencies
include(${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/CMakeLists.txt)

add_dependencies(pylimer_tools igraph)
include_directories(${igraph_INCLUDE_DIRS})
target_link_libraries(pylimer_tools PUBLIC ${igraph_LIBRARIES})

message("pylimer_tools includes igraph from ${igraph_INCLUDE_DIRS}")

# optional additional dependencies
find_package(OpenMP)
if (OpenMP_FOUND AND OpenMP_INCLUDE_DIRS)
  include_directories(${OpenMP_INCLUDE_DIRS})
	target_link_libraries(pylimer_tools PRIVATE ${OpenMP_LIBRARIES})
	add_definitions(-DOPENMP_FOUND=true)
	message(INFO "Found OpenMP")
else()
	message("Did not find OpenMP")
endif()

target_sources(pylimer_tools
	PRIVATE
		${CMAKE_CURRENT_LIST_DIR}/entities/AtomGraphParent.cpp
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

if (NOT DEFINED CODE_COVERAGE)
option(CODE_COVERAGE "Collect coverage from test library" OFF)
endif()

if(CODE_COVERAGE)
	if (CMAKE_BUILD_TYPE EQUAL "RELEASE")
		message(WARNING "Using coverage reduces build type to debug from release")
	endif()
	set(ENABLE_COVERAGE ON)
	set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})
	find_package(codecov)

	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-inline-functions -fno-default-inline")
	set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-inline-functions -fno-default-inline")
	target_compile_options(pylimer_tools PRIVATE -O0 -g -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-inline-functions -fno-default-inline)

	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-abs-path")
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fprofile-abs-path")
		target_compile_options(pylimer_tools PRIVATE -fprofile-abs-path)
	endif()

	target_link_libraries(
		pylimer_tools
		PRIVATE
		--coverage
	)
	add_coverage(pylimer_tools)
else()
	message("Code Coverage is disabled.")
endif()

# add other needed packages
# target_link_libraries(pylimer_tools 
# 	PRIVATE
# )
# target_link_libraries(pylimer_tools ${igraph_LIBRARIES})
