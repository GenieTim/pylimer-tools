cmake_minimum_required(VERSION 3.15.0)
project(pylimer_tests VERSION "0.1.15" LANGUAGES CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../vendor/CMakeLists.txt)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

option(CODE_COVERAGE "Collect coverage from test library" OFF)
option(LEAK_ANALYSIS "Check for memory leaks in the test" OFF)

if (CODE_COVERAGE)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../vendor/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})
set(ENABLE_COVERAGE ON)
find_package(codecov)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
# memory leak analysis:
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=leak")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=leak")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Ofast")
endif()

# load test package
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.0-preview4
)

FetchContent_MakeAvailable(Catch2)

# and other needed packages
# setupIgraph()
# find_package(igraph REQUIRED)
# if(igraph_FOUND)
# else()
# 	message(WARNING "DID NOT FIND igraph")
# endif()

# add sources & targets
add_subdirectory(../src/pylimer_tools_cpp pylimer_tools.out)

add_dependencies(pylimer_tools igraphLib)
include_directories(${igraph_INCLUDE_DIRS})
message("Found igraph for tests. Including ${igraph_INCLUDE_DIRS}, linking to ${igraph_LIBRARIES}")

add_library(test_sources 
	OBJECT 
		pylimer_tools/test_generator.cpp 
		pylimer_tools/test_force_relaxation.cpp
		pylimer_tools/test_universe.cpp 
		pylimer_tools/test_universesequence.cpp 
		pylimer_tools/test_molecule.cpp 
		pylimer_tools/test_atom.cpp 
)

target_link_libraries(test_sources 
PRIVATE
	Catch2::Catch2
	pylimer_tools
	${igraph_LIBRARIES}
)

add_executable(pylimer_tests pylimer_tools/test.cpp)
target_link_libraries(pylimer_tests 
	PRIVATE 
	Catch2::Catch2WithMain
	test_sources
	pylimer_tools
	${igraph_LIBRARIES}
)

if (LEAK_ANALYSIS)
	if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
		target_link_libraries(pylimer_tests 
		PRIVATE
		-fsanitize=address
		)
	else()
		message(WARNING "Leak analysis only supported for Clang")
	endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(pylimer_tests)

if(CODE_COVERAGE)
	
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline")
	set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline")

	target_link_libraries(
		pylimer_tests
		PRIVATE
		--coverage
	)
	
	# add_coverage(pylimer_tools)
	add_coverage(test_sources)
	add_coverage(pylimer_tests)
endif()
