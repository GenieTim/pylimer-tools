cmake_minimum_required(VERSION 3.15.0)
project(pylimer_tests VERSION "0.1.0")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast")
endif()

# load test package
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.0-preview3
)

FetchContent_MakeAvailable(Catch2)

# and other needed packages
find_package(igraph REQUIRED)
if(igraph_FOUND)
  include_directories(${igraph_INCLUDE_DIRS})
	message("Found igraph")
else()
	message(WARNING "DID NOT FIND igraph")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# add sources & targets
add_subdirectory(../src/pylimer_tools_cpp pylimer_tools.out)

add_library(test_sources 
	OBJECT 
		pylimer_tools/test_universe.cpp 
		pylimer_tools/test_universesequence.cpp 
		pylimer_tools/test_molecule.cpp 
)

target_link_libraries(test_sources 
PRIVATE
	Catch2::Catch2
	pylimer_tools
)

add_executable(pylimer_tests pylimer_tools/test.cpp)
target_link_libraries(pylimer_tests 
	PRIVATE 
	Catch2::Catch2WithMain
	pylimer_tools
	igraph
	test_sources
)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(pylimer_tests)
