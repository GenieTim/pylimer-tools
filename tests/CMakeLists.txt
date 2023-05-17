cmake_minimum_required(VERSION 3.15.0)

cmake_policy(SET CMP0077 NEW)

project(pylimer_tests VERSION "0.1.15" LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_CURRENT_SOURCE_DIR}/../vendor/CMakeLists.txt)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

option(CODE_COVERAGE "Collect coverage from test library" OFF)
option(LEAK_ANALYSIS "Check for memory leaks in the test" OFF)

if (CODE_COVERAGE)
	set(ENABLE_COVERAGE ON)
	set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../vendor/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})
	find_package(codecov)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -Weffc++") # -Wsign-conversion

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
# memory leak analysis:
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=leak")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=leak")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g")
endif()

# load test package
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.3.2
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

add_library(test_sources 
	OBJECT 
		pylimer_tools/test_dpd_simulator.cpp
		pylimer_tools/test_force_balance.cpp
		pylimer_tools/test_force_relaxation.cpp
		pylimer_tools/test_generator.cpp 
		pylimer_tools/test_molecule.cpp 
		pylimer_tools/test_neighbourlist.cpp 
		pylimer_tools/test_readers_writers.cpp 
		pylimer_tools/test_universe.cpp 
		pylimer_tools/test_universesequence.cpp 
)
# pylimer_tools/test_force_balance2.cpp

add_dependencies(test_sources pylimer_tools)
include_directories(${igraph_INCLUDE_DIRS})
include_directories(${nlopt_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS}/eigen3)
message("Found igraph for tests. Including ${igraph_INCLUDE_DIRS}, linking to ${igraph_LIBRARIES}")

target_link_libraries(test_sources 
PRIVATE
	Catch2::Catch2
	pylimer_tools
	${igraph_LIBRARIES}
	${nlopt_LIBRARIES}
)

add_executable(pylimer_tests pylimer_tools/test.cpp)

add_dependencies(pylimer_tests pylimer_tools)
add_dependencies(pylimer_tests test_sources)

target_link_libraries(pylimer_tests 
	PRIVATE 
	Catch2::Catch2WithMain
	test_sources
	pylimer_tools
	${igraph_LIBRARIES}
	${nlopt_LIBRARIES}
)

add_library(header_test_sources
	OBJECT
		pylimer_tools/test_calc.cpp
		pylimer_tools/test_vector_utils.cpp
		pylimer_tools/test_eigen_neighbourlist_benchmark.cpp 
)
add_dependencies(header_test_sources pylimer_tools)
target_link_libraries(header_test_sources PUBLIC 
	Catch2::Catch2
)

# link Eigen libraries. Difficulty due to Windows GitHub Action, somehow.
if (TARGET Eigen3::Eigen)
	message(STATUS "Target Eigen3::Eigen found")
	target_link_libraries(header_test_sources PUBLIC Eigen3::Eigen)
else()
	include_directories(${EIGEN3_INCLUDE_DIRS})
	include_directories(${EIGEN3_INCLUDE_DIRS}/eigen3)
	target_link_libraries(header_test_sources PUBLIC ${EIGEN3_LIBRARIES})
	message(STATUS "Eigen include directories: ${EIGEN3_INCLUDE_DIRS}, libraries ${EIGEN3_LIBRARIES}")
endif()

# separate the header tests out, for the coverage generation
add_executable(header_tests pylimer_tools/test_headers.cpp)

add_dependencies(header_tests igraphLib)
add_dependencies(header_tests header_test_sources)

target_link_libraries(header_tests 
	PRIVATE 
	Catch2::Catch2WithMain
	header_test_sources
	${igraph_LIBRARIES}
	${EIGEN3_LIBRARIES}
)

if (LEAK_ANALYSIS)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
	set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
	target_link_libraries(test_sources 
		PRIVATE
		-fsanitize=address
	)
	target_link_libraries(pylimer_tests 
		PRIVATE
		-fsanitize=address
	)
	target_link_libraries(header_tests 
		PRIVATE
		-fsanitize=address
	)
endif()

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(pylimer_tests)

catch_discover_tests(header_tests)

if(CODE_COVERAGE)
	
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
	set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
	target_compile_options(pylimer_tests PRIVATE -O0 -g -fprofile-arcs -ftest-coverage )

	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-abs-path")
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fprofile-abs-path")
		target_compile_options(pylimer_tests PRIVATE -fprofile-abs-path)
	endif()

	add_coverage(pylimer_tests)
	target_link_libraries(
		pylimer_tests
		PRIVATE
		--coverage
	)
	
	# add_coverage(pylimer_tools)
	add_coverage(test_sources)

	target_link_libraries(
		test_sources
		PRIVATE
		--coverage
	)

	target_compile_options(header_tests PRIVATE -O0 -g -fprofile-arcs -ftest-coverage )
	if(CMAKE_COMPILER_IS_GNUCXX)
		target_compile_options(header_tests PRIVATE -fprofile-abs-path)
	endif()

	add_coverage(header_tests)
	target_link_libraries(header_tests 
		PRIVATE 
		--coverage
	)

	coverage_evaluate()
endif()
