cmake_minimum_required(VERSION 3.15.0)

cmake_policy(SET CMP0077 NEW)

project(pylimer_tests VERSION "0.1.15" LANGUAGES CXX C)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_CURRENT_SOURCE_DIR}/../vendor/CMakeLists.txt)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

option(CODE_COVERAGE "Collect coverage from test library" OFF)
option(LEAK_ANALYSIS "Check for memory leaks in the test" OFF)
option(WITH_ERROR "Throw errors when compiling" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pedantic")

if (WITH_ERROR)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wshadow -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -Weffc++") # -Wsign-conversion
else()
 	#-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Weffc++") # -Wsign-conversion
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
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
		pylimer_tools/test_normal_mode_analysis.cpp
		pylimer_tools/test_readers_writers.cpp 
		pylimer_tools/test_universe.cpp 
		pylimer_tools/test_universesequence.cpp 
		pylimer_tools/test_eigen_neighbourlist_benchmark.cpp 
)
# pylimer_tools/test_force_balance2.cpp

IF (WIN32)
		message(AUTHOR_WARNING "On Windows, OpenMP will not be used")
# elseif(CODE_COVERAGE)
# 		message(AUTHOR_WARNING "Code-coverage is enabled, openmp will not be used")
else()
	find_package(OpenMP)
	if (OpenMP_CXX_FOUND)
		include_directories(${OpenMP_CXX_INCLUDE_DIRS})
		target_link_libraries(test_sources PUBLIC OpenMP::OpenMP_CXX)
		add_definitions(-DOPENMP_FOUND=true)
		message(STATUS "Found OpenMP")
	else()
		message(AUTHOR_WARNING "Did not find OpenMP")
	endif()
endif()

add_dependencies(test_sources pylimer_tools)
include_directories(${igraph_INCLUDE_DIRS})
include_directories(${Spectra_INCLUDE_DIRS})
include_directories(${nlopt_INCLUDE_DIRS})
include_directories(${cereal_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS}/Eigen3)
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
)
add_dependencies(header_test_sources igraph::igraph)
add_dependencies(header_test_sources cerealLib)
#add_dependencies(header_test_sources pylimer_tools)
target_link_libraries(header_test_sources PUBLIC 
	Catch2::Catch2
	${igraph_LIBRARIES}
)

target_link_libraries(header_test_sources PUBLIC Eigen3::Eigen)

# separate the header tests out, for the coverage generation
add_executable(header_tests pylimer_tools/test_headers.cpp)

add_dependencies(header_tests igraph::igraph)
add_dependencies(header_tests header_test_sources)

target_link_libraries(header_tests 
	PRIVATE 
	Catch2::Catch2WithMain
	header_test_sources
	${igraph_LIBRARIES}
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

# code coverage setup in the end
if (CODE_COVERAGE)
	set(ENABLE_COVERAGE ON)
	include(CodeCoverage)
	append_coverage_compiler_flags()
	setup_target_for_coverage_lcov(
		NAME pylimer_tests-gcov
		EXECUTABLE pylimer_tests
		DEPENDENCIES pylimer_tests
		LCOV_ARGS --ignore-errors inconsistent,source
		EXECUTABLE_ARGS --skip-benchmarks --durations yes "~[long]"
		EXCLUDE "${CMAKE_CURRENT_LIST_DIR}/build/_deps/*" "${CMAKE_CURRENT_LIST_DIR}/_deps/*" "${CMAKE_CURRENT_LIST_DIR}/../_deps/*"
	)
	setup_target_for_coverage_lcov(
		NAME header_tests-gcov
		EXECUTABLE header_tests
		DEPENDENCIES header_tests
		LCOV_ARGS --ignore-errors inconsistent,source
		EXECUTABLE_ARGS --skip-benchmarks --durations yes "~[long]"
		EXCLUDE "${CMAKE_CURRENT_LIST_DIR}/build/_deps/*" "${CMAKE_CURRENT_LIST_DIR}/_deps/*" "${CMAKE_CURRENT_LIST_DIR}/../_deps/*"
	)
	# BASE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/.."
endif()
