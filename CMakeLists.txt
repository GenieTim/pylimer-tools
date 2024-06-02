cmake_minimum_required(VERSION 3.15.0)

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

# start project
if (NOT DEFINED VERSION_NR)
project(pylimer_tools VERSION "0.2.4")
else()
project(pylimer_tools VERSION ${VERSION_NR})
endif()
add_definitions( -DOVERALL_PROJECT_VERSION="${CMAKE_PROJECT_VERSION}")

include(${CMAKE_CURRENT_LIST_DIR}/vendor/myConfigureFile.cmake)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fPIC")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -fPIC")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/vendor/CMakeLists.txt)

# prepare python bindings
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/vendor/pybind11)

message("Using toolchain file ${CMAKE_TOOLCHAIN_FILE}")

# add sources & targets
add_subdirectory(src/pylimer_tools_cpp)
add_dependencies(pylimer_tools igraph)
add_dependencies(pylimer_tools nlopt)
add_dependencies(pylimer_tools cerealLib)

pybind11_add_module(pylimer_tools_cpp MODULE 
  src/pylimer_tools_cpp/pybind11/readers.cpp
  src/pylimer_tools_cpp/pybind11/entities.cpp
  src/pylimer_tools_cpp/pybind11/writers.cpp
  src/pylimer_tools_cpp/pybind11/generators.cpp
  src/pylimer_tools_cpp/pybind11/calculations.cpp
  src/pylimer_tools_cpp/pybind11/pylimer_tools_cpp.cpp
)

# add other dependencies
include_directories(${igraph_INCLUDE_DIRS})
include_directories(${nlopt_INCLUDE_DIRS})
include_directories(${cereal_INCLUDE_DIRS})

message("pylimer_tools_cpp includes igraph from ${igraph_INCLUDE_DIRS}, nlopt from ${nlopt_INCLUDE_DIRS} and cereal from ${cereal_INCLUDE_DIRS}")

# link
target_link_libraries(pylimer_tools_cpp 
  PUBLIC 
    pylimer_tools
    ${igraph_LIBRARIES}
    ${nlopt_LIBRARIES}
)

if(CMAKE_COMPILER_IS_GNUCXX)
target_link_libraries(pylimer_tools_cpp 
  PUBLIC 
    stdc++fs
)
endif()

# link Eigen libraries. Difficulty due to Windows GitHub Action, somehow.
if (TARGET Eigen3::Eigen)
	message(STATUS "Target Eigen3::Eigen found")
	target_link_libraries(pylimer_tools PUBLIC Eigen3::Eigen)
else()
	include_directories(${EIGEN3_INCLUDE_DIRS})
	include_directories(${EIGEN3_INCLUDE_DIRS}/eigen3)
	target_link_libraries(pylimer_tools PUBLIC ${EIGEN3_LIBRARIES})
	message(STATUS "Eigen include directories: ${EIGEN3_INCLUDE_DIRS}, libraries ${EIGEN3_LIBRARIES}")
endif()

# compile & install
target_compile_definitions(pylimer_tools_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS pylimer_tools_cpp LIBRARY DESTINATION .)
