cmake_minimum_required(VERSION 3.15.0)
project(pylimer_tools VERSION "0.1.0")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast")

# prepare python bindings
find_package(PythonExtensions REQUIRED)

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific directory
  # here.
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

# Now we can find pybind11
find_package(pybind11 CONFIG REQUIRED)

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
add_subdirectory(src/_pylimer_tools)
pybind11_add_module(pylimer_bound_entities MODULE src/_pylimer_tools/pybind11/entities.cpp)
pybind11_add_module(pylimer_bound_calc MODULE src/_pylimer_tools/pybind11/calculations.cpp)

target_link_libraries(pylimer_bound_entities 
    PRIVATE 
      igraph
      pylimer_tools
  )
target_compile_definitions(pylimer_bound_entities PRIVATE VERSION_INFO=${PROJECT_VERSION})

target_link_libraries(pylimer_bound_calc 
    PRIVATE 
      igraph
      pylimer_tools
  )
target_compile_definitions(pylimer_bound_calc PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS pylimer_bound_entities LIBRARY DESTINATION .)
install(TARGETS pylimer_bound_calc LIBRARY DESTINATION .)
