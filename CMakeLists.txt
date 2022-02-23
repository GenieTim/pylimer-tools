cmake_minimum_required(VERSION 3.15.0)

# start project
project(pylimer_tools VERSION "0.1.0")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Ofast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fPIC")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -fPIC")
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/vendor/CMakeLists.txt)

# prepare python bindings
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH /usr/local/anaconda3/lib/python3.7/site-packages/skbuild/resources/cmake/)
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

message("Using toolchain file ${CMAKE_TOOLCHAIN_FILE}")

# add sources & targets
add_subdirectory(src/pylimer_tools_cpp)
add_dependencies(pylimer_tools igraphLib)

pybind11_add_module(pylimer_tools_cpp MODULE 
  src/pylimer_tools_cpp/pybind11/entities.cpp
  src/pylimer_tools_cpp/pybind11/readers.cpp
  src/pylimer_tools_cpp/pybind11/writers.cpp
  src/pylimer_tools_cpp/pybind11/generators.cpp
  src/pylimer_tools_cpp/pybind11/calculations.cpp
  src/pylimer_tools_cpp/pybind11/pylimer_tools_cpp.cpp
)

# add other dependencies
include_directories(${igraph_INCLUDE_DIRS})

message("pylimer_tools_cpp includes igraph from ${igraph_INCLUDE_DIRS}")

# link
target_link_libraries(pylimer_tools_cpp 
  PUBLIC 
    pylimer_tools
    ${igraph_LIBRARIES}
)

target_compile_definitions(pylimer_tools_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS pylimer_tools_cpp LIBRARY DESTINATION .)
