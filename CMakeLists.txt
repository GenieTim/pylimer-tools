cmake_minimum_required(VERSION 3.15.0)
project(pylimer_tools VERSION "0.1.0")
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
find_package(igraph CONFIG REQUIRED)

add_subdirectory(src/_pylimer_tools)
pybind11_add_module(_pylimer_tools MODULE src/_pylimer_tools/pybind11/entities.cpp)

target_link_libraries(_pylimer_tools 
    PRIVATE 
      igraph
      pylimer_tools
  )
target_compile_definitions(_pylimer_tools PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _pylimer_tools LIBRARY DESTINATION pylimer_tools)
