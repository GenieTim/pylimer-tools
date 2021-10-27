cmake_minimum_required(VERSION 3.10)

# set the project name
project(Test VERSION 1.0)

include_directories(/usr/local/Cellar/igraph/0.9.4/include/)
link_directories(/usr/local/Cellar/igraph/0.9.4/lib)
link_directories(/usr/local/opt/igraph/lib/)

# add the executable
add_executable(Test test.cpp)

# add other needed packages
find_library(igraph SHARED IMPORTED)
target_link_libraries(Test PUBLIC igraph)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
