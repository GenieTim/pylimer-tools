include( ExternalProject )
add_library(pylimer_entities)

# add dependencies

find_package(igraph)
if(NOT igraph_FOUND)
	message(WARNING "DID NOT FIND igraph. Compiling as external project...")
	ExternalProject_Add (igraph
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/igraph
    GIT_REPOSITORY git@github.com:igraph/igraph.git
    BUILD_IN_SOURCE TRUE
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
		BUILD_COMMAND ${CMAKE_COMMAND} --build
    GIT_PROGRESS 1
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
 	)
	SET(igraph_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/igraph)
	SET(igraph_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/igraph)
endif()
target_link_libraries(pylimer_entities ${igraph_LIBRARIES})
include_directories(${igraph_INCLUDE_DIRS})

# BOOST, vgl. https://gitlab.kitware.com/cmake/cmake/-/issues/18744
if(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    set(BOOST_BUILD_COMMAND ${CMAKE_COMMAND} -E time <SOURCE_DIR>/b2 --prefix=<INSTALL_DIR> --build-dir=${BOOST_BUILD_DIR} -j $ENV{CMAKE_BUILD_PARALLEL_LEVEL} variant=release link=static cxxflags=-fPIC cxxstd=14 install)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL Android)
    set(BOOST_BUILD_COMMAND ${CMAKE_COMMAND} -E time <SOURCE_DIR>/b2 --without-python --abbreviate-paths -d2 --debug-configuration --prefix=<INSTALL_DIR> --build-dir=${BOOST_BUILD_DIR} -j $ENV{CMAKE_BUILD_PARALLEL_LEVEL} architecture=arm address-model=32 binary-format=elf abi=aapcs target-os=android variant=release link=static cxxflags=-fPIC cxxstd=14 install)
else()
    message(WARNING "Building Boost on ${CMAKE_SYSTEM_NAME} is not currently supported.")
    set(BOOST_BUILD_COMMAND ${CMAKE_COMMAND} -E time <SOURCE_DIR>/b2 --prefix=<INSTALL_DIR> --build-dir=${BOOST_BUILD_DIR} -j $ENV{CMAKE_BUILD_PARALLEL_LEVEL} variant=release link=static cxxflags=-fPIC cxxstd=17 install)
endif()

find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
	message(WARNING "DID NOT FIND boost. Compiling as external project...")
	ExternalProject_Add (boost
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/boost
    GIT_REPOSITORY git@github.com:boostorg/boost.git
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E time  ${CMAKE_CURRENT_BINARY_DIR}/boost/bootstrap.sh --with-libraries=filesystem,program_options,regex
              COMMAND ${CMAKE_COMMAND} -E make_directory ${BOOST_BUILD_DIR}
    BUILD_COMMAND ${BOOST_BUILD_COMMAND}
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
 	)
	SET(Boost_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/boost)
endif()
target_link_libraries(pylimer_entities ${Boost_LIBRARIES})

target_sources(pylimer_entities
	PRIVATE
		Molecule.cpp
		Universe.cpp
		UniverseSequence.cpp
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/Box.h
		${CMAKE_CURRENT_LIST_DIR}/Atom.h
		${CMAKE_CURRENT_LIST_DIR}/Universe.h
		${CMAKE_CURRENT_LIST_DIR}/Molecule.h
		${CMAKE_CURRENT_LIST_DIR}/UniverseSequence.h
)
