cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0077 NEW)

# set the project name
project(pylimer_tools VERSION 1.0 LANGUAGES CXX)
add_library(pylimer_tools pylimer_tools.cpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# other options we offer
if (NOT DEFINED CODE_COVERAGE)
	option(CODE_COVERAGE "Collect coverage from test library" OFF)
endif()

if (NOT DEFINED CEREALIZABLE)
	option(CEREALIZABLE "Enable serialisation of various classes" ON)
endif()
if (CEREALIZABLE)
	add_definitions(-DCEREALIZABLE)
	add_compile_definitions(CEREALIZABLE)
endif()

# add dependencies
include(${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/CMakeLists.txt)

include(${CMAKE_CURRENT_LIST_DIR}/../../vendor/myConfigureFile.cmake)

add_dependencies(pylimer_tools igraph)
include_directories(${igraph_INCLUDE_DIRS})
add_dependencies(pylimer_tools nlopt)
include_directories(${nlopt_INCLUDE_DIRS})
add_dependencies(pylimer_tools cerealLib)
include_directories(${cereal_INCLUDE_DIRS})
add_dependencies(pylimer_tools Spectra)
include_directories(${Spectra_INCLUDE_DIRS})

target_link_libraries(pylimer_tools PUBLIC 
	${igraph_LIBRARIES} 
	${nlopt_LIBRARIES}
)

# add_dependencies(pylimer_tools Eigen3::Eigen)

target_link_libraries(pylimer_tools PUBLIC Eigen3::Eigen)

if (LAPACKE_FOUND)
	target_link_libraries(pylimer_tools PUBLIC ${LAPACKE_LIBRARIES})
endif()

message(STATUS "pylimer_tools includes igraph from ${igraph_INCLUDE_DIRS}, nlopt from ${nlopt_INCLUDE_DIRS}, Spectra from ${Spectra_INCLUDE_DIRS} and cereal from ${cereal_INCLUDE_DIRS}")

# optional additional dependencies
IF (WIN32)
		message(AUTHOR_WARNING "On Windows, OpenMP will not be used")
elseif(CODE_COVERAGE)
		message(AUTHOR_WARNING "Code-coverage is enabled, openmp will not be used")
else()
	find_package(OpenMP)
	if (OpenMP_CXX_FOUND)
		include_directories(${OpenMP_CXX_INCLUDE_DIRS})
		# target_link_libraries(pylimer_tools PRIVATE ${OpenMP_CXX_LIBRARIES})
		target_link_libraries(pylimer_tools PUBLIC OpenMP::OpenMP_CXX)
		add_definitions(-DOPENMP_FOUND=true)
		message(STATUS "Found OpenMP")
	else()
		message(AUTHOR_WARNING "Did not find OpenMP")
	endif()
endif()

target_sources(pylimer_tools
	PRIVATE
		${CMAKE_CURRENT_LIST_DIR}/entities/AtomGraphParent.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/Molecule.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/Universe.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/UniverseSequence.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/NeighbourList.cpp
		${CMAKE_CURRENT_LIST_DIR}/entities/EigenNeighbourList.cpp

		${CMAKE_CURRENT_LIST_DIR}/io/AveFileReader.cpp
		${CMAKE_CURRENT_LIST_DIR}/io/DataFileParser.cpp
		${CMAKE_CURRENT_LIST_DIR}/io/DumpFileParser.cpp

		${CMAKE_CURRENT_LIST_DIR}/calc/Correlator.cpp
		${CMAKE_CURRENT_LIST_DIR}/calc/NormalModeAnalyzer.cpp
		
		${CMAKE_CURRENT_LIST_DIR}/topo/EntanglementDetector.cpp

		${CMAKE_CURRENT_LIST_DIR}/sim/OutputSupportingSimulation.cpp
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceEvaluator.cpp
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceRelaxation.cpp
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceBalance.cpp
		${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp
	PUBLIC
		${CMAKE_CURRENT_LIST_DIR}/entities/Box.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Atom.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Molecule.h
		${CMAKE_CURRENT_LIST_DIR}/entities/Universe.h
		${CMAKE_CURRENT_LIST_DIR}/entities/UniverseSequence.h
		${CMAKE_CURRENT_LIST_DIR}/entities/NeighbourList.h
		${CMAKE_CURRENT_LIST_DIR}/entities/EigenNeighbourList.h

		${CMAKE_CURRENT_LIST_DIR}/io/AveFileReader.h
		${CMAKE_CURRENT_LIST_DIR}/io/DataFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/io/DumpFileParser.h
		${CMAKE_CURRENT_LIST_DIR}/utils/CerealUtils.h
		${CMAKE_CURRENT_LIST_DIR}/utils/ExtraEigenTypes.h
		${CMAKE_CURRENT_LIST_DIR}/utils/GraphUtils.h
		${CMAKE_CURRENT_LIST_DIR}/utils/PerformanceTimer.h
		${CMAKE_CURRENT_LIST_DIR}/utils/VectorUtils.h

		${CMAKE_CURRENT_LIST_DIR}/calc/Correlator.h
		${CMAKE_CURRENT_LIST_DIR}/calc/MEHPanalysis.h
		${CMAKE_CURRENT_LIST_DIR}/calc/MMTanalysis.h
		${CMAKE_CURRENT_LIST_DIR}/calc/NormalModeAnalyzer.h

		${CMAKE_CURRENT_LIST_DIR}/topo/EntanglementDetector.h

		${CMAKE_CURRENT_LIST_DIR}/sim/OutputSupportingSimulation.h
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPUtilityStructures.h
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceEvaluator.h
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceRelaxation.h
		${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceBalance.h
		${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.h
)
# ${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceBalance2.cpp
# ${CMAKE_CURRENT_LIST_DIR}/sim/MEHPForceBalance2.h

# DEBUG whether vectorization works or not
# if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#   # using Clang
# 	# set_source_files_properties( ${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp PROPERTIES COMPILE_FLAGS "-O3 -Rpass-missed -Rpass-analysis" )
# 	set_source_files_properties( ${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp PROPERTIES COMPILE_FLAGS "-O3 -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize" )
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#   # using GCC
# 	set_source_files_properties( ${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp PROPERTIES COMPILE_FLAGS "-O3 -fopt-info-vec-missed" )
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#   # using Intel C++
# 	set_source_files_properties( ${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp PROPERTIES COMPILE_FLAGS "-O3 -qopt-report" )
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#   # using Visual Studio C++
# 	set_source_files_properties( ${CMAKE_CURRENT_LIST_DIR}/sim/DPDSimulator.cpp PROPERTIES COMPILE_FLAGS "-O3 -Qvec-report:2" )
# endif()

set_target_properties(pylimer_tools PROPERTIES VERSION ${PROJECT_VERSION})

if (CODE_COVERAGE)
	set(ENABLE_COVERAGE ON)
	if (CMAKE_BUILD_TYPE EQUAL "RELEASE")
		message(WARNING "Using coverage reduces build type to debug from release")
	endif()
	include(CodeCoverage)
	append_coverage_compiler_flags()
else()
	message("Code Coverage is disabled.")
endif()

if (LEAK_ANALYSIS)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
	set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
	target_link_libraries(pylimer_tools 
		PRIVATE
		-fsanitize=address
	)
endif()

OUTPUT_FLAGS(pylimer_tools)

# add other needed packages
# target_link_libraries(pylimer_tools 
# 	PRIVATE
# )
# target_link_libraries(pylimer_tools ${igraph_LIBRARIES})
