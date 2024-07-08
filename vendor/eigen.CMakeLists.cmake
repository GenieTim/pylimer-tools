# load the Eigen library
if (NOT DEFINED eigen_LOADED)
	find_package(Eigen3 3.4 NO_MODULE) # 3.4
	if(${Eigen3_FOUND}) # AND (${Eigen3_VERSION} VERSION_GREATER_EQUAL 3.4)
			message(STATUS "Found Eigen3 Version: ${Eigen3_VERSION} Path: ${Eigen3_DIR}")
	else()
			include(FetchContent)
			FetchContent_Declare(
				eigen3 
				GIT_REPOSITORY https://gitlab.com/libeigen/eigen
				GIT_TAG 3.4.0
			)
			FetchContent_MakeAvailable(eigen3)
	endif()

	set(eigen_LOADED ON)
endif()

# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake)
find_package(LAPACKE)
# include(FindLAPACKLibs)

# if (LAPACKLIBS_FOUND)
if(LAPACKE_FOUND)
	include_directories(${LAPACKE_INCLUDE_DIRS})
	add_definitions(-DEIGEN_USE_LAPACKE)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		set(BLA_VENDOR Intel10_64lp)
    add_definitions(-DEIGEN_USE_MKL_ALL)
endif()
