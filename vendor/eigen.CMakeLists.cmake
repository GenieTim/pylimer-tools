# load the Eigen library
if (NOT DEFINED eigen_LOADED)
	find_package(Eigen3 3.4 NO_MODULE) # 3.4
	if(TARGET Eigen3::Eigen)  #(${Eigen3_FOUND}) # AND (${Eigen3_VERSION} VERSION_GREATER_EQUAL 3.4)
			message(STATUS "Found Eigen3 Version: ${Eigen3_VERSION} Path: ${Eigen3_DIR}")
	else()
			include(FetchContent)
			FetchContent_Declare(
				Eigen3 
				GIT_REPOSITORY https://gitlab.com/libeigen/eigen
				GIT_TAG 3.4.0
				GIT_SHALLOW TRUE
				GIT_PROGRESS TRUE
			)
			set(EIGEN_BUILD_DOC OFF)
			set(EIGEN_BUILD_PKGCONFIG OFF)
			FetchContent_MakeAvailable(Eigen3)
	endif()

	set(eigen_LOADED ON)
	message(STATUS "Eigen include directories: ${EIGEN3_INCLUDE_DIRS}, libraries ${EIGEN3_LIBRARIES}")
endif()

# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake)
# include(${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake)
find_package(LAPACKE)
# include(FindLAPACKLibs)

# if (LAPACKLIBS_FOUND)
if(LAPACKE_FOUND)
	include_directories(${LAPACKE_INCLUDE_DIRS})
	include_directories(${CMAKE_CURRENT_LIST_DIR}/lapacke-extra)
	add_definitions(-DEIGEN_USE_LAPACKE)
	add_definitions(-DHAVE_LAPACK_CONFIG_H)
	add_definitions(-DLAPACK_COMPLEX_CPP)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
		set(BLA_VENDOR Intel10_64lp)
    add_definitions(-DEIGEN_USE_MKL_ALL)
endif()
