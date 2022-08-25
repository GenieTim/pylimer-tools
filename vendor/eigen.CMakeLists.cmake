include(FetchContent)

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
