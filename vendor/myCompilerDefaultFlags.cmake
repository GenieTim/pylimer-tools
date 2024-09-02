# improve a few defaults for the compiler flags

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

if (NOT DEFINED HIGH_PERFORMANCE)
	option(HIGH_PERFORMANCE "Add even more compiler flags to improve performance" OFF)
endif()

if (HIGH_PERFORMANCE)
	if (MSVC)
		add_compile_options(/03)
	else()
		add_compile_options(-O3 -march=native)
	endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fPIC")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -fPIC")
endif()

function(OUTPUT_FLAGS target_name)
get_target_property(COMPILE_DEFS ${target_name} COMPILE_DEFINITIONS)
message(STATUS "Compile definitions of ${target_name} are ${COMPILE_DEFS} for build type ${CMAKE_BUILD_TYPE}")
endfunction()
