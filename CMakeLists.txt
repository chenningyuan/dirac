cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

project (Dirac C Fortran)

option(DIRAC_ENABLE_PROFILING "foo" OFF)
option(DIRAC_ENABLE_BOUNDS_CHECK "oof" OFF)
option(DIRAC_ENABLE_CODE_COVERAGE "bar" OFF)
#option(DIRAC_ENABLE_XCFUN "Enable the XCFun module" OFF)

find_package(BLAS)
find_package(LAPACK)

if (CMAKE_COMPILER_IS_GNUFortran)
	set (CMAKE_Fortran_FLAGS "-Wall")
	set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g3")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
	if (ENABLE_CODE_COVERAGE)
		set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --coverage")
	endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set (CMAKE_Fortran_FLAGS "-fpp1 -w -assume byterecl")
	set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
	add_definitions(-DVAR_IFC -DSYS_LINUX)
endif ()

if (CMAKE_COMPILER_IS_GNUCC)
	set (CMAKE_CC_FLAGS "-Wall -Wno-sign-compare")
	set (CMAKE_CC_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
	if (ENABLE_CODE_COVERAGE)
		set (CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} --coverage")
	endif()
elseif (CMAKE_CC_COMPILER_ID MATCHES Intel)
	set (CMAKE_CC_FLAGS "")
	set (CMAKE_CC_FLAGS_DEBUG "-O0 -g")
	set (CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG")
endif ()

INCLUDE(build/cmake.definitions)

#if (DIRAC_ENABLE_XCFUN)
#        add_definitions(-DMOD_XCFUN)
#endif()

include_directories(include ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

INCLUDE(build/cmake.sources)

target_link_libraries(dirac.x blas lapack)
