cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(xcint Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

option(ENABLE_PROFILING      "Enable profiling"          OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"       OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"      OFF)
option(ENABLE_INTERNAL_MATH  "Enable internal math library"       OFF)
option(ENABLE_MPI            "Use MPI parallelization"   OFF)

include(cmake/Archs.cmake)
include(FCompilers)
include(CCompilers)

set(LIBS)

add_definitions(-DMOD_XCFUN)
add_subdirectory(xcfun)
set(LIBS
    ${LIBS}
    xcfun
    )

include(Math)

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        add_definitions(-DVAR_MPI2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

include_directories(
    ${PROJECT_BINARY_DIR}/modules
    ${PROJECT_BINARY_DIR}/xcfun/fortran
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

include(XCint_Sources)

add_library(
    xcint
    ${SOURCES}
    )

add_executable(
    xcint.x
    src/test_driver.F90
    )

target_link_libraries(
    xcint.x
    xcint
    ${LIBS}
    )

install(TARGETS xcint ARCHIVE DESTINATION lib)

include(Tests)
include(CTest)
enable_testing()

include(CompilationInfo)
include(Documentation)
