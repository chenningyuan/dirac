cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

project(Dirac C Fortran)

option(WITH_PROFILING     "Enable profiling"          OFF)
option(WITH_BOUNDS_CHECK  "Enable bounds check"       OFF)
option(WITH_CODE_COVERAGE "Enable code coverage"      OFF)
option(WITH_EXTERNAL_MATH "Use external math library" ON)
option(WITH_MPI           "Use MPI parallelization"   OFF)

add_definitions(-DPRG_DIRAC)

# now only Linux, later add more machine info
add_definitions(-DSYS_LINUX)

# later get from user with a good default
add_definitions(-DINSTALL_WRKMEM=64000000)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
#   set(CMAKE_Fortran_FLAGS         "-Wall")
    set(CMAKE_Fortran_FLAGS         "     ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    add_definitions(-DVAR_GFORTRAN)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    add_definitions(-DVAR_IFC)
endif ()

if (CMAKE_COMPILER_IS_GNUCC)
        set(CMAKE_CC_FLAGS         "-Wall -Wno-sign-compare")
        set(CMAKE_CC_FLAGS_DEBUG   "-O0 -g3")
        set(CMAKE_CC_FLAGS_RELEASE "-O3 -Wno-unused")
        if (WITH_CODE_COVERAGE)
                set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} --coverage")
        endif ()
elseif (CMAKE_CC_COMPILER_ID MATCHES Intel)
       set(CMAKE_CC_FLAGS         "")
       set(CMAKE_CC_FLAGS_DEBUG   "-O0 -g")
       set(CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG")
endif ()

if (WITH_EXTERNAL_MATH)
        find_package(BLAS)
        find_package(LAPACK)
endif ()

if (WITH_MPI)
        message("mpi active")
        find_package(MPI)
        add_definitions(-DVAR_MPI)
        if (MPI_FOUND)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
                include_directories(${MPI_INCLUDE_PATH})
        endif ()
endif ()

include_directories(include ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include(build/cmake.sources)
include(build/cmake.definitions)

if (WITH_EXTERNAL_MATH)
        if (WITH_MPI)
                target_link_libraries(dirac.x blas lapack)
        else ()
                target_link_libraries(dirac.x blas lapack mpi)
        endif ()
endif ()
