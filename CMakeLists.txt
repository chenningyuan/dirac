cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

project(Dirac Fortran C CXX)

option(WITH_PROFILING     "Enable profiling"          OFF)
option(WITH_BOUNDS_CHECK  "Enable bounds check"       OFF)
option(WITH_CODE_COVERAGE "Enable code coverage"      OFF)
option(WITH_EXTERNAL_MATH "Use external math library" ON)
option(WITH_MPI           "Use MPI parallelization"   OFF)
option(WITH_XCFUN         "Use XCFun library"         OFF)

add_definitions(-DPRG_DIRAC)

set(DIRAC_WORK_MEM_WORDS "64000000" CACHE STRING "DIRAC work memory in words")
add_definitions(-DINSTALL_WRKMEM=${DIRAC_WORK_MEM_WORDS})

# radovan: note to myself for later
# mark_as_advanced(VALGRIND_EXECUTABLE)

# now only Linux, later add more machine info
add_definitions(-DSYS_LINUX)

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

if (CMAKE_COMPILER_IS_GNUCXX)
        set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings")
        set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
        set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
        if (WITH_CODE_COVERAGE)
                set (CMAKE_CXX_FLAGS 
                        "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
                set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
        endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
        set (CMAKE_CXX_FLAGS "-Wno-unknown-pragmas")
        set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
        set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
        set (CMAKE_CXX_LINK_FLAGS "-shared-intel")
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

find_package(Git)
include(ExternalProject)

if (WITH_XCFUN)
        ExternalProject_Add(xcfun
                PREFIX ${CMAKE_SOURCE_DIR}/import
                GIT_REPOSITORY git://repo.ctcc.no/xcfun.git
                GIT_TAG master
                INSTALL_DIR ${CMAKE_BINARY_DIR}/import
                CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/import
                )
endif ()


include_directories(include ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_custom_command (
        OUTPUT
        ${CMAKE_SOURCE_DIR}/include/fortran_interface.h
        DEPENDS
        gen_fortran_interface_h
        COMMAND
        gen_fortran_interface_h > ${CMAKE_SOURCE_DIR}/include/fortran_interface.h
        )

add_custom_command (
        OUTPUT
        ${CMAKE_SOURCE_DIR}/gp/compinfo.h
        DEPENDS
        ${CMAKE_SOURCE_DIR}/utils/generate_compinfo.sh
        COMMAND
        ${CMAKE_SOURCE_DIR}/utils/generate_compinfo.sh > ${CMAKE_SOURCE_DIR}/gp/compinfo.h
        )

add_executable(gen_fortran_interface_h ${CMAKE_SOURCE_DIR}/build/gen_fortran_interface_h.c)

include(cmake/Sources.cmake)
include(cmake/Definitions.cmake)

set (GENERATED_FILES
        ${CMAKE_SOURCE_DIR}/include/fortran_interface.h
        ${CMAKE_SOURCE_DIR}/gp/compinfo.h
        )

add_custom_target(interfaces ALL
        DEPENDS
        ${GENERATED_FILES}
)

add_executable(dirac.x
        ${DIRAC_SOURCES}
        ${GENERATED_FILES}
        )

set_target_properties(dirac.x
        PROPERTIES LINKER_LANGUAGE Fortran
        )

target_link_libraries(dirac.x blas lapack)
