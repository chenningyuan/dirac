cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

set(CMAKE_BUILD_TYPE Debug)

project(Dirac Fortran C)

option(DIRAC_ENABLE_PROFILING     "Enable profiling"          OFF)
option(DIRAC_ENABLE_BOUNDS_CHECK  "Enable bounds check"       OFF)
option(DIRAC_ENABLE_CODE_COVERAGE "Enable code coverage"      OFF)
option(DIRAC_ENABLE_EXTERNAL_MATH "Use external math library" ON)
option(DIRAC_ENABLE_MPI           "Use MPI parallelization"   OFF)
option(DIRAC_ENABLE_XCFUN         "Use XCFun library"         ON)

add_definitions(-DPRG_DIRAC)

set(DIRAC_WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "DIRAC work memory in words"
    )
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
    if (DIRAC_ENABLE_CODE_COVERAGE)
        set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} --coverage")
    endif ()
elseif (CMAKE_CC_COMPILER_ID MATCHES Intel)
    set(CMAKE_CC_FLAGS         "")
    set(CMAKE_CC_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG")
endif ()

set(DIRAC_LIBRARIES_TO_LINK)

if (DIRAC_ENABLE_EXTERNAL_MATH)
    find_package(BLAS)
    find_package(LAPACK)
    set(DIRAC_LIBRARIES_TO_LINK
        ${DIRAC_LIBRARIES_TO_LINK}
        blas
        lapack
        )
endif ()

if (DIRAC_ENABLE_MPI)
    find_package(MPI)
    add_definitions(-DVAR_MPI)
    if (MPI_FOUND)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif ()
endif ()

add_custom_command(
    OUTPUT
    ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
    DEPENDS
    gen_fortran_interface_h
    COMMAND
    gen_fortran_interface_h > ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
    )

add_executable(
    gen_fortran_interface_h
    ${PROJECT_SOURCE_DIR}/build/gen_fortran_interface_h.c
    )

add_custom_command(
    OUTPUT
    ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    DEPENDS
    ${PROJECT_SOURCE_DIR}/utils/generate_compinfo.sh
    COMMAND
    ${PROJECT_SOURCE_DIR}/utils/generate_compinfo.sh > ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    )

set(GENERATED_FILES
    ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
    ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    )

set(DIRAC_INTERFACES_LIST
    ${GENERATED_FILES}
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

include(cmake/Definitions.cmake)
include(cmake/Sources.cmake)

if (DIRAC_ENABLE_XCFUN)
    add_subdirectory(xcfun)
    set(DIRAC_INTERFACES_LIST
        ${DIRAC_INTERFACES_LIST}
        xcfun
        )
    set(DIRAC_LIBRARIES_TO_LINK
        ${DIRAC_LIBRARIES_TO_LINK}
        xcfun
        )
endif ()

add_custom_target(dirac_interfaces
    ALL
    DEPENDS
    ${DIRAC_INTERFACES_LIST}
)

add_executable(dirac.x
    ${DIRAC_SOURCES}
    ${GENERATED_FILES}
    )

target_link_libraries(
    dirac.x
    ${DIRAC_LIBRARIES_TO_LINK}
    )
