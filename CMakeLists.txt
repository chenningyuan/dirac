cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

set(CMAKE_BUILD_TYPE Debug)

project(Dirac Fortran C)

option(DIRAC_ENABLE_PROFILING     "Enable profiling"          OFF)
option(DIRAC_ENABLE_BOUNDS_CHECK  "Enable bounds check"       OFF)
option(DIRAC_ENABLE_CODE_COVERAGE "Enable code coverage"      OFF)
option(DIRAC_ENABLE_EXTERNAL_MATH "Use external math library" ON)
option(DIRAC_ENABLE_MPI           "Use MPI parallelization"   OFF)
option(DIRAC_ENABLE_XCFUN         "Use XCFun library"         ON)
option(DIRAC_ENABLE_UTILITIES     "Compile utilities"         OFF)
option(DIRAC_ENABLE_MAIN          "Compile main binary"       ON)

add_definitions(-DPRG_DIRAC)

set(DIRAC_WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "DIRAC work memory in words"
    )
add_definitions(-DINSTALL_WRKMEM=${DIRAC_WORK_MEM_WORDS})

# radovan: note to myself for later
# mark_as_advanced(VALGRIND_EXECUTABLE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYS_LINUX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "i686")
        add_definitions(-DARCH32BIT)
    elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    else()
    endif()
else()
endif()

include(ConfigCompilerFlags)
if (CMAKE_COMPILER_IS_GNUC)
    set (CMAKE_C_FLAGS "-Wall")
    set (CMAKE_C_FLAGS_DEBUG "-O0 -g3")
    set (CMAKE_C_FLAGS_RELEASE "-g -O2 -Wno-unused")
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
    set (CMAKE_C_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set (CMAKE_C_FLAGS_DEBUG "-g -O0")
    set (CMAKE_C_FLAGS_RELEASE "-g -O2")
    set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif ()
config_compiler_flags(C)

if (CMAKE_COMPILER_IS_GNUFORTRAN)
#   set(CMAKE_Fortran_FLAGS         "-Wall")
    set(CMAKE_Fortran_FLAGS         "     ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    add_definitions(-DVAR_GFORTRAN)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    add_definitions(-DVAR_IFC)
endif ()
config_compiler_flags(Fortran)

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

if(DIRAC_ENABLE_MAIN)
    add_executable(dirac.x
        ${DIRAC_SOURCES}
        ${GENERATED_FILES}
        )

    target_link_libraries(
        dirac.x
        ${DIRAC_LIBRARIES_TO_LINK}
        )
endif()

if(DIRAC_ENABLE_UTILITIES)
    add_executable(dirac_mointegral_export.x
        ${PROJECT_SOURCE_DIR}/utils/dirac_mointegral_export.F90
        )
    add_executable(rspread.x
        ${PROJECT_SOURCE_DIR}/utils/rspread.F
        )
    add_executable(labread.x
        ${PROJECT_SOURCE_DIR}/utils/labread.F
        )
    add_library(utils_dependencies
        ${PROJECT_SOURCE_DIR}/abacus/hergp.F
        ${PROJECT_SOURCE_DIR}/dirac/codata.F90
        ${PROJECT_SOURCE_DIR}/gp/execinfo.c
        ${PROJECT_SOURCE_DIR}/gp/gpc.c
        ${PROJECT_SOURCE_DIR}/gp/gphjj.F
        ${PROJECT_SOURCE_DIR}/gp/gpjeppe.F
        ${PROJECT_SOURCE_DIR}/gp/gpsaue.F
        ${PROJECT_SOURCE_DIR}/gp/gptrygve.F
        ${PROJECT_SOURCE_DIR}/gp/memory_allocator.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_control.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_errorhandler.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_internal_ada_1d.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_internal_ada_2d.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_internal_ada_3d.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_internal_ada_4d.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_internal_init.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_parameters.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_track_if.F90
        ${PROJECT_SOURCE_DIR}/gp/memory_tracker.F90
        ${PROJECT_SOURCE_DIR}/gp/mempkg.F
        ${PROJECT_SOURCE_DIR}/gp/qpack.F
        ${PROJECT_SOURCE_DIR}/gp/stacktraverse.c
        ${PROJECT_SOURCE_DIR}/gp/tracing.c
        ${PROJECT_SOURCE_DIR}/gp/tracker_util.c
        ${PROJECT_SOURCE_DIR}/pdpack/arhpack.F
        ${PROJECT_SOURCE_DIR}/pdpack/dge.F
        ${PROJECT_SOURCE_DIR}/pdpack/dsp.F
        ${PROJECT_SOURCE_DIR}/pdpack/eispack.F
        ${PROJECT_SOURCE_DIR}/pdpack/linextra.F
        ${PROJECT_SOURCE_DIR}/pdpack/printpkg.F
        )
    add_executable(twofit.x
        ${PROJECT_SOURCE_DIR}/utils/twofit.F
        )
    target_link_libraries(
        twofit.x
        utils_dependencies
        ${DIRAC_LIBRARIES_TO_LINK}
        )
    add_executable(cfread.x
        ${PROJECT_SOURCE_DIR}/utils/twofit.F
        )
    target_link_libraries(
        cfread.x
        utils_dependencies
        ${DIRAC_LIBRARIES_TO_LINK}
        )
    add_executable(vibcal.x
        ${PROJECT_SOURCE_DIR}/utils/vibcal.F
        )
    target_link_libraries(
        vibcal.x
        utils_dependencies
        ${DIRAC_LIBRARIES_TO_LINK}
        )
    add_executable(polfit.x
        ${PROJECT_SOURCE_DIR}/utils/polfit.F
        )
    target_link_libraries(
        polfit.x
        utils_dependencies
        ${DIRAC_LIBRARIES_TO_LINK}
        )
    add_executable(mx2fit.x
        ${PROJECT_SOURCE_DIR}/utils/mx2fit.F
        )
    target_link_libraries(
        mx2fit.x
        utils_dependencies
        ${DIRAC_LIBRARIES_TO_LINK}
        )

endif()
