cmake_minimum_required(VERSION 2.6)

project(DIRAC Fortran C)

set(VERSION_MAJOR 10)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

# do not rebuild if rules (compiler flags) change
# this has advantages and disadvantages
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Sources.cmake")
    include(Sources)
else()
    message(FATAL_ERROR "could not find cmake/Sources.cmake; run pamadm first")
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Definitions.cmake")
    include(Definitions)
else()
    message(FATAL_ERROR "could not find cmake/Definitions.cmake; run pamadm first")
endif()

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"      OFF)
option(ENABLE_PROFILING      "Enable profiling"            OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"         OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"        OFF)
option(ENABLE_EXTERNAL_MATH  "Use external math library"   ON)
option(ENABLE_MPI            "Use MPI parallelization"     OFF)
option(ENABLE_OMP            "Use OpenMP parallelization"  OFF)
option(ENABLE_XCFUN          "Use XCFun library"         ON)
option(ENABLE_UTILITIES      "Compile utilities"         OFF)
option(ENABLE_TESTS          "Enable make test"          ON)
option(USE_LOCAL_DISKS             "Use local disks"           ON)

if(ENABLE_TESTS)
    SET(BUILDNAME "${BUILDNAME}" 
            CACHE STRING "Name of build on the dashboard"
        )
    message("-- DIRAC: BUILDNAME is \"${BUILDNAME}\"")
    if(NOT "${RUNTEST_TIMEOUT}" STREQUAL "")
        set(DART_TESTING_TIMEOUT
            "${RUNTEST_TIMEOUT}"
            CACHE STRING
            "Set timeout in seconds for every single test."
            FORCE
            )
        message("-- DIRAC: DART_TESTING_TIMEOUT (time limit for single test) set to ${DART_TESTING_TIMEOUT}")
    endif()
    if(NOT "${RUNTEST_NR_PROC}" STREQUAL "")
        message("-- DIRAC: RUNTEST_NR_PROC set to ${RUNTEST_NR_PROC}")
        set(RUNTEST_FLAGS "--mpi=${RUNTEST_NR_PROC}" )
        message("-- DIRAC: RUNTEST_FLAGS set to \"${RUNTEST_FLAGS}\"")
    endif()
    if(NOT "${RUNTEST_OWN_TIMEOUT_FLAG}" STREQUAL "")
        set(RUNTEST_TIMEOUT_FLAG "--timeout=${RUNTEST_OWN_TIMEOUT_FLAG}")
        message("-- DIRAC: RUNTEST_TIMEOUT_FLAG set to \"${RUNTEST_TIMEOUT_FLAG}\"")
    endif()
    include(CTest)
    enable_testing()
    if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Tests.cmake")
        include(Tests)
    else()
        message(FATAL_ERROR "could not find cmake/Tests.cmake; run pamadm first")
    endif()
endif()

add_definitions(-DPRG_DIRAC)
if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DINT_STAR8)
endif()

set(WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "DIRAC work memory in words"
    FORCE
    )

# find a good default scratch base directory
# first try scratch/, then work/, then tmp/
if(EXISTS "/scratch/")
    set(DEFAULT_SCRATCH_BASE_DIR "/scratch")
else()
    if(EXISTS "/work/")
        set(DEFAULT_SCRATCH_BASE_DIR "/work")
    else()
        set(DEFAULT_SCRATCH_BASE_DIR "/tmp")
    endif()
endif()

if(NOT SCRATCH_DIR)
    # use default
    set(SCRATCH_DIR "${DEFAULT_SCRATCH_BASE_DIR}/scratch_$ENV{USER}")
else()
    # use the one chosen by user
    set(SCRATCH_DIR "${SCRATCH_DIR}")
endif()

# set scratch directory
set(SCRATCH_DIR
    "${SCRATCH_DIR}"
    CACHE STRING
    "DIRAC scratch directory"
    FORCE
    )

message("-- DIRAC: SCRATCH_DIR is \"${SCRATCH_DIR}\"")
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile.config "SCRATCH_DIR=${SCRATCH_DIR}")

find_program(PROFILER gprof)
if(PROFILER)
    file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nPROFILER=${PROFILER}")
endif()

find_program(VALGRIND_EXECUTABLE valgrind)
if(VALGRIND_EXECUTABLE)
    file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nVALGRIND=${VALGRIND_EXECUTABLE}")
endif()

if(USE_LOCAL_DISKS)
    add_definitions(-DSYS_LOCALDISKS)
    file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nLOCALDISKS=yes")
else()
    file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nLOCALDISKS=no")
endif()

find_program(PBS_ENVIRONMENT pbsnodes)
if(PBS_ENVIRONMENT)
    # this is a PBS batch system
    file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nMACH_FILE=$PBS_NODEFILE")
endif()

add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

include(Archs)
include(Compilers)

set(LIBS)

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

if(ENABLE_OMP)
    find_package(OpenMP)
    if(OpenMP_FOUND)
        add_definitions(-DVAR_OMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_COMPILE_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

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
    ${PROJECT_SOURCE_DIR}/cbuild/gen_fortran_interface_h.c
    )

add_custom_command(
    OUTPUT
    ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    DEPENDS
    ${PROJECT_SOURCE_DIR}/utils/generate_compinfo.py
    COMMAND
    python ${PROJECT_SOURCE_DIR}/utils/generate_compinfo.py --fc=${CMAKE_Fortran_COMPILER} --cc=${CMAKE_C_COMPILER} --libs="${LIBS}" > ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    )

set(GENERATED_FILES
    ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
    ${PROJECT_SOURCE_DIR}/gp/compinfo.h
    )

set(INTERFACES_LIST
    ${GENERATED_FILES}
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/modules
    ${PROJECT_BINARY_DIR}/xcfun/fortran
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

if(ENABLE_XCFUN)
    add_subdirectory(xcfun)
    set(INTERFACES_LIST
        ${INTERFACES_LIST}
        xcfun
        )
    set(LIBS
        ${LIBS}
        xcfun
        )
endif()

add_custom_target(dirac_interfaces
    ALL
    DEPENDS
    ${INTERFACES_LIST}
)

add_executable(dirac.x
    ${C_SOURCES}
    ${FREE_FORTRAN_SOURCES}
    ${FIXED_FORTRAN_SOURCES}
    ${GENERATED_FILES}
    )

# radovan: Miro, you can add link flags like this: 
# SET_TARGET_PROPERTIES(dirac.x PROPERTIES LINK_FLAGS "-Wl,-E")
# radovan: i have commented it - does not work on AIX, do we really need it?

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    SET_TARGET_PROPERTIES(dirac.x PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
endif()

target_link_libraries(
    dirac.x
    ${LIBS}
    )

include(Install)

if(ENABLE_UTILITIES)
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
        ${LIBS}
        )
    add_executable(cfread.x
        ${PROJECT_SOURCE_DIR}/utils/twofit.F
        )
    target_link_libraries(
        cfread.x
        utils_dependencies
        ${LIBS}
        )
    add_executable(vibcal.x
        ${PROJECT_SOURCE_DIR}/utils/vibcal.F
        )
    target_link_libraries(
        vibcal.x
        utils_dependencies
        ${LIBS}
        )
    add_executable(polfit.x
        ${PROJECT_SOURCE_DIR}/utils/polfit.F
        )
    target_link_libraries(
        polfit.x
        utils_dependencies
        ${LIBS}
        )
    add_executable(mx2fit.x
        ${PROJECT_SOURCE_DIR}/utils/mx2fit.F
        )
    target_link_libraries(
        mx2fit.x
        utils_dependencies
        ${LIBS}
        )
endif()

include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "T. Saue, L. Visscher, H. J. Aa. Jensen, and contributors")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "R. Bast")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL")
set(CPACK_RESOURCE_FILE_LICENSE    "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR    "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR    "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH    "${VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DIRAC ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
    )

include(CPack)
