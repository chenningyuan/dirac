cmake_minimum_required(VERSION 2.6)

project(DIRAC Fortran C)

set(DIRAC_VERSION_MAJOR 10)
set(DIRAC_VERSION_MINOR 1)
set(DIRAC_VERSION_PATCH 0)

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

option(DIRAC_ENABLE_64BIT_INTEGERS "Enable 64-bit integers"    OFF)
option(DIRAC_ENABLE_PROFILING      "Enable profiling"          OFF)
option(DIRAC_ENABLE_BOUNDS_CHECK   "Enable bounds check"       OFF)
option(DIRAC_ENABLE_CODE_COVERAGE  "Enable code coverage"      OFF)
option(DIRAC_ENABLE_EXTERNAL_MATH  "Use external math library" ON)
option(DIRAC_ENABLE_MPI            "Use MPI parallelization"   OFF)
option(DIRAC_ENABLE_RSH            "Use rsh/rcp for distribution over nodes "   OFF)
option(DIRAC_ENABLE_XCFUN          "Use XCFun library"         ON)
option(DIRAC_ENABLE_UTILITIES      "Compile utilities"         OFF)
option(DIRAC_ENABLE_TESTS          "Enable make test"          ON)
option(USE_LOCAL_DISKS             "Use local disks"           ON)

if(DIRAC_ENABLE_TESTS)
    SET(BUILDNAME "${BUILDNAME}" 
            CACHE STRING "Name of build on the dashboard"
        )
    message("-- DIRAC: BUILDNAME is \"${BUILDNAME}\"")
    if(NOT "${DIRAC_RUNTEST_TIMEOUT}" STREQUAL "")
        set(DART_TESTING_TIMEOUT
            "${DIRAC_RUNTEST_TIMEOUT}"
            CACHE STRING
            "Set timeout in seconds for every single test."
            FORCE
            )
        message("-- DIRAC: DART_TESTING_TIMEOUT (time limit for single test) set to ${DART_TESTING_TIMEOUT}")
    endif()
    if(NOT "${DIRAC_RUNTEST_NR_PROC}" STREQUAL "")
        message("-- DIRAC: DIRAC_RUNTEST_NR_PROC set to ${DIRAC_RUNTEST_NR_PROC}")
        set(DIRAC_RUNTEST_FLAGS "--mpi=${DIRAC_RUNTEST_NR_PROC}" )
        message("-- DIRAC: RUNTEST_FLAGS set to \"${DIRAC_RUNTEST_FLAGS}\"")
    endif()
    if(NOT "${DIRAC_RUNTEST_OWN_TIMEOUT_FLAG}" STREQUAL "")
        set(DIRAC_RUNTEST_TIMEOUT_FLAG "--timeout=${DIRAC_RUNTEST_OWN_TIMEOUT_FLAG}")
        message("-- DIRAC: RUNTEST_TIMEOUT_FLAG set to \"${DIRAC_RUNTEST_TIMEOUT_FLAG}\"")
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
if(DIRAC_ENABLE_64BIT_INTEGERS)
    add_definitions(-DINT_STAR8)
endif()

set(DIRAC_WORK_MEM_WORDS
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

if(NOT DIRAC_SCRATCH_DIR)
    # use default
    set(DIRAC_SCRATCH_DIR "${DEFAULT_SCRATCH_BASE_DIR}/DIRAC_scratch_$ENV{USER}")
else()
    # use the one chosen by user
    set(DIRAC_SCRATCH_DIR "${DIRAC_SCRATCH_DIR}")
endif()

# set scratch directory
set(DIRAC_SCRATCH_DIR
    "${DIRAC_SCRATCH_DIR}"
    CACHE STRING
    "DIRAC scratch directory"
    FORCE
    )

message("-- DIRAC: DIRAC_SCRATCH_DIR is \"${DIRAC_SCRATCH_DIR}\"")
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile.config "SCRATCH_DIR=${DIRAC_SCRATCH_DIR}")


if(NOT DIRAC_ENABLE_RSH)
#  assign/find ssh/scp
     find_program(SCP_EXECUTABLE scp)
     if(SCP_EXECUTABLE)
         file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nRCP=${SCP_EXECUTABLE}")
     else()
          message(FATAL_ERROR "could not find scp executable !")
     endif()
     find_program(SSH_EXECUTABLE ssh)
     if(SSH_EXECUTABLE)
         file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nRSH=${SSH_EXECUTABLE}")
     else()
         message(FATAL_ERROR "could not find ssh executable !")
     endif()
     message("-- DIRAC: Nodes distribution tools are \"${SCP_EXECUTABLE}\" & \"${SSH_EXECUTABLE}\" ")
else()
#  assign/find rsh/rcp
     find_program(RCPEXECUTABLE rcp)
     if(RCPEXECUTABLE)
         file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nRCP=${RCPEXECUTABLE}")
     else()
         message(FATAL_ERROR "could not find rcp executable !")
     endif()
     find_program(RSHEXECUTABLE rsh)
     if(RSHEXECUTABLE)
         file(APPEND ${PROJECT_SOURCE_DIR}/Makefile.config "\nRSH=${RSHEXECUTABLE}")
     else()
         message(FATAL_ERROR "could not find rsh executable !")
     endif()
     message("-- DIRAC: Nodes distribution tools are \"${RCPEXECUTABLE}\" & \"${RSHEXECUTABLE}\" ")
endif()

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

add_definitions(-DINSTALL_WRKMEM=${DIRAC_WORK_MEM_WORDS})

include(Archs)
include(Compilers)

set(LIBS)

include(Math)

if(DIRAC_ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        add_definitions(-DVAR_MPI2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
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
    ${PROJECT_SOURCE_DIR}/build/gen_fortran_interface_h.c
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

set(DIRAC_INTERFACES_LIST
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

if(DIRAC_ENABLE_XCFUN)
    add_subdirectory(xcfun)
    set(DIRAC_INTERFACES_LIST
        ${DIRAC_INTERFACES_LIST}
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
    ${DIRAC_INTERFACES_LIST}
)

add_executable(dirac.x
    ${DIRAC_C_SOURCES}
    ${DIRAC_FREE_FORTRAN_SOURCES}
    ${DIRAC_FIXED_FORTRAN_SOURCES}
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
set(CPACK_PACKAGE_VERSION_MAJOR    "${DIRAC_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR    "${DIRAC_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH    "${DIRAC_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DIRAC ${DIRAC_VERSION_MAJOR}.${DIRAC_VERSION_MINOR}.${DIRAC_VERSION_PATCH}"
    )

include(CPack)
