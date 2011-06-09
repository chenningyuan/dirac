cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(DIRAC Fortran C)
add_definitions(-DPRG_DIRAC)

# set version
set(VERSION_MAJOR 11)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(THIS_IS_CMAKE_BUILD 1)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

# set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
		AND NOT cmake_build_type_tolower STREQUAL "release"
		AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
	message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
	Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

# do not rebuild if rules (compiler flags) change
# this has advantages and disadvantages
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# check whether pamadm-generated files exist, if not exit
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Sources.cmake"
    AND EXISTS "${CMAKE_SOURCE_DIR}/cmake/Definitions.cmake"
    AND EXISTS "${CMAKE_SOURCE_DIR}/cmake/Tests.cmake")
else()
    message(FATAL_ERROR "could not find cmake/Sources.cmake and/or cmake/Definitions.cmake and/or cmake/Tests.cmake; run pamadm first")
endif()

include(Definitions)

option(ENABLE_EXTERNAL_MATH  "Enable external math library"  ON)
option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"        OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int module"         OFF)
option(ENABLE_XCFUN          "Enable XCFun library"          ON)
option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)

# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest own timeout

set(DART_TESTING_TIMEOUT
    "1200"
    CACHE STRING
    "Set default timeout in seconds for every single test"
    )

#
# set additional runtest flags - number of CPU and runtest timeout
#
# Miro: for CTest you can not have more than one $-variable in the add_test $-parameter !
#
if(NOT "${RUNTEST_NR_PROC}" STREQUAL "")
    set(RUNTEST_FLAG_NR_PROC "--mpi=${RUNTEST_NR_PROC}")
    message("-- runtest extra flag active, RUNTEST_FLAG_NR_PROC=\"${RUNTEST_FLAG_NR_PROC}\" ")
endif()
if(NOT "${RUNTEST_TIMEOUT}" STREQUAL "")
    set(RUNTEST_FLAG_TIMEOUT "--timeout=${RUNTEST_TIMEOUT}")
    message("-- runtest extra flag active, RUNTEST_FLAG_TIMEOUT=\"${RUNTEST_FLAG_TIMEOUT}\" ")
endif()

include(Tests)
include(CTest)
enable_testing()

# set work memory
set(WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "Work memory in words"
    )
add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

include(Archs)
include(Compilers)

set(LIBS)

include(Math)
if(NOT MATH_FOUND)
    message("-- Using own BLAS and LAPACK math")
    add_definitions(-DUSE_DIRAC_BLAS)
    add_definitions(-DUSE_DIRAC_DLAMCH)
    add_definitions(-DUSE_DIRAC_DLAPACK)
    add_definitions(-DUSE_DIRAC_ZLAPACK)
endif()

option(FORCE_DIRAC_BLAS    "Force use of DIRAC BLAS"    OFF)
option(FORCE_DIRAC_DLAMCH  "Force use of DIRAC DLAMCH"  OFF)
option(FORCE_DIRAC_DLAPACK "Force use of DIRAC DLAPACK" OFF)
option(FORCE_DIRAC_ZLAPACK "Force use of DIRAC ZLAPACK" OFF)
mark_as_advanced(FORCE_DIRAC_BLAS)
mark_as_advanced(FORCE_DIRAC_DLAMCH)
mark_as_advanced(FORCE_DIRAC_DLAPACK)
mark_as_advanced(FORCE_DIRAC_ZLAPACK)
if(FORCE_DIRAC_BLAS)
    add_definitions(-DUSE_DIRAC_BLAS)
    message(" -- Force use of DIRAC BLAS")
endif()
if(FORCE_DIRAC_DLAMCH)
    add_definitions(-DUSE_DIRAC_DLAMCH)
    message("-- Force use of DIRAC DLAMCH")
endif()
if(FORCE_DIRAC_DLAPACK)
    add_definitions(-DUSE_DIRAC_DLAPACK)
    message("-- Force use of DIRAC DLAPACK")
endif()
if(FORCE_DIRAC_ZLAPACK)
    add_definitions(-DUSE_DIRAC_ZLAPACK)
    message("-- Force use of DIRAC ZLAPACK")
endif()

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DINT_STAR8)
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})

        # test whether MPI-2 is available
        include(MPI-2-test)
        if(MPI_2_SUPPORT_OK)
            add_definitions(-DVAR_MPI2)
            message("-- MPI-2 support found, define -DVAR_MPI2")
        else()
            message("-- no MPI-2 support found, will try with MPI-1")
        endif()
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

if(ENABLE_CUDA)
    find_package(CUDA)
    if(CUDA_FOUND)
        set(LIBS
            ${LIBS}
            ${CUDA_CUBLAS_LIBRARIES}
            )
        add_definitions(-DUSE_CUBLAS)
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

set(GENERATED_FILES
    ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/modules
    ${PROJECT_BINARY_DIR}/xcint/xcfun/fortran
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

if(ENABLE_XCFUN)
    add_subdirectory(xcint/xcfun)
    set(LIBS
        ${LIBS}
        xcfun
        )
endif()

# get svn or git revision number
set(SVN_REVISION "")
set(GIT_REVISION "")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    execute_process(
        COMMAND ./get-git-revision
        TIMEOUT 1
        OUTPUT_VARIABLE GIT_REVISION
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND ./get-svn-revision
        TIMEOUT 1
        OUTPUT_VARIABLE SVN_REVISION
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

#fixme: it's cmake time, not compile time
#       also revision numbers are cmake time
# get compilation time
set(COMPILATION_TIME "")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    execute_process(
        COMMAND date
        TIMEOUT 1
        OUTPUT_VARIABLE COMPILATION_TIME
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

# configure compinfo.h
configure_file(
    ${CMAKE_SOURCE_DIR}/gp/compinfo.h.in
    ${CMAKE_SOURCE_DIR}/gp/compinfo.h
    )

add_custom_target(dirac_interfaces
    ALL
    DEPENDS
    ${GENERATED_FILES}
)

include(Sources)

set(DIRAC_SOURCES
    ${C_SOURCES}
    ${FREE_FORTRAN_SOURCES}
    ${FIXED_FORTRAN_SOURCES}
    )

add_library(
    dirac
    ${DIRAC_SOURCES}
    ${GENERATED_FILES}
    )

add_executable(
    dirac.x
    main/main.F90
    )

# custom commands on created dirac.x executable
add_custom_command(
TARGET dirac.x 
POST_BUILD COMMAND 
../utils/diagn.py
)

add_executable(
    dirac_mointegral_export.x
    ${PROJECT_SOURCE_DIR}/utils/dirac_mointegral_export.F90
    )
add_executable(
    rspread.x
    ${PROJECT_SOURCE_DIR}/utils/rspread.F
    )
add_executable(
    labread.x
    ${PROJECT_SOURCE_DIR}/utils/labread.F
    )
add_executable(
    twofit.x
    ${PROJECT_SOURCE_DIR}/utils/twofit.F
    )
add_executable(
    cfread.x
    ${PROJECT_SOURCE_DIR}/utils/twofit.F
    )
add_executable(
    vibcal.x
    ${PROJECT_SOURCE_DIR}/utils/vibcal.F
    )
add_executable(
    polfit.x
    ${PROJECT_SOURCE_DIR}/utils/polfit.F
    )
add_executable(
    mx2fit.x
    ${PROJECT_SOURCE_DIR}/utils/mx2fit.F
    )
foreach(
    BINARY
    dirac.x
    dirac_mointegral_export.x
    rspread.x
    labread.x
    twofit.x
    cfread.x
    vibcal.x
    polfit.x
    mx2fit.x
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${BINARY} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    target_link_libraries(
        ${BINARY}
        dirac
        ${LIBS}
        )
endforeach()

include(Install)
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

SET(EXPORT_DIR ${CMAKE_BINARY_DIR}/export)
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${EXPORT_DIR};/"
    )

include(CPack)

file(GLOB_RECURSE H_FILES   RELATIVE ${CMAKE_SOURCE_DIR} "*.h")
file(GLOB_RECURSE INC_FILES RELATIVE ${CMAKE_SOURCE_DIR} "*.inc")

set(INCLUDE_FILES
    ${H_FILES}
    ${INC_FILES}
    )

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#miro TODO: adapt under Windows
#radovan: don't agree - i don't see why we should ever
#         want to prepare the release on a windows machine
    add_custom_target(release
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND tar cf ${EXPORT_DIR}/sources.tar -C ${CMAKE_SOURCE_DIR} ${DIRAC_SOURCES} ${INCLUDE_FILES}
        COMMAND tar xf ${EXPORT_DIR}/sources.tar -C ${EXPORT_DIR}
        COMMAND rm     ${EXPORT_DIR}/sources.tar
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/cmake              ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/cbuild            ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis             ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis_dalton      ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/test              ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/utils             ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/modules           ${EXPORT_DIR}
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/xcint/xcfun       ${EXPORT_DIR}/xcint
        COMMAND cp    ${CMAKE_SOURCE_DIR}/gp/compinfo.h.in  ${EXPORT_DIR}/gp
        COMMAND cp    ${CMAKE_SOURCE_DIR}/pam               ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/runtest           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/setup             ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/setup.examples    ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/CMakeLists.txt    ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/COPYING           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/CTestConfig.cmake ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/INSTALL           ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.guess      ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.status     ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/config.sub        ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/configure         ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/Makefile          ${EXPORT_DIR}
        COMMAND cp    ${CMAKE_SOURCE_DIR}/README            ${EXPORT_DIR}
        COMMAND make package_source
        COMMENT "Packaging source files"
        VERBATIM
        )
endif()

# add short make target for the executable's diagnostics 
# with the  external python script
add_custom_target(dgn
    COMMAND ${CMAKE_SOURCE_DIR}/utils/diagn.py
)

include(Report)

# print out enviroment variables
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message("Actual Unix environment variable PATH=$ENV{PATH}")
    message("Actual Unix environment variable LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}")
else()
    message("Actual Windows environment variable PATH=$ENV{PATH}")
	message("Actual Windows environment variable LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}")
endif()

# handle the .diracrc file containing pam flags
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(DIRACRC_FILE "$ENV{HOME}/.diracrc")
else()
    set(DIRACRC_FILE "$ENV{UserProfile}/.diracrc")
endif()	
	
if(EXISTS "${DIRACRC_FILE}")
    message("Configuration file ${DIRACRC_FILE} exists.")
    # TODO: print content of this file
else()
    configure_file(
        ${CMAKE_SOURCE_DIR}/cmake/diracrc.in
        ${DIRACRC_FILE}
        )
    message("Configuration file ${DIRACRC_FILE} did not exist.")
    message("This template file was created by cmake. You should verify/edit it during/after the compilation.")
    message("You need to set at least --scratch in ${DIRACRC_FILE} for runtest and make test to work.")
endif()
# read and print out content of the .diracrc file
file (READ "${DIRACRC_FILE}" DIRACRC_CONTENT)
message("Content of the ${DIRACRC_FILE} configuration file: \n${DIRACRC_CONTENT}" )


