cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(DIRAC Fortran C)
add_definitions(-DPRG_DIRAC)

# set version
set(VERSION_MAJOR 11)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

add_definitions(-DTHIS_IS_CMAKE_BUILD)

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

option(ENABLE_INTERNAL_MATH  "Enable internal math library"       OFF)
option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"             OFF)
option(ENABLE_MPI            "Enable MPI parallelization"         OFF)
option(ENABLE_MPI2_DETECTION "Enable detection of MPI-2 standard" ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"      OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"       OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int module"              OFF)
option(ENABLE_XCFUN          "Enable XCFun library"               ON)
option(ENABLE_PROFILING      "Enable profiling"                   OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"               OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking"    OFF)
option(ENABLE_UTILITIES      "Enable utitilites"                  OFF)

# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest's own timeout
#miro: this is not properly working as no stdout/err is displayed on CDashboard!
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

option(FORCE_DIRAC_BLAS   "Force use of DIRAC BLAS"   OFF)
option(FORCE_DIRAC_LAPACK "Force use of DIRAC LAPACK" OFF)
mark_as_advanced(FORCE_DIRAC_BLAS)
mark_as_advanced(FORCE_DIRAC_LAPACK)

include(Math)
if(NOT BLAS_FOUND OR FORCE_DIRAC_BLAS)
    message("-- Using own BLAS library (slow)")
    add_definitions(-DUSE_DIRAC_BLAS)
endif()
if(NOT LAPACK_FOUND OR FORCE_DIRAC_LAPACK)
    message("-- Using own LAPACK library (slow)")
    add_definitions(-DUSE_DIRAC_DLAMCH)
    add_definitions(-DUSE_DIRAC_DLAPACK)
    add_definitions(-DUSE_DIRAC_ZLAPACK)
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

        if(ENABLE_MPI2_DETECTION)
            # test whether MPI-2 is available
            include(MPI-2-test)
            if(MPI_2_SUPPORT_OK)
                add_definitions(-DVAR_MPI2)
                message("-- MPI-2 support found, define -DVAR_MPI2")
            else()
                message("-- no MPI-2 support found, will try with MPI-1")
            endif()
        else()
            message("-- MPI-2 support disabled, will try with MPI-1")
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
    ${PROJECT_SOURCE_DIR}/utils/gen_fortran_interface_h.c
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

include(CompilationInfo)

# copy pam to build dir
configure_file(
    ${CMAKE_SOURCE_DIR}/pam
    ${CMAKE_BINARY_DIR}/pam
    COPYONLY
    )

set(GENERATED_FILES
    ${PROJECT_SOURCE_DIR}/include/fortran_interface.h
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
    gp/compilation_info.F90
    gp/git_revision_info.F90
    )

add_library(
    dirac 
    ${DIRAC_SOURCES} 
    ${GENERATED_FILES} 
    )


set(LIST_OF_EXECUTABLES)
macro(compile_standalone _executable _source LIST_OF_EXECUTABLES)
    add_executable(
        ${_executable}
        ${_source}
    )
    set(LIST_OF_EXECUTABLES
        ${LIST_OF_EXECUTABLES}
        ${_executable}
    )
endmacro()
compile_standalone(dirac.x                   main/main.F90                     "${LIST_OF_EXECUTABLES}")
if(ENABLE_UTILITIES)
    compile_standalone(dirac_mointegral_export.x utils/dirac_mointegral_export.F90 "${LIST_OF_EXECUTABLES}")
    compile_standalone(rspread.x                 utils/rspread.F                   "${LIST_OF_EXECUTABLES}")
    compile_standalone(labread.x                 utils/labread.F                   "${LIST_OF_EXECUTABLES}")
    compile_standalone(twofit.x                  utils/twofit.F                    "${LIST_OF_EXECUTABLES}")
    compile_standalone(cfread.x                  utils/cfread.F                    "${LIST_OF_EXECUTABLES}")
    compile_standalone(vibcal.x                  utils/vibcal.F                    "${LIST_OF_EXECUTABLES}")
    compile_standalone(polfit.x                  utils/polfit.F                    "${LIST_OF_EXECUTABLES}")
    compile_standalone(mx2fit.x                  utils/mx2fit.F                    "${LIST_OF_EXECUTABLES}")
endif()
foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    target_link_libraries(
        ${EXECUTABLE}
        dirac
        ${LIBS}
        )
endforeach()

include(Install)
include(Release)

if(VERBOSE_OUTPUT)
    include(Report)
endif()

include(Documentation)

# add short make target for the executable's diagnostics 
# with the  external python script
add_custom_target(dgn
    COMMAND ${CMAKE_SOURCE_DIR}/utils/diagn.py
)

# custom commands on created dirac.x executable
add_custom_command(
    TARGET Nightly Experimental dirac.x 
    if(VERBOSE_OUTPUT)
        POST_BUILD COMMAND ../utils/diagn.py
    endif()
)
