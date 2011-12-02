cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#miro: to suppress Windows warning about policy...
cmake_policy(SET CMP0010 OLD)

project(DIRAC Fortran C)
add_definitions(-DPRG_DIRAC)


# set version
set(VERSION_MAJOR 11)
set(VERSION_MINOR 0)
set(VERSION_PATCH 2)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

include(CheckFortranSourceCompiles)

# for windows need substitute (type)
set(cat_program cat)

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

# check whether pamadm-generated files exist
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Sources.cmake"
    AND EXISTS "${CMAKE_SOURCE_DIR}/cmake/Definitions.cmake"
    AND EXISTS "${CMAKE_SOURCE_DIR}/cmake/Tests.cmake")
else()
    # they do not exist: run pamadm -a
    execute_process(
        COMMAND python pamadm -a
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_QUIET
    )
    message("-- You are lazy. Running pamadm -a for you ...")
endif()

option(ENABLE_INTERNAL_MATH   "Enable internal math library"                              OFF)
option(ENABLE_64BIT_INTEGERS  "Enable 64-bit integers"                                    OFF)
option(ENABLE_MPI             "Enable MPI parallelization"                                OFF)
option(ENABLE_MPI2_DETECTION  "Enable detection of MPI-2 standard"                        ON)
option(ENABLE_OMP             "Enable OpenMP parallelization"                             OFF)
option(ENABLE_CUDA            "Enable CUDA GPU acceleration"                              OFF)
option(ENABLE_PROFILING       "Enable profiling"                                          OFF)
option(ENABLE_BOUNDS_CHECK    "Enable bounds check"                                       OFF)
option(ENABLE_CODE_COVERAGE   "Enable code coverage"                                      OFF)
option(ENABLE_STATIC_LINKING  "Enable static libraries linking"                           OFF)
option(ENABLE_UTILITIES       "Enable utitilites"                                         OFF)
option(ENABLE_LARGE_MOLECULES "Increase common block parameters for large molecules"      OFF)
option(ENABLE_MATLAB_LOG      "Log interesting matrices in matlab format (very verbose!)" OFF)
option(ENABLE_SANITY_CHECKS   "Enable configure-time sanity checks"                       ON)

if(ENABLE_STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    if (ENABLE_MPI OR ENABLE_OMP)
        message(FATAL_ERROR "Static linking only works in sequential mode")
    endif()
endif()

option(FORCE_DIRAC_BLAS   "Force use of DIRAC BLAS"   OFF)
option(FORCE_DIRAC_LAPACK "Force use of DIRAC LAPACK" OFF)
mark_as_advanced(FORCE_DIRAC_BLAS)
mark_as_advanced(FORCE_DIRAC_LAPACK)

# set work memory
set(WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "Work memory in words"
    )
add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

if(ENABLE_LARGE_MOLECULES)
    add_definitions(-DLARGE_MOLECULES)
endif()

if(ENABLE_MATLAB_LOG)
    add_definitions(-DMOD_MATLAB_LOG)
endif()

include(Testing)
include(Definitions)
include(Archs)
include(SaveCompilerFlags)
include(FCompilers)
include(CCompilers)
include(Parallel)

set(LIBS)

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
    ${PROJECT_BINARY_DIR}/xcfun/fortran
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

add_subdirectory(xcfun)
set(LIBS
    ${LIBS}
    xcfun
    )

include(CompilationInfo)

# copy pam to build dir as pam and pam-dirac
configure_file(
    ${CMAKE_SOURCE_DIR}/pam
    ${CMAKE_BINARY_DIR}/pam
    COPYONLY
    )
configure_file(
    ${CMAKE_SOURCE_DIR}/pam
    ${CMAKE_BINARY_DIR}/pam-dirac
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
    )

add_library(
    dirac 
    ${DIRAC_SOURCES} 
    ${GENERATED_FILES} 
    )

# get size of static allocations in dirac.x
add_custom_target(
    info
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/static_size.py dirac.x
    )

# get size of static allocations in dirac.x and print sorted list of object 
add_custom_target(
    info_ext
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/static_size_ext.py dirac.x
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
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-E")
    endif()
    target_link_libraries(
        ${EXECUTABLE}
        dirac
        ${LIBS}
        )
endforeach()

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

include(Install)
include(Release)
include(Report)
include(Documentation)
