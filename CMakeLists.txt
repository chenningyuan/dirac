cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# suppress Windows warning about policy
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
    ${CMAKE_SOURCE_DIR}/cmake/pamadm-generated
    ${CMAKE_SOURCE_DIR}/cmake/architecture
    ${CMAKE_SOURCE_DIR}/cmake/parallel-environment
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/packaging
    ${CMAKE_SOURCE_DIR}/cmake/testing
    )

include(ConfigSafeGuards)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# check whether pamadm-generated files exist
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/SourcesList.cmake"      AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/DefinitionsList.cmake"  AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/IncludeFilesList.cmake" AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/TestsList.cmake")
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
option(ENABLE_PCMMOD          "Enable Polarizable Continuum Model through PCMSolver"      ON)

if(ENABLE_STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
endif()

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

include(ConfigTesting)
include(DefinitionsList)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigParallelEnvironment)

set(LIBS)

include(Math)

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
    xcfun
    ${LIBS}
    )

if(ENABLE_PCMMOD)
    include(ExternalProject)
    add_custom_target(external)
                                              
    set (ExternalProjectCMakeArgs  
    -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCBLAS_ROOT=${CBLAS_ROOT}
    -DEIGEN3_ROOT=${EIGEN3_ROOT}
    )
                                                            
    ExternalProject_Add(pcmsolver
    PREFIX ${PROJECT_SOURCE_DIR}/external
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pcmsolver
    BINARY_DIR ${PROJECT_SOURCE_DIR}/external/pcmsolver-build
    INSTALL_DIR ${PROJECT_SOURCE_DIR}/external
    CMAKE_ARGS ${ExternalProjectCMakeArgs}
    )                                                          
    add_dependencies(external pcmsolver)                       
endif()

include(BinaryInfo)

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

include(SourcesList)

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

# get size of static allocations
add_custom_target(
    info
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_static_size.py dirac.x
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

add_dependencies(dirac.x external)

if(ENABLE_PCMMOD)
    #  link_directories(${CMAKE_SOURCE_DIR}/external/lib)
    #  link_directories(${CMAKE_SOURCE_DIR}/external/pcmsolver/external/lib)
set(LIBS
    solver
	green
	cavity
	metal
	pedra
	wem
	wavcav
	getkw
    ${LIBS}
    )
# We should configure appropriately the pam script to
# read also PCMSolver input throught its separate
# input parsing script.
endif()

if(ENABLE_STATIC_LINKING AND CMAKE_Fortran_COMPILER_ID MATCHES GNU AND BLAS_FOUND AND LAPACK_FOUND)
    # Only newer versions of GNU-linkers have this flag ! Check "ld --help"
    # does not work with gfortran44 : http://mailman.ucar.edu/pipermail/wrf-users/2011/002172.html
    target_link_libraries(dirac.x -Wl,--no-export-dynamic)
    if (BLAS_TYPE MATCHES ATLAS)
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=30471#c7
# http://gcc.gnu.org/bugzilla/show_bug.cgi?format=multiple&id=42477
        target_link_libraries(dirac.x -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
    endif()
endif()

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

include(ConfigPackaging)

add_custom_target(
    html
    COMMAND sphinx-build -b html -d _build/doctrees ${CMAKE_SOURCE_DIR}/doc/sphinx _build/html
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
add_custom_target(
    latex
    COMMAND sphinx-build -b latex -d _build/doctrees ${CMAKE_SOURCE_DIR}/doc/sphinx _build/latex
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )

# run pamadm -a in source root
add_custom_target(
    adm
    COMMAND python pamadm -a
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
