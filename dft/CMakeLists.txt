cmake_minimum_required(VERSION 2.6)

project(DFT Fortran)

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

option(ENABLE_PROFILING      "Enable profiling"          OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"       OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"      OFF)
option(ENABLE_EXTERNAL_MATH  "Use external math library" ON)
option(ENABLE_MPI            "Use MPI parallelization"   OFF)

include(cmake/Archs.cmake)
include(cmake/Compilers.cmake)

set(LIBS)

add_definitions(-DMOD_XCFUN)
add_subdirectory(xcfun)
set(LIBS
    ${LIBS}
    xcfun
    )

if(ENABLE_EXTERNAL_MATH)
    set(BLAS_LANG "Fortran")
    find_package(BLAS)
    if(BLAS_FOUND)
        set(LIBS ${LIBS} ${BLAS_LIBRARIES})
    endif()
    set(LAPACK_LANG "Fortran")
    find_package(LAPACK)
    if(LAPACK_FOUND)
        set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    endif()
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        add_definitions(-DVAR_MPI2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

include_directories(
    ${PROJECT_BINARY_DIR}/modules
    ${PROJECT_BINARY_DIR}/xcfun/fortran
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

set(SOURCES
    ${PROJECT_SOURCE_DIR}/../grid/num_grid_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval.F90
    ${PROJECT_SOURCE_DIR}/src/sigma_prefactor_setting.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_block.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_sub.F90
    ${PROJECT_SOURCE_DIR}/src/ao_xyz_power.F90
    ${PROJECT_SOURCE_DIR}/src/density_eval.F90
    ${PROJECT_SOURCE_DIR}/src/dft_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/external/character_processing.F90
    ${PROJECT_SOURCE_DIR}/src/external/keyword.F90
    ${PROJECT_SOURCE_DIR}/src/external/quaternion_algebra.F90
    ${PROJECT_SOURCE_DIR}/src/external/dirac_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/extra.F90
    ${PROJECT_SOURCE_DIR}/src/interface_ao.F90
    ${PROJECT_SOURCE_DIR}/src/interface_functional_read.F90
    ${PROJECT_SOURCE_DIR}/src/interface_grid.F90
    ${PROJECT_SOURCE_DIR}/src/interface_input.F90
    ${PROJECT_SOURCE_DIR}/src/interface_mo.F90
    ${PROJECT_SOURCE_DIR}/src/interface_file_io.F90
    ${PROJECT_SOURCE_DIR}/src/main.F90
    ${PROJECT_SOURCE_DIR}/src/xc_ac.F90
    ${PROJECT_SOURCE_DIR}/src/xc_blas3.F90
    ${PROJECT_SOURCE_DIR}/src/xc_derv.F90
    ${PROJECT_SOURCE_DIR}/src/xc_driver.F90
    ${PROJECT_SOURCE_DIR}/src/xc_max_block_length.F90
    ${PROJECT_SOURCE_DIR}/src/xc_mpi.F90
    ${PROJECT_SOURCE_DIR}/xcfun/fortran/xcfun_module.F90
    ${PROJECT_SOURCE_DIR}/xcfun/fortran/xcfun_autogen.F90
    )

add_executable(
    dft.x
    ${SOURCES}
    )

target_link_libraries(
    dft.x
    ${LIBS}
    )

include(Tests)
include(CTest)
enable_testing()
