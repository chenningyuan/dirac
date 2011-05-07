cmake_minimum_required(VERSION 2.6)

project(DFT Fortran C)

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

option(DFT_ENABLE_PROFILING      "Enable profiling"          OFF)
option(DFT_ENABLE_BOUNDS_CHECK   "Enable bounds check"       OFF)
option(DFT_ENABLE_CODE_COVERAGE  "Enable code coverage"      OFF)
option(DFT_ENABLE_EXTERNAL_MATH  "Use external math library" ON)
option(DFT_ENABLE_MPI            "Use MPI parallelization"   OFF)
option(DFT_ENABLE_XCFUN          "Use XCFun library"         ON)

include(cmake/Archs.cmake)
include(cmake/Compilers.cmake)

set(LIBS)

if(DFT_ENABLE_EXTERNAL_MATH)
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

if(DFT_ENABLE_MPI)
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
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

set(DFT_SOURCES
    ${PROJECT_SOURCE_DIR}/src/quaternion_algebra.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_block.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_point.F90
    ${PROJECT_SOURCE_DIR}/src/ao_xyz_power.F90
    ${PROJECT_SOURCE_DIR}/src/interface_ao.F90
    ${PROJECT_SOURCE_DIR}/src/density_eval.F90
    ${PROJECT_SOURCE_DIR}/src/dft_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/dft_input.F90
    ${PROJECT_SOURCE_DIR}/src/xc_ac.F90
    ${PROJECT_SOURCE_DIR}/src/dirac_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/extra.F90
    ${PROJECT_SOURCE_DIR}/src/main.F90
    ${PROJECT_SOURCE_DIR}/src/num_grid_io.F90
    ${PROJECT_SOURCE_DIR}/src/interface_mo.F90
    ${PROJECT_SOURCE_DIR}/src/xc_blas3.F90
    ${PROJECT_SOURCE_DIR}/src/xc_derv.F90
    ${PROJECT_SOURCE_DIR}/src/xc_driver.F90
    ${PROJECT_SOURCE_DIR}/src/xc_max_block_length.F90
    ${PROJECT_SOURCE_DIR}/src/xc_mpi.F90
    ${PROJECT_SOURCE_DIR}/../grid/num_grid_cfg.F90
    )

add_executable(dft.x
    ${DFT_SOURCES}
    )

target_link_libraries(
    dft.x
    ${LIBS}
    )
