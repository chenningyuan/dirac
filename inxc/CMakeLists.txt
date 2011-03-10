cmake_minimum_required(VERSION 2.6)

project(INXC Fortran C)

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

option(INXC_ENABLE_PROFILING      "Enable profiling"          OFF)
option(INXC_ENABLE_BOUNDS_CHECK   "Enable bounds check"       OFF)
option(INXC_ENABLE_CODE_COVERAGE  "Enable code coverage"      OFF)
option(INXC_ENABLE_EXTERNAL_MATH  "Use external math library" ON)
option(INXC_ENABLE_MPI            "Use MPI parallelization"   OFF)
option(INXC_ENABLE_XCFUN          "Use XCFun library"         ON)

#add_definitions(-DPRG_DIRAC)

include(cmake/Archs.cmake)
include(cmake/Compilers.cmake)

set(LIBS)

if(INXC_ENABLE_EXTERNAL_MATH)
    find_package(Math)
else()
    set(MATH_FOUND FALSE)
endif()

if(MATH_FOUND)
    set(LIBS
        ${LIBS}
        ${MATH_LIBS}
        )
else()
endif()

if(INXC_ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        add_definitions(-DVAR_MPI2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/../include
    ${PROJECT_SOURCE_DIR}/../dft
    ${PROJECT_SOURCE_DIR}/../gp
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

#file(GLOB INXC_SOURCES
#    ${PROJECT_SOURCE_DIR}/src/*.F90
#    ${PROJECT_SOURCE_DIR}/src/*.F
#    ${PROJECT_SOURCE_DIR}/../dft/*.c
#    )

set(INXC_SOURCES
    ${PROJECT_SOURCE_DIR}/src/ao_eval.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_block.F90
    ${PROJECT_SOURCE_DIR}/src/ao_eval_point.F90
    ${PROJECT_SOURCE_DIR}/src/ao_xyz_power.F90
    ${PROJECT_SOURCE_DIR}/src/basis_info.F90
    ${PROJECT_SOURCE_DIR}/src/basis_info_ng.F90
    ${PROJECT_SOURCE_DIR}/src/block_config.F90
    ${PROJECT_SOURCE_DIR}/src/density_eval.F90
    ${PROJECT_SOURCE_DIR}/src/dft_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/dft_functional_derv.F90
    ${PROJECT_SOURCE_DIR}/src/dft_input.F90
    ${PROJECT_SOURCE_DIR}/src/dftac.F
    ${PROJECT_SOURCE_DIR}/src/dirac_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/extra.F
    ${PROJECT_SOURCE_DIR}/src/main.F90
    ${PROJECT_SOURCE_DIR}/src/occupation_info.F90
    ${PROJECT_SOURCE_DIR}/src/xc_driver.F90
    ${PROJECT_SOURCE_DIR}/src/xc_mpi.F90
    ${PROJECT_SOURCE_DIR}/../dft/emb-functional.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-becke.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-cam.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-gga.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-gllb_hole.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-hjsx.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-kin_pw91k.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-kin_tf.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-kin_vw.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-kt.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-lb94.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-lbalpha.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-lyp.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-optx.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-p86c.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pbec.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pbex.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pw86x.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pw91c.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pw91x.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-pz81.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-revpbex.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-rpbex.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-slater.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-srpbec.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-srvwn5.c
    ${PROJECT_SOURCE_DIR}/../dft/fun-vwn.c
    ${PROJECT_SOURCE_DIR}/../dft/functionals.c
    ${PROJECT_SOURCE_DIR}/../dft/general.c
    )

add_executable(inxc.x
    ${INXC_SOURCES}
    )

target_link_libraries(
    inxc.x
    ${LIBS}
    )
