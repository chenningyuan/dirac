include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For librelccsd.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_RELCCSD_FORTRAN_SOURCES
    ccdriv.F
    cceqns.F
    ccfspc.F
    ccgenp.F
    ccgetv.F
    ccgrad.F
    ccintm.F
    ccints.F
    ccmain.F
    ccpam.F
    ccrpa.F
    ccsort.F
    ccxblas.F
    waio.F
    ${PROJECT_SOURCE_DIR}/src/prp/mp2lag.F
    ${PROJECT_SOURCE_DIR}/src/gp/gpsaue.F
)

set(FREE_RELCCSD_FORTRAN_SOURCES
    relcc_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/gp/machine_parameters.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_allocator.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_control.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_errorhandler.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_internal_ada_1d.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_internal_ada_2d.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_internal_ada_3d.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_internal_ada_4d.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_internal_init.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_legacy_dirac.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_parameters.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_track_if.F90
    ${PROJECT_SOURCE_DIR}/src/gp/memory_tracker.F90
    ${PROJECT_SOURCE_DIR}/src/gp/xmlout.F90
    ${PROJECT_SOURCE_DIR}/src/gp/exception_handling.F90
    ${PROJECT_SOURCE_DIR}/src/gp/quaternion_algebra.F90
    ${PROJECT_SOURCE_DIR}/src/dirac/codata.F90
    ${PROJECT_SOURCE_DIR}/src/interface_mpi/integer_kind_mpilib.F90
    ${PROJECT_SOURCE_DIR}/src/interface_mpi/integer_model.F90
    ${PROJECT_SOURCE_DIR}/src/interface_mpi/interface_to_mpi.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

#miro: all these manually picked up sources are necessary for building up the relccsd library, "make relccsd"
# -  the standalone executable, relccsd.x, however, depends on the huge libdirac.a library
# - the quit.F file is not listed due to replacement by DIRAC QUIT routine
# - the ccunix.c is not included
add_library(
    relccsd
    ${FIXED_RELCCSD_FORTRAN_SOURCES}
    ${FREE_RELCCSD_FORTRAN_SOURCES}
)
