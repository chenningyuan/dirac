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
    ${PROJECT_SOURCE_DIR}/src/main/dirac.F
# needs fix missing object pamgeo
   ${PROJECT_SOURCE_DIR}/src/dirac/dirgeo.F
)

set(FREE_RELCCSD_FORTRAN_SOURCES
    relcc_cfg.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

#miro: all these manually picked up sources are necessary for building up the relccsd library, "make relccsd"
# - the quit.F file is not listed due to replacement by DIRAC QUIT routine
# - the ccunix.c is not included
add_library(
    relccsd
    ${FIXED_RELCCSD_FORTRAN_SOURCES}
    ${FREE_RELCCSD_FORTRAN_SOURCES}
)

# add intermodule dependencies
add_dependencies(relccsd gp)
add_dependencies(relccsd dirac)
add_dependencies(relccsd main)
