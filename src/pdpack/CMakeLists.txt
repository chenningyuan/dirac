include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libpdpack.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_PDPACK_FORTRAN_SOURCES
    arhpack.F
    dge.F
    dlamch.F
    dsp.F
    eispack.F
    gp_auxlapack.F
    gp_dlapack.F
    gp_zlapack.F
    gp_zlinpack.F
    gpblas.F
    linextra.F
    printpkg.F

)

set(FREE_PDPACK_FORTRAN_SOURCES

)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_PDPACK_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_PDPACK_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    pdpack
    OBJECT
    ${FIXED_PDPACK_FORTRAN_SOURCES}
    ${FREE_PDPACK_FORTRAN_SOURCES}
)
