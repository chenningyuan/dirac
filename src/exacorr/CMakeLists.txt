include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )


set(FIXED_EXACORR_FORTRAN_SOURCES
    exacorr_dirac.F
)

set(FREE_EXACORR_FORTRAN_SOURCES
   exacorr_datatypes.F90
   exacorr_ao_to_mo.F90
   exacorr_cc.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_EXACORR_FORTRAN_SOURCES}  EXACORR COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_EXACORR_FORTRAN_SOURCES} EXACORR COMPILE_FLAGS "-qfixed")
endif()

add_library(
    exacorr
    OBJECT
    ${FIXED_EXACORR_FORTRAN_SOURCES}
    ${FREE_EXACORR_FORTRAN_SOURCES}
)

# deal with dependencies
add_dependencies(exacorr  main)
add_dependencies(exacorr  exatensor)
