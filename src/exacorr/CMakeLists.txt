include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )


set(FIXED_EXACORR_FORTRAN_SOURCES
    qcorr_dirac.F
)

set(FREE_EXACORR_FORTRAN_SOURCES
   dirac_exatensor.F90
   qcorr_aomo_tranform.F90
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
