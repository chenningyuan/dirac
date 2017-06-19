include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libreladc.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_RELADC_FORTRAN_SOURCES
    adccore.F
    adcdble.F
    adcdiag.F
    adcsort.F
)

set(FREE_RELADC_FORTRAN_SOURCES
    adc_cfg.F90
    adc_fano_exchange.F90
    adc_fano_routines.F90
    adc_fano_real_routines.F90
    adc_fano_complex_routines.F90
    adc_fano_diag.F90
    adc_mat.F90
    adc_fano_matmul.F90
    adc_fano.F90
    adcutil.F90
    polprp_buildtm.F
    polprp_cfg.F90
    polprp_cousat.F
    polprp_davidson.F
    polprp_driver.F
    polprp_helpers.F
    polprp_mainblk.F
    polprp_mpimaster.F
    polprp_mpiserver.F
    polprp_prpint.F
    polprp_qstackm.F90
)

set(C_RELADC_SOURCES
    polprp_qstack.c
)

if (ENABLE_STIELTJES)
  #miro: needs extra sources !
   set(FREE_RELADC_FORTRAN_SOURCES
        ${PROJECT_SOURCE_DIR}/external/stieltjes/stieltjes_ql_diag.f90
        ${PROJECT_SOURCE_DIR}/external/stieltjes/stieltjesmod.f90
        ${FREE_RELADC_FORTRAN_SOURCES}
   )
   set(C_RELADC_SOURCES
        ${PROJECT_SOURCE_DIR}/external/stieltjes/herm_spline.c
   )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_RELADC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_RELADC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    reladc
    OBJECT
    ${FREE_RELADC_FORTRAN_SOURCES}
    ${FIXED_RELADC_FORTRAN_SOURCES}
    ${C_RELADC_SOURCES}
)

# add intermodule dependencies
add_dependencies(reladc gp)
if (ENABLE_STIELTJES)
    add_dependencies(reladc stieltjes)
endif()
