include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For librelccsd.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_RELCCSD_FORTRAN_SOURCES
   cceqn_diis.F
   cceqn_driver_amplitudes.F
   cceqn_driver_initialization.F
   cceqn_driver_lambda.F
   cceqn_mp2_amplitudes.F
   cceqn_t1_amplitudes.F
   cceqn_t1_utils.F
   cceqn_t2_amplitudes.F
   cceqn_t3_energy_correction.F
   cceqn_utils.F
   ccfspc.F
   ccgenp.F
   ccgetv.F
   ccgrad.F
   ccintm.F
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

set(FIXED_RELCCSD_FORTRAN_SOURCES_noruntimecheck
   ccints.F
   ccdriv.F # Attempt to fetch from allocatable variable A_HEFF when it is not allocated
)

set(FREE_RELCCSD_FORTRAN_SOURCES
   relcc_cfg.F90
   auto_offset.F90
   ccsort_exp.F90
   contraction.F90
   sigma_eom.F90
   davidson/davidson.F90
   density.F90
   relcc_cfg.F90
   trial_vectors_eom.F90
   eom_driver.F90 
   # needs dependency
   ${PROJECT_SOURCE_DIR}/src/main/mpi_framework.F90
# Miro: previously "problematic" source files wrt runtime check
   #intermediates.F90
   lambda_equation.F90
   #wrapper_sorting.F90 
)

set(FREE_RELCCSD_FORTRAN_SOURCES_noruntimecheck
  # Miro: indeed, I had to return the file back, unfortunately
   intermediates.F90
   wrapper_sorting.F90  
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES_noruntimecheck}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'relccsd' module with some exceptions")
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

#miro: all these manually picked up sources are necessary for building up the relccsd library, "make relccsd"
# - the quit.F file is not listed due to replacement by DIRAC QUIT routine
# - the ccunix.c file is not included
add_library(
    relccsd
    OBJECT
    ${FREE_RELCCSD_FORTRAN_SOURCES}
    ${FREE_RELCCSD_FORTRAN_SOURCES_noruntimecheck}
    ${FIXED_RELCCSD_FORTRAN_SOURCES}
    ${FIXED_RELCCSD_FORTRAN_SOURCES_noruntimecheck}
)

# add intermodule dependencies
add_dependencies(relccsd gp)
add_dependencies(relccsd dirac)
add_dependencies(relccsd main)
add_dependencies(relccsd interface_mpi)
