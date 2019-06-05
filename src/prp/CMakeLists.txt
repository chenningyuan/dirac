include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libprp.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_PRP_FORTRAN_SOURCES
   dcbexp.h  
   dcbstex.h  
   dcbxpp.h  
   dcbxqr.h  
   gauss_legendre.h  
   ibndxdef.h 
   thrldp.h
   excprp.F
   mp2lag.F
   pamexp.F
   pamfck.F
   pamgrd.F
   pammag.F
   pamogrd.F
   pamprp.F
   pamrsp.F
   pamset.F
   pamtpa.F
   pamxlr.F
   rspmpg.F
   van_der_waals.F
   bed_tm.F
   bed_osc.F90
   # fix missing rmolgrd
   ${PROJECT_SOURCE_DIR}/src/dirac/dirgrd.F
)

#Miro: volunteers are welcome to fix runtimecheck flaws !
set(FIXED_PRP_FORTRAN_SOURCES_noruntimecheck
   pamrvc.F
   pamxpp.F
   pamstex.F
)

set(FREE_PRP_FORTRAN_SOURCES
   conf_parameters.F90
   london_helper.F90
   london_utils.F90
   london_direct.F90
   paminp.F90
)

#Miro: volunteers are welcome to fix runtimecheck flaws !
set(FREE_PRP_FORTRAN_SOURCES_noruntimecheck
   london_reorth.F90
   pamxqr.F90 # Attempt to fetch from allocatable variable VECAEP when it is not allocated
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_PRP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FREE_PRP_FORTRAN_SOURCES_noruntimecheck}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_PRP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${FIXED_PRP_FORTRAN_SOURCES_noruntimecheck} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

# enable runtimechecking flags only to SELECTED source-files
if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'prp' module, with some source files excluded from checking")
    set_source_files_properties(${FREE_PRP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_PRP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    prp
    OBJECT
    ${FREE_PRP_FORTRAN_SOURCES}
    ${FREE_PRP_FORTRAN_SOURCES_noruntimecheck}
    ${FIXED_PRP_FORTRAN_SOURCES}
    ${FIXED_PRP_FORTRAN_SOURCES_noruntimecheck}
)

# deal with dependencies
add_dependencies(prp  gp)
add_dependencies(prp  dirac)
add_dependencies(prp  main)
add_dependencies(prp  krmc)
add_dependencies(prp  relccsd)
if (ENABLE_ESR)
    add_dependencies(prp  esr)
endif()
if (ENABLE_PCMSOLVER)
    add_dependencies(prp  pcm_dirac)
endif()
