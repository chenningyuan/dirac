include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libluciarel.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_LUCIAREL_FORTRAN_SOURCES
cands.inc  
 cintfo.inc 
 crun.inc  
  files_r.inc 
 implicit.inc   
    memory_r.inc 
  oper.inc   
 strinp.inc
cgas.inc   
 clunit.inc 
 cstate.inc  
gasstr.inc 
  integrals_off.inc 
 mxpdim.inc   
  stinf.inc 
  wrkspc.inc
cicisp.inc 
 cprnt.inc 
  ctcc.inc  
  glbbas.inc 
  ipoist8.inc    
    noccn_inf.inc 
 strbas.inc
   ciinfo_r.F
   diagonal_r.F
   diagonalize_r.F
   gasci_main_driver_r.F
   interface_r.F
   io_r.F
   memory_r.F
   par_lucia_r.F
   par_lushmem_r.F
   par_solver_r.F
   property_r.F
   resolution_utils_r.F
   sigden_r.F
   strings_r.F
   utils_r.F
   vecmat_r.F
)

set(FREE_LUCIAREL_FORTRAN_SOURCES
   io_ttss_r.F90
   krci_cfg.F90
   krci_prop.F90
   mospinor_info.F90
   string_info.F90
   symmetry_setup.F90
)

set(LUCIAREL_C_SOURCES
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_LUCIAREL_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_LUCIAREL_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    luciarel
    ${FIXED_LUCIAREL_FORTRAN_SOURCES}
    ${FREE_LUCIAREL_FORTRAN_SOURCES}
    ${LUCIAREL_C_SOURCES}
)

# deal with intermodule dependencies
add_dependencies(luciarel gp)
