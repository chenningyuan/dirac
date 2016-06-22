include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For liblucita.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_LUCITA_FORTRAN_SOURCES
blksize.inc 
 cgas.inc  
  cprnt.inc 
  cstate.inc  
gasstr.inc  
 irat.inc   
  multd2h.inc  
orbinf_lucita.inc 
 stinf.inc  
 units.inc
cands.inc  
  cicisp.inc  
crun.inc  
  files.inc  
 glbbas.inc
   luci.inc   
  mxpdim.inc  
 orbinp.inc     
    strbas.inc
  vrsn.inc
cc_exc.inc 
  clunit.inc 
 csfbas.inc  
fio.inc    
 gshelp.inc  
 lucinp.inc  
 oper.inc   
  rasdim.inc   
      strinp.inc  
wrkspc.inc
cecore.inc   
comjep.inc
  csm.inc   
  fnc.inc  
   intform.inc
  machine.inc 
 orbinf.inc  
 spinfo.inc    
     SysDef.inc
    ciinfo.F
    density.F
    diagonal.F
    doubledir.F
    fock.F
    interface.F
    io.F
    lucia_cc.F
    lucia_icci.F
    lucita_diagonalize.F
    memory.F
    par_lucia.F
    par_solver.F
    program.F
    property.F
    sigma.F
    strings.F
    transform.F
    utils.F
    vecmat.F
)

set(FREE_LUCITA_FORTRAN_SOURCES
)

set(LUCITA_C_SOURCES
 flush_stdout.c
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_LUCITA_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_LUCITA_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    lucita
    OBJECT
    ${FIXED_LUCITA_FORTRAN_SOURCES}
    ${FREE_LUCITA_FORTRAN_SOURCES}
    ${LUCITA_C_SOURCES}
)

add_dependencies(lucita qcorr)
