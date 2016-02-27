include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdensfit.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_DENSFIT_FORTRAN_SOURCES
   aobtch.h 
   denfit_noddy.F 
   df_dirac.F 
   df_fck3.F 
   df_utils.F 
   df_vxc.F 
   iteration.h
)

set(FREE_DENSFIT_FORTRAN_SOURCES
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_DENSFIT_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_DENSFIT_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
   densfit
   ${FIXED_DENSFIT_FORTRAN_SOURCES}
   ${FREE_DENSFIT_FORTRAN_SOURCES}
)

add_dependencies(densfit xcint)
