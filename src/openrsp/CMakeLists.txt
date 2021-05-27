include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libopenrsp.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_OPENRSP_FORTRAN_SOURCES
    VIBCTL_interface.F
)

set(FREE_OPENRSP_FORTRAN_SOURCES
     openrsp_interface_response.F90
     vib_prop.F90
     openrsp_interface_2el.F90
     prop_contribs.F90
     fc.F90
     openrsp_cfg.F90
     openrsp_interface_1el.F90
     prop_test.F90
     rsp_equations.F90
     birefring.F90
     openrsp_output.F90
     openrsp_main.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_OPENRSP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_OPENRSP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'openrsp' module without exceptions")
    set_source_files_properties(${FREE_OPENRSP_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_OPENRSP_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    openrsp
    OBJECT
    ${FIXED_OPENRSP_FORTRAN_SOURCES}
    ${FREE_OPENRSP_FORTRAN_SOURCES}
)

# add intermodule dependencies
#add_dependencies(openrsp gp)

