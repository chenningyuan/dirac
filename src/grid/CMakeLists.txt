include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libgrid.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_GRID_FORTRAN_SOURCES
    dftwrk.h
    Lebedev-Laikov.F
    dftgrd.F
    grtest.F
)

set(FREE_GRID_FORTRAN_SOURCES
  num_grid_cfg.F90
  num_grid_gen.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_GRID_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_GRID_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    grid
    ${FIXED_GRID_FORTRAN_SOURCES}
    ${FREE_GRID_FORTRAN_SOURCES}
)

#
# deal with intedependencies
#
## modules/interface_grid.mod
add_dependencies(grid xcint)
## modules/dirac_cfg.mod
add_dependencies(grid gp)
