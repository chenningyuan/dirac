include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For liblocalization.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_LOCALIZATION_FORTRAN_SOURCES
)

set(FREE_LOCALIZATION_FORTRAN_SOURCES
   func_PipekMezey.F90
   localization_interface.F90
   second_order_minimization.f90
   select_coeff_own_bas.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_LOCALIZATION_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_LOCALIZATION_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    localization
    OBJECT
    ${FIXED_LOCALIZATION_FORTRAN_SOURCES}
    ${FREE_LOCALIZATION_FORTRAN_SOURCES}
)

add_dependencies(localization gp)
