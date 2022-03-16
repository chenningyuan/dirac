include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libx-amfi.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_XAMFI_FORTRAN_SOURCES
)

set(FREE_XAMFI_FORTRAN_SOURCES
     x_fock.F90
     amfDataStorage.F90
     amfH5Interface.F90
     xamfi_utils.F90
     xamfi_driver.F90
     xamfi_environment.F90
     xamfi_global_parameters.F90
     xamfi_internal_parameters.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_XAMFI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_XAMFI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    x-amfi
    OBJECT
    ${FIXED_XAMFI_FORTRAN_SOURCES}
    ${FREE_XAMFI_FORTRAN_SOURCES}
)

#
# deal with intedependencies
#
if (ENABLE_INTEREST)
    add_dependencies(x-amfi interest)
endif()
add_dependencies(x-amfi x2c)
add_dependencies(x-amfi grid)
