include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/main
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdirac_main.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_MAIN_FORTRAN_SOURCES
   dirac.F
)

set(FREE_MAIN_FORTRAN_SOURCES
   logo.F90
   mpi_framework.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_MAIN_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_MAIN_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'main' module without exceptions")
    set_source_files_properties(${FREE_MAIN_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_MAIN_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    dirac_main
    STATIC
    ${FIXED_MAIN_FORTRAN_SOURCES}
    ${FREE_MAIN_FORTRAN_SOURCES}
)

# deal with intermodule dependencies
add_dependencies(dirac_main gp)
add_dependencies(dirac_main visual)
add_dependencies(dirac_main fde)
add_dependencies(dirac_main interest)
if (ENABLE_STIELTJES)
    add_dependencies(dirac_main stieltjes)
endif()

add_dependencies(dirac_main laplace-minimax)

if(ENABLE_PCMSOLVER)
    add_dependencies(dirac_main pcm_dirac)
endif()
if(ENABLE_PELIB)
    add_dependencies(dirac_main pelib_interface)
endif()
add_dependencies(dirac_main interface_mpi)
