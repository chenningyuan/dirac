include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_BINARY_DIR}/laplace-minimax/src/laplace-minimax/inc
)

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdummy.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

# source files
set(FREE_F95_DUMMY_FORTRAN_SOURCES
   dummy_fixed.F90
)

set(FREE_F03_DUMMY_FORTRAN_SOURCES
   dummy_free.F90
)

# compile with 2003 standard
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set_source_files_properties(${FREE_F03_DUMMY_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-std=f2003")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set_source_files_properties(${FREE_F03_DUMMY_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-stand f03")
endif()

# check XL
if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_DUMMY_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_DUMMY_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'dummy' module without exceptions")
    set_source_files_properties(${FREE_DUMMY_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_DUMMY_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()
add_library(
    dummy
    OBJECT
    ${FREE_F95_DUMMY_FORTRAN_SOURCES}
    ${FREE_F03_DUMMY_FORTRAN_SOURCES}
)

if (ENABLE_LAPLACE)
    add_dependencies(dummy laplace-minimax)
endif()
