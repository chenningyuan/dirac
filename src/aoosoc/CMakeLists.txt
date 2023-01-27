include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libaoosoc.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_AOOSOC_FORTRAN_SOURCES
)

set(FREE_AOOSOC_FORTRAN_SOURCES
     atomic_oo_order_so_correction.F90
     module_aoosoc_pfg.F90
     atomic_oo_order_so_correction_cfg.F90
     atomic_oo_order_so_correction_utils.F90
     module_aoosoc_init.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_AOOSOC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_AOOSOC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'aoosoc' module without exceptions")
    set_source_files_properties(${FREE_AOOSOC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_AOOSOC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    aoosoc
    OBJECT
    ${FIXED_AOOSOC_FORTRAN_SOURCES}
    ${FREE_AOOSOC_FORTRAN_SOURCES}
)

#
# deal with intedependencies
#
add_dependencies(aoosoc interest)
add_dependencies(aoosoc x2c)
