include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libpcm_dirac.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_PCM_DIRAC_FORTRAN_SOURCES
)

set(FREE_PCM_DIRAC_FORTRAN_SOURCES
    pcm_gp.F90
    pcm_integrals.F90
    pcm_linear_response.F90
    pcm_scf.F90
    pcm_write.F90
    pcmmod_cfg.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_PCM_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_PCM_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    pcm_dirac
    ${FIXED_PCM_DIRAC_FORTRAN_SOURCES}
    ${FREE_PCM_DIRAC_FORTRAN_SOURCES}
)

add_dependencies(pcm_dirac interface_mpi)
add_dependencies(pcm_dirac pcmsolver)
