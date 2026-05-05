include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/dirac
)

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libdirac.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_DIRAC_FORTRAN_SOURCES
   atom2molecule_utils.F
   dirana.F
   diraos.F
   dirbss.F
   dircnv.F
   dirden.F
   dir_effQED.F
   dirgeo.F
   dirgp.F
   dirgrd.F
   dirmp2no.F
   dirone.F
   diropen.F
   dirout.F
   dirrdn.F
   dirscf.F
   dirset.F
   dirsol.F
   dirspf.F
   dirtra.F
   dirtwo.F
   dirv1c.F
   huckel.F
)

set(FREE_DIRAC_FORTRAN_SOURCES
   dirimpmos.F90
   # dirmp2no.F needs fortran module from relcc_cfg.F90
   # (we cannot use add_dependencies to relccsd, because relccsd also depends on dirac)
   ${PROJECT_SOURCE_DIR}/src/relccsd/relcc_cfg.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()


if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'dirac' module without exceptions")
    set_source_files_properties(${FREE_DIRAC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_DIRAC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    dirac
    STATIC
    ${FIXED_DIRAC_FORTRAN_SOURCES}
    ${FREE_DIRAC_FORTRAN_SOURCES}
)

#
# Oh, intermodule dependendencies, worth to reduce if possible
#
add_dependencies(dirac dirac_main)
add_dependencies(dirac abacus)
add_dependencies(dirac gp)
add_dependencies(dirac xcfun_fortran_bindings)
add_dependencies(dirac fde)
add_dependencies(dirac localization)
add_dependencies(dirac visual)
add_dependencies(dirac x2c)
add_dependencies(dirac interface_mpi)
add_dependencies(dirac interest)
if (ENABLE_PCMSOLVER)
    add_dependencies(dirac pcm_dirac)
endif()
if (ENABLE_XAMFI)
    add_dependencies(dirac x-amfi)
endif()
if (ENABLE_SRDFT)
    add_dependencies(dirac srdft)
endif()
add_dependencies(dirac reladc)

target_link_libraries(dirac PRIVATE qcorr)

if (ENABLE_KRCC)
    add_dependencies(dirac krcc)
endif()
if(ENABLE_EXACORR)
    add_dependencies(dirac exacorr)
endif()


