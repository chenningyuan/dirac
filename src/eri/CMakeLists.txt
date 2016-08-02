include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For liberi.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_ERI_FORTRAN_SOURCES
       eri2car2.F
       eri2odb.F
       fck3mpi.F
       eri2vec.F
       fck3int.F
       eri2her.F
       eri2csi.F
       eri2aob.F
       eri2fck.F
       eri2pro.F
       eri2sph.F
       eri2drv.F
       eri2cft.F
       eri2car1.F
       fck3drv.F
       fck3aob.F
       fck3ftuv.F
       eri2r12.F
       eri2ave.F
       eri2out.F
       fck3odb.F
       eri2ftuv.F
       fck3box.F
       eri2fck3.F
       fck3mom.F
       erisph.F
       eri2dis.F
       eri2sym.F
       eri2exc.F
)

set(FREE_ERI_FORTRAN_SOURCES
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_ERI_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_ERI_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
    eri
    OBJECT
    ${FIXED_ERI_FORTRAN_SOURCES}
    ${FREE_ERI_FORTRAN_SOURCES}
)

#
# deal with intermodule dependencies 
# 
add_dependencies(eri gp)
