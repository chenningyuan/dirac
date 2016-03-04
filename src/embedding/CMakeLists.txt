include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libembedding.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_EMBEDDING_FORTRAN_SOURCES
  emb_quit.F
  get_vemb_mat.F
  # emb_dirac_integration.F90 needs GTDOAV routine :
  ${PROJECT_SOURCE_DIR}/src/dirac/dftaux.F
)

set(FREE_EMBEDDING_FORTRAN_SOURCES
   get_vemb_mat.F
   emb_quit.F
   emb_mpi.F90
   emb_mag_cfg.F90
   emb_dirac_integration.F90
   emb_types.F90
   emb_nadd_derv.F90
   emb_export.F90
   emb_io.F90
   emb_input.F90
   emb_dirac_input.F90
   emb_dirac_evaluators.F90
   emb_data.F90
   emb_cfg.F90
   emb_xcfun_interface.F90
   emb_max_block_length.F90
   emb_main.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_EMBEDDING_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_EMBEDDING_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

add_library(
   embedding 
   ${FIXED_EMBEDDING_FORTRAN_SOURCES}
   ${FREE_EMBEDDING_FORTRAN_SOURCES}
)

# add intermodule dependencies
add_dependencies(embedding xcfun)
add_dependencies(embedding gp)
