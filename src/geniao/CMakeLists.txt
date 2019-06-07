if(ENABLE_GENIAO)
  include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

  set(FREE_GENIAO_FORTRAN_SOURCES
   generate_iao.F90
   get_overlap.F90
   make_iao.F90
   utils.F90
   atomic_fragments.F90
   read_xyz.F90
   write_basis.F
   make_avas.F90
  )

  if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_GENIAO_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfree=f90")
  endif()

  add_library(
    geniao
    OBJECT
    ${FREE_GENIAO_FORTRAN_SOURCES}
  )

# deal with dependencies
  add_dependencies(geniao exacorr)
endif()
