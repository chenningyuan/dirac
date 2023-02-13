add_library(embedding)

target_sources(embedding
  PRIVATE
    emb_quit.F
    get_vemb_mat.F
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
    emb_dirac_evaluators.F90
  )

target_include_directories(embedding
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(embedding
  PRIVATE
    gp
    xcint
  )
