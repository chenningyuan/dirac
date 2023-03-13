add_library(qcorr)

target_sources(qcorr
  PRIVATE
    qcorr_common_block_interface.F90
    qcorr_main.F90
    qcorr_cfg.F90
  )

target_link_libraries(qcorr
  PRIVATE
    luciarel
  )

target_include_directories(qcorr
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )
