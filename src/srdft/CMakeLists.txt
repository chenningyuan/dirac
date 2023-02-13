add_library(srdft)

target_sources(srdft
  PRIVATE
    srdft_cfg.F90
    srdft.F90
  )

target_include_directories(srdft
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(srdft
  PRIVATE
    gp
    xcint
    grid
  )
