add_library(grid)

target_sources(grid
  PRIVATE
    num_grid_cfg.F90
    dftgrd.F
    num_grid_gen.F90
    grtest.F
    Lebedev-Laikov.F
  )

target_include_directories(grid
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
  )

target_link_libraries(grid
  PRIVATE
    gp
    xcint
  )
