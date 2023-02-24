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
    ../include
    ../cfun
  )

# FIXME this will be refactored later
add_dependencies(grid xcint)
add_dependencies(grid gp)
