add_library(visual)

target_sources(visual
  PRIVATE
    integration.F90
    visual_cfg.F90
    visual.F90
    visual_in_point.F90
    visual_london.F90
  )

target_include_directories(visual
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(visual
  PRIVATE
    gp
    matrix_operations
  )
