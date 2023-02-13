add_library(laplace)

target_sources(laplace
  PRIVATE
    pamlaplce.F90
    lap_driver.F90
  )

target_include_directories(laplace
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(laplace
  PRIVATE
    laplace-minimax
    gp
  )
