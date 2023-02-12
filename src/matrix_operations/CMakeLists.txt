add_library(matrix_operations)

target_sources(matrix_operations
  PRIVATE
    dirac_interface.F90
    matrix_defop.F90
    matrix_genop.F90
    matrix_module.F90
    matrix_operations.F90
    matrix_stat.F90
    memory_bouncer.F90
  )

target_include_directories(matrix_operations
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(matrix_operations
  PRIVATE
    gp
    grid
  )
