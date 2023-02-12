add_subdirectory(gen1int)

add_library(gen1int_interface)

target_sources(gen1int_interface
  PRIVATE
    gen1int_api.F90
    gen1int_cube.F90
    gen1int_host.F90
    gen1int_matrix.F90
    gen1int_shell.F90
  )

target_compile_definitions(gen1int_interface
  PRIVATE
    PRG_DIRAC
  )

target_include_directories(gen1int_interface
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(gen1int_interface
  PRIVATE
    gen1int
  )
