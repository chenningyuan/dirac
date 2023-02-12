add_library(interface_mpi)

target_sources(interface_mpi
  PRIVATE
    interface_to_mpi_f77.F
    integer_kind_mpilib.F90
    integer_model.F90
    interface_to_mpi.F90
  )

target_include_directories(interface_mpi
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

# target_link_libraries(interface_mpi
#   PRIVATE
#     gp
#     xcint
#     visual
#     embedding
#   )
