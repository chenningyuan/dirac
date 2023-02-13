add_library(pelib)

target_sources(pelib
  PRIVATE
    src/pe_mpi.F90
    src/pe_precision.F90
    src/pe_constants.F90
    src/pe_variables.F90
    src/pe_utils.F90
    src/pe_operators.F90
    src/pe_multipole_tools.F90
    src/pe_analysis_tools.F90
    src/pe_potential_derivatives.F90
    src/pe_induced_moments.F90
    src/interfaces/pe_blas_interfaces.F90
    src/interfaces/pe_lapack_interfaces.F90
    src/interfaces/pe_integral_interfaces.F90
    src/polarizable_embedding.F90
  )

target_link_libraries(pelib
  PRIVATE
    gen1int_interface
  )

target_compile_definitions(pelib
  PRIVATE
    GEN1INT_INTEGRALS
    PRG_DIRAC
  )
