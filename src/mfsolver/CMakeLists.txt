include_directories(
  ${PROJECT_SOURCE_DIR}/src/include
  )

set(FREE_EXACORR_FORTRAN_SOURCES
 mfsolver_cfg.F90
 mfsolver_dirac_input.F90
)

add_library(
  mfsolver
  STATIC
  ${FREE_EXACORR_FORTRAN_SOURCES}
)

add_dependencies(mfsolver interface_mpi)
