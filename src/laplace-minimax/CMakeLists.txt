add_library(laplace-minimax)

target_sources(laplace-minimax
  PRIVATE
    src/dd128_arithmetics.F90
    src/dd128_linalg.F90
    src/lap_rderror.F90
    src/lap_rddata.F90
    src/lap_paraopt.F90
    src/lap_maehly.F90
    src/lap_rmsd.F90
    src/lap_numlap.F90
    src/laplace_minimax.F90
  )

target_include_directories(laplace-minimax
  PRIVATE
    inc
  )

file(COPY data/init_para.txt  DESTINATION ${PROJECT_BINARY_DIR}/laplace-minimax/data)
file(COPY data/init_error.txt DESTINATION ${PROJECT_BINARY_DIR}/laplace-minimax/data)

install(
  FILES
    data/init_para.txt
  DESTINATION
    share/dirac/laplace-minimax/data
  PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
  )

install(
  FILES
    data/init_error.txt
  DESTINATION
    share/dirac/laplace-minimax/data
  PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
  )
