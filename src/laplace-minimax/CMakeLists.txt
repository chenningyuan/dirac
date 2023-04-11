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

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  target_compile_options(laplace-minimax
    PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:${Fortran_FLAGS}> -O1 -H -ffree-form -Wall
    )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  target_compile_options(laplace-minimax
    PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:${Fortran_FLAGS}> -O3 -free -warn all -assume protect_parens -nogen-interfaces
    )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  target_compile_options(laplace-minimax
    PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:${Fortran_FLAGS}> -O4 -fast -fastsse -Mfree -Mlarge_arrays -Mcache_align -Msmart -Msmartalloc -Minform=inform
    )
endif()

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
