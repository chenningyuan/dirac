add_library(stieltjes)

target_sources(stieltjes
  PRIVATE
    stieltjes_ql_diag.f90
    stieltjesmod.f90
    driver.f90
    input_routines.f90
    herm_spline.c
  )

# FIXME what to do with the test
