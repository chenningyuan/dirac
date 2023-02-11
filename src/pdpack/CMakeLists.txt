add_library(pdpack)

target_sources(pdpack
  PRIVATE
    arhpack.F
    dge.F
    dlamch.F
    dsp.F
    eispack.F
    gp_auxlapack.F
    gp_dlapack.F
    gp_zlapack.F
    gpblas.F
    linextra.F
    printpkg.F
    gp_zlinpack.F
 )

target_include_directories(pdpack
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )
