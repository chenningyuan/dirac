add_subdirectory(exatensor)

add_library(exacorr)

target_sources(exacorr
  PRIVATE
    exacorr_dirac.F
    exacorr_main.F90
    exacorr_datatypes.F90
    exacorr_global.F90
    exacorr_cfg.F90
    exacorr_utils.F90
    exacorr_eri.F90
    exacorr_tensor_methods.F90
    intermediates.F90
    exacorr_ao_to_mo.F90
    exacorr_gradient.F90
    exacorr_mo.F90
    exacorr_respect.F90
    exacorr_cc.F90
    exacorr_mp2no.F90
    talsh_cholesky_routines.F90
    talsh_ao_to_mo.F90
    talsh_cholesky_integrals.F90
    talsh_common_routines.F90
    talsh_gradient.F90
    talsh_triples.F90
    talsh_cc.F90
    talsh_mp2lap.F90
    talsh_mp2no.F90
    mp2no_driver.F90
    mp2no_matrixtransforms.F90
    mp2no_fileinterfaces.F90
  )

target_include_directories(exacorr
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/include
  )

target_link_libraries(exacorr
  PRIVATE
    gp
    interest
  )

add_dependencies(exacorr exatensor)
