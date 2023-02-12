add_library(reladc)

target_sources(reladc
  PRIVATE
    adcdble.F
    adcdiag.F
    adcsort.F
    polprp_buildtm.F
    polprp_cousat.F
    polprp_driver.F
    polprp_helpers.F
    polprp_mpimaster.F
    polprp_mpiserver.F
    polprp_prpint.F
    adccore.F
    polprp_mainblk.F
    adc_mat.F90
    adc_cfg.F90
    adc_fano_exchange.F90
    adc_fano_diag.F90
    adc_fano_routines.F90
    adc_fano_complex_routines.F90
    adc_fano_real_routines.F90
    adc_fano_matmul.F90
    adcutil.F90
    adc_fano.F90
    polprp_cfg.F90
    polprp_qstackm.F90
    polprp_blockdav.F90
    polprp_blockdav_c.F90
    polprp_qstack.c
  )

target_include_directories(reladc
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(reladc
  PRIVATE
    gp
  )

if(CMAKE_Fortran_COMPILER_ID MATCHES XL OR CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  # xlf compiler does not support isnan used
  # pgfortran compiler does not suppore real*16
  option(ENABLE_STIELTJES "Enable stieltjes library" OFF)
else()
  option(ENABLE_STIELTJES "Enable stieltjes library" ON)
endif()

if(ENABLE_STIELTJES)
  add_subdirectory(stieltjes)

  target_compile_definitions(reladc PRIVATE -DHAS_STIELTJES)

  target_link_libraries(reladc
    PRIVATE
      stieltjes
    )
endif()
