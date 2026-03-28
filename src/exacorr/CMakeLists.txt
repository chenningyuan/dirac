if(ENABLE_EXATENSOR)
  include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    )

  set(FIXED_EXACORR_FORTRAN_SOURCES
    exacorr_dirac.F
  )

  if (EXATENSOR_TALSH)
   set(FREE_EXACORR_FORTRAN_SOURCES exatensor_placeholders.F90)
   message(STATUS "ExaTensor not available or requested: activating placeholder")
  endif()

  set(FREE_EXACORR_FORTRAN_SOURCES
   ${FREE_EXACORR_FORTRAN_SOURCES}
   mbauto/eti_auto_ccd_spinor.F90
   mbauto/eti_auto_ccsd_spinor.F90
   mbauto/eti_auto_ccsdt_spinor.F90
   mbauto/eti_auto_ccsdt-1a_spinor.F90
   mbauto/eti_auto_ccsdt-1b_spinor.F90
   mbauto/eti_auto_ccsdt-2_spinor.F90
   mbauto/eti_auto_ccsdt-3_spinor.F90
   mbauto/eti_auto_ccsdtq_spinor.F90
   mbauto/eti_auto_ccsdtqp_spinor.F90
   mbauto/exatensor_auto_ccd_spinor.F90
   mbauto/exatensor_auto_ccsd_spinor.F90
   mbauto/exatensor_auto_ccsdt_spinor.F90
   mbauto/exatensor_auto_ccsdtq_spinor.F90
   mbauto/exatensor_auto_ccsdtqp_spinor.F90
   tapp.F90
   tapp_top.F90
   exacorr_tensor_interface.F90
   tapp_tenpi_ccd.F90
   tapp_tenpi_ccsd.F90
   tapp_tenpi_ccsdt.F90
   tapp_tenpi_ccsdtq.F90
   eti_tenpi_ccd.F90
   eti_tenpi_ccsd.F90
   eti_tenpi_ccsdt.F90
   eti_tenpi_ccsdtq.F90
   tapp_tenpi_lambda_ccd.F90
   tapp_tenpi_lambda_ccsd.F90
   tapp_tenpi_lambda_ccsdt.F90
   tapp_tenpi_lambda_ccsdtq.F90
   eti_tenpi_lambda_ccd.F90
   eti_tenpi_lambda_ccsd.F90
   eti_tenpi_lambda_ccsdt.F90
   eti_tenpi_lambda_ccsdtq.F90
   eti_tenpi_1rdm_ccd.F90
   tapp_tenpi_1rdm_ccd.F90
   tapp_tenpi_1rdm_ccsd.F90
   tapp_tenpi_1rdm_ccsdt.F90
   tapp_tenpi_1rdm_ccsdtq.F90
   eti_tenpi_1rdm_ccsd.F90
   eti_tenpi_1rdm_ccsdt.F90
   eti_tenpi_1rdm_ccsdtq.F90
   exatensor_tenpi_ccd.F90
   exatensor_tenpi_ccsd.F90
   exatensor_tenpi_ccsdt.F90
   exatensor_tenpi_ccsdtq.F90
   exacorr_datatypes.F90
   exacorr_1el.F90
   exacorr_global.F90
   exacorr_cfg.F90
   exacorr_utils.F90
   exacorr_eri.F90
   exacorr_tensor_methods.F90
   intermediates_cc.F90
   intermediates_lambda.F90
   intermediates_cc2.F90
   intermediates_eom.F90
   intermediates_eom_rsp.F90
   exacorr_ao.F90
   exacorr_ao_to_mo.F90
   exacorr_gradient.F90
   exacorr_mo.F90
   exacorr_respect.F90
   exacorr_cc.F90
   exacorr_mp2no.F90
   exacorr_printout.F90
   exacorr_scaling_parameters.F90
   eti_cholesky_routines.F90
   eti_ao_to_mo.F90
   eti_integrals.F90
   eti_cholesky_integrals.F90
   eti_common_routines.F90
   eti_sigma_eomccsd.F90
   eti_trial_eomccsd.F90
   eti_drivers_eomccsd.F90
   eti_density_matrices.F90
   eti_1rdm_2rdm_ccsd_kernels.F90
   eti_1rdm_2rdm_automatic_kernels.F90
   eti_drivers_gs_ccsd.F90
   eti_drivers_gs_automatic_kernels.F90
   eti_triples.F90
   eti_cc.F90
   eti_cis.F90
   eti_cis_density.F90
   eti_cis_property.F90
   eti_ccrsp.F90
   eti_ccrsp_common.F90
   eti_eom_QR.F90
   eti_cc_QR.F90
   eti_ccci_QR.F90
   eti_eom_LR.F90
   eti_eom_rspeq_solver.F90
   eti_QR_TPA.F90
   solvers/eti_davidson.F90
   solvers/eti_diis.F90
   eti_mp2lap.F90
   eti_mp2no.F90
   mp2no_driver.F90
   mp2no_matrixtransforms.F90
   mp2no_fileinterfaces.F90
   eti_rt_solver.F90
  )

  if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_EXACORR_FORTRAN_SOURCES}  EXACORR COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_EXACORR_FORTRAN_SOURCES} EXACORR COMPILE_FLAGS "-qfixed")
  endif()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
    set_source_files_properties(${PROJECT_SOURCE_DIR}/src/exacorr/tapp_top.F90 PROPERTY COMPILE_FLAGS "-O0")
    message(STATUS "src/exacorr/tapp_top.F90: disabling optimization for Intel OneAPI compiler")
  endif()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "22")
    if (CMAKE_BUILD_TYPE MATCHES "release") 
      set_source_files_properties(${PROJECT_SOURCE_DIR}/src/exacorr/exacorr_cc.F90 PROPERTY COMPILE_FLAGS "-O0")
      message(STATUS "src/exacorr/exacorr_cc.F90: disabling inteprocedural optimization for Intel compiler")
    endif()
  endif()


  add_library(
    exacorr
    STATIC
    ${FIXED_EXACORR_FORTRAN_SOURCES}
    ${FREE_EXACORR_FORTRAN_SOURCES}
  )

# deal with dependencies
  add_dependencies(exacorr  dirac_main)
  add_dependencies(exacorr  exatensor)
  add_dependencies(exacorr  interest)
  add_dependencies(exacorr  laplace) # uses module type_laplace
#  add_dependencies(exacorr  mfsolver)
endif()
