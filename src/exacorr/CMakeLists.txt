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
   mbauto/talsh_auto_ccd_spinor.F90
   mbauto/talsh_auto_ccsd_spinor.F90
   mbauto/talsh_auto_ccsdt_spinor.F90
   mbauto/talsh_auto_ccsdt-1a_spinor.F90
   mbauto/talsh_auto_ccsdt-1b_spinor.F90
   mbauto/talsh_auto_ccsdt-2_spinor.F90
   mbauto/talsh_auto_ccsdt-3_spinor.F90
   mbauto/talsh_auto_ccsdtq_spinor.F90
   mbauto/talsh_auto_ccsdtqp_spinor.F90
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
   talsh_tenpi_ccd.F90
   talsh_tenpi_ccsd.F90
   talsh_tenpi_ccsdt.F90
   talsh_tenpi_ccsdtq.F90
   tapp_tenpi_lambda_ccd.F90
   tapp_tenpi_lambda_ccsd.F90
   tapp_tenpi_lambda_ccsdt.F90
   tapp_tenpi_lambda_ccsdtq.F90
   talsh_tenpi_lambda_ccd.F90
   talsh_tenpi_lambda_ccsd.F90
   talsh_tenpi_lambda_ccsdt.F90
   talsh_tenpi_lambda_ccsdtq.F90
   talsh_tenpi_1rdm_ccd.F90
   tapp_tenpi_1rdm_ccd.F90
   tapp_tenpi_1rdm_ccsd.F90
   tapp_tenpi_1rdm_ccsdt.F90
   tapp_tenpi_1rdm_ccsdtq.F90
   talsh_tenpi_1rdm_ccsd.F90
   talsh_tenpi_1rdm_ccsdt.F90
   talsh_tenpi_1rdm_ccsdtq.F90
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
   talsh_cholesky_routines.F90
   talsh_ao_to_mo.F90
   talsh_cholesky_integrals.F90
   talsh_common_routines.F90
   talsh_sigma_eomccsd.F90
   talsh_trial_eomccsd.F90
   talsh_drivers_eomccsd.F90
   talsh_density_matrices.F90
   talsh_1rdm_2rdm_ccsd_kernels.F90
   talsh_triples.F90
   talsh_cc.F90
   talsh_cis.F90
   talsh_cis_density.F90
   talsh_cis_property.F90
   talsh_ccrsp.F90
   talsh_ccrsp_common.F90
   talsh_eom_QR.F90
   talsh_cc_QR.F90
   talsh_ccci_QR.F90
   talsh_eom_LR.F90
   talsh_eom_rspeq_solver.F90
   talsh_QR_TPA.F90
   davidson/talsh_davidson.F90
   talsh_mp2lap.F90
   talsh_mp2no.F90
   mp2no_driver.F90
   mp2no_matrixtransforms.F90
   mp2no_fileinterfaces.F90
   talsh_rt_solver.F90
  )

  if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_EXACORR_FORTRAN_SOURCES}  EXACORR COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_EXACORR_FORTRAN_SOURCES} EXACORR COMPILE_FLAGS "-qfixed")
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
