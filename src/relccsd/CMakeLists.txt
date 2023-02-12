add_library(relccsd)

target_sources(relccsd
  PRIVATE
    cceqn_diis.F
    cceqn_driver_amplitudes.F
    cceqn_driver_initialization.F
    cceqn_mp2_amplitudes.F
    cceqn_t1_amplitudes.F
    cceqn_t1_utils.F
    cceqn_t2_amplitudes.F
    cceqn_t3_energy_correction.F
    cceqn_utils.F
    ccfspc.F
    ccgenp.F
    ccgetv.F
    ccintm.F
    ccints.F
    ccpam.F
    ccrpa.F
    ccxblas.F
    waio.F
    ccdriv.F
    cceqn_driver_lambda.F
    ccgrad.F
    relcc_cfg.F90
    spinor_indexing.F90
    auto_offset.F90
    ccsort.F90
    ccsort_exp.F90
    wrapper_sorting.F90
    contraction.F90
    intermediates.F90
    lambda_equation.F90
    sigma_eom.F90
    davidson/davidson.F90
    trial_vectors_eom.F90
    density.F90
    projectors.F90
    eom_driver.F90
)

target_include_directories(relccsd
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(relccsd
  PRIVATE
    prp
  )

#-------------------------------------------------------------------------------

# FIXME does not link and dependencies should be pushed one level down
option(ENABLE_RELCCSD_STANDALONE "Enable standalone RELCCSD" OFF)

if(ENABLE_RELCCSD_STANDALONE)
  add_executable(relccsd.x)

  target_sources(relccsd.x
    PRIVATE
      ccmain.F
    )

  target_include_directories(relccsd.x
    PRIVATE
      ${PROJECT_SOURCE_DIR}/src/include
    )

  target_link_libraries(relccsd.x
    PRIVATE
      abacus
      amfi
      cfun
      dirac
      dirrci
      ecp
      embedding
      esr
      gp
      grid
      input
      interest
      interest
      interface_mpi
      krmc
      laplace
      localization
      luciarel
      lucita
      main
      matrix_operations
      moltra
      prp
      qcorr
      reladc
      relccsd
      visual
      x2c
      xcfun_fortran_bindings
      xcint
    )
endif()
