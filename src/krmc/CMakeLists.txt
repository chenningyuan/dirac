add_library(krmc)

target_sources(krmc
  PRIVATE
    krmc_gradient.F
    krmccan.F
    krmcgascip.F
    krmcgp.F
    krmcinp.F
    krmcmlf.F
    krmcopt.F
    krmcout.F
    krmcsav.F
    krmctra.F
    krmcvec.F
    krmcluci.F
    krci_ctl.F
    mc_energies.F90
    mcscf_routines.F90
    gendet.c
  )

target_include_directories(krmc
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(krmc
  PRIVATE
    luciarel
    qcorr
  )
