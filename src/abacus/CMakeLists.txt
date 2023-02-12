add_library(hergp)

target_sources(hergp
  PRIVATE
    hergp.F
  )

target_include_directories(hergp
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

#-------------------------------------------------------------------------------

add_library(abacus)

target_sources(abacus
  PRIVATE
    abadip.F
    abalnr.F
    abamolden.F
    abanuc.F
    abaop2.F
    abaopt.F
    abaortvec.F
    abarint.F
    abatro.F
    abavib.F
    abavrml.F
    abawalk.F
    her1ave.F
    her1car.F
    her1cat.F
    her1drv.F
    her1int.F
    her1odc.F
    her1out.F
    her1pro.F
    her1sym.F
    her2ave.F
    her2dir.F
    her2drv.F
    her2el1.F
    her2el2.F
    her2fck.F
    her2gab.F
    her2her.F
    her2odc.F
    her2out.F
    her2sym.F
    herbas.F
    herdirect.F
    herdrv.F
    hergam.F
    hergeo.F
    hergp.F
    hergroup.F
    hermem.F
    herrdn.F
    herrdn_dirac.F
    hersol.F
    hersym.F
    hertst.F
    vibnor.F
    herpar.F
    carsph.F90
    hermodules.F90
    filemodule.F90
    gab_settings.F90
    aobasis_collector.F90
  )

target_include_directories(abacus
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(abacus
  PRIVATE
    ecp
    gen1int_interface
  )

# FIXME
add_dependencies(abacus pelib)
add_dependencies(abacus pcmsolver)
