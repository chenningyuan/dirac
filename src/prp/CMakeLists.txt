add_library(prp)

target_sources(prp
  PRIVATE
    excprp.F
    mp2lag.F
    pamexp.F
    pamfck.F
    pamgrd.F
    pammag.F
    pamogrd.F
    pamprp.F
    pamrsp.F
    pamrvc.F
    pamset.F
    pamstex.F
    pamnqcc.F
    pamtpa.F
    pamxlr.F
    pamxpp.F
    rspmpg.F
    van_der_waals.F
    bed_tm.F
    bed_osc.F90
    conf_parameters.F90
    london_helper.F90
    london_utils.F90
    london_direct.F90
    london_reorth.F90
    paminp.F90
    pamxqr.F90
)

target_link_libraries(prp
  PRIVATE
    esr
)

target_include_directories(prp
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

option(ENABLE_LAO_REARRANGED "LAO properties without connection matrices" ON)

if(ENABLE_LAO_REARRANGED)
  target_compile_definitions(prp PRIVATE MOD_LAO_REARRANGED)

  dirac_test(lao_shield_rearranged "long" "")
endif()
