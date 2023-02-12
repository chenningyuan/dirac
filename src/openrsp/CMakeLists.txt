add_library(openrsp)

target_sources(openrsp
  PRIVATE
    VIBCTL_interface.F
    openrsp_interface_response.F90
    vib_prop.F90
    openrsp_interface_2el.F90
    prop_contribs.F90
    fc.F90
    openrsp_cfg.F90
    openrsp_interface_1el.F90
    prop_test.F90
    rsp_equations.F90
    birefring.F90
    openrsp_output.F90
    openrsp_main.F90
  )

target_include_directories(openrsp
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(openrsp
  PRIVATE
    xcint
    matrix_operations
  )

dirac_test(legacy_openrsp_cars "openrsp" "")
dirac_test(legacy_openrsp_cme "openrsp" "")
dirac_test(legacy_openrsp_delta "openrsp" "")
dirac_test(legacy_openrsp_efgb "openrsp" "")
dirac_test(legacy_openrsp_efgb_dft "openrsp" "")
dirac_test(legacy_openrsp_electric "openrsp" "")
dirac_test(legacy_openrsp_jones "openrsp" "")
dirac_test(legacy_openrsp_magnetizability "openrsp" "")
dirac_test(legacy_openrsp_molecular_gradient "openrsp" "")
dirac_test(legacy_openrsp_pv "openrsp" "")
dirac_test(legacy_openrsp_quadrupole "openrsp" "")
dirac_test(legacy_openrsp_roa "openrsp" "")
dirac_test(legacy_openrsp_pnc_gradient "openrsp" "")
