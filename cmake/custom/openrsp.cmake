option(ENABLE_OPENRSP "Enable OpenRSP" OFF)

if(ENABLE_OPENRSP)
    add_definitions(-DMOD_OPENRSP)
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/openrsp/openrsp_interface_response.F90
        src/openrsp/vib_prop.F90
        src/openrsp/openrsp_interface_2el.F90
        src/openrsp/prop_contribs.F90
        src/openrsp/fc.F90
        src/openrsp/openrsp_cfg.F90
        src/openrsp/openrsp_interface_1el.F90
        src/openrsp/prop_test.F90
        src/openrsp/rsp_equations.F90
        src/openrsp/birefring.F90
        src/openrsp/openrsp_output.F90
        src/openrsp/openrsp_main.F90
        )
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        src/openrsp/VIBCTL_interface.F
        )
    dirac_test(legacy_openrsp_cars "openrsp")
    dirac_test(legacy_openrsp_cme "openrsp")
    dirac_test(legacy_openrsp_delta "openrsp")
    dirac_test(legacy_openrsp_efgb "openrsp")
    dirac_test(legacy_openrsp_efgb_dft "openrsp")
    dirac_test(legacy_openrsp_electric "openrsp")
    dirac_test(legacy_openrsp_jones "openrsp")
    dirac_test(legacy_openrsp_magnetizability "openrsp")
    dirac_test(legacy_openrsp_molecular_gradient "openrsp")
    dirac_test(legacy_openrsp_pv "openrsp")
    dirac_test(legacy_openrsp_quadrupole "openrsp")
    dirac_test(legacy_openrsp_roa "openrsp")
    dirac_test(legacy_openrsp_pnc_gradient "openrsp")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_OPENRSP)
endif()

message(STATUS "OpenRSP library: ${ENABLE_OPENRSP}")
