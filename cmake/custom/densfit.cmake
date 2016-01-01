cmake_dependent_option(ENABLE_DENSFIT "Density fitting module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_DENSFIT)
    add_definitions(-DMOD_DNF)
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        src/densfit/df_dofit.F90
        src/densfit/df_vxc.F
        src/densfit/df_utils.F
        src/densfit/df_fck3.F
        src/densfit/df_dirac.F
        )
    dirac_test(integrals_densfit "")
    dirac_test(integrals_jengine "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_DNF)
endif()

message(STATUS "Density fitting module: ${ENABLE_DENSFIT}")
