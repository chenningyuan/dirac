cmake_dependent_option(ENABLE_AOOSOC "Atomic oo-order spin-orbit correction module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_AOOSOC)
    add_definitions(-DMOD_AOOSOC)
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/aoosoc/atomic_oo_order_so_correction.F90
        src/aoosoc/module_aoosoc_pfg.F90
        src/aoosoc/atomic_oo_order_so_correction_cfg.F90
        src/aoosoc/atomic_oo_order_so_correction_utils.F90
        src/aoosoc/module_aoosoc_init.F90
        )
    dirac_test(atomic_x2c_2esoc_models_energy "medium")
    dirac_test(atomic_x2c_spinfree_energy "medium")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_AOOSOC)
endif()

message(STATUS "Atomic oo-order spin-orbit correction module: ${ENABLE_AOOSOC}")
