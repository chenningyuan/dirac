cmake_dependent_option(ENABLE_SRDFT "srDFT module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_SRDFT)
    add_definitions(-DMOD_SRDFT)
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/srdft/srdft.F90
        src/srdft/srdft_cfg.F90
        )
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_SRDFT)
endif()

message(STATUS "srDFT module: ${ENABLE_SRDFT}")
