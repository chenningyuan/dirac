cmake_dependent_option(ENABLE_ESR "ESR property module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_ESR)
    add_definitions(-DMOD_ESR)
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        src/prp/esr/pamesr.F
        src/prp/esr/pamesrci.F
        )
    dirac_test(esr "short;gascip;ci")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_ESR)
endif()

message(STATUS "ESR property module: ${ENABLE_ESR}")
