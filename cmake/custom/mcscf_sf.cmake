cmake_dependent_option(ENABLE_MCSCF_spinfree "Spinfree MCSCF module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_MCSCF_spinfree)
    add_definitions(-DMOD_MCSCF_spinfree)
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_MCSCF_spinfree)
endif()

message(STATUS "Spinfree MCSCF module: ${ENABLE_MCSCF_spinfree}")
