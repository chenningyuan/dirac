cmake_dependent_option(ENABLE_LAO_REARRANGED "LAO properties without connection matrices" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_LAO_REARRANGED)
    add_definitions(-DMOD_LAO_REARRANGED)
    dirac_test(lao_shield_rearranged "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_LAO_REARRANGED)
endif()

message(STATUS "LAO properties without connection matrices: ${ENABLE_LAO_REARRANGED}")
