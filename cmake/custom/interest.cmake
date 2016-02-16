cmake_dependent_option(ENABLE_INTEREST "Enable Interest" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_INTEREST)
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
    add_definitions(-DMOD_INTEREST)
#miro: interest library might need some math routines from "${EXTERNAL_LIBS}"
    set(EXTERNAL_LIBS interest ${EXTERNAL_LIBS})
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/interest-interface/module_interest_interface.F90
        )
    dirac_test(mdirac_closed-shell_dft "interest")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_INTEREST)
endif()

message(STATUS "Interest library: ${ENABLE_INTEREST}")
