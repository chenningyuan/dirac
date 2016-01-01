option(ENABLE_QCORR "Davidson-type +Q corrections for (MR)CISD" ON)

if(ENABLE_QCORR)
    add_definitions(-DMOD_QCORR)
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/qcorr/qcorr_common_block_interface.F90
        src/qcorr/qcorr_main.F90
        src/qcorr/qcorr_cfg.F90
        )
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_QCORR)
endif()

message(STATUS "Davidson-type +Q corrections for (MR)CISD: ${ENABLE_QCORR}")
