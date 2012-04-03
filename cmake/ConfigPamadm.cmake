if(EXISTS ${CMAKE_SOURCE_DIR}/pamadm)
    # check whether all pamadm-generated files exist
    if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/SourcesList.cmake"      AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/DefinitionsList.cmake"  AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/IncludeFilesList.cmake" AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/TestsList.cmake"        AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags")
        # they exist, read previously used flags
        execute_process(
            COMMAND cat ${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags
            OUTPUT_VARIABLE PAMADM_FLAGS
            )
    else()
        # set flags to '-a'
        set(PAMADM_FLAGS "-a")
    endif()

    execute_process(
        COMMAND python pamadm ${PAMADM_FLAGS}
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_QUIET
    )

    add_custom_target(
        run_pamadm
        COMMAND python pamadm ${PAMADM_FLAGS} > /dev/null
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )

    message("-- Packaging sources and tests using pamadm ${PAMADM_FLAGS}")
else()
    # pamadm does not exist
    # this must then be released code
    # in this case don't bother
endif()
