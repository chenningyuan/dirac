# check whether pamadm-generated files exist
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/SourcesList.cmake"      AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/DefinitionsList.cmake"  AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/IncludeFilesList.cmake" AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/TestsList.cmake"        AND
   EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags")
    # they exist, read flags
    execute_process(
        COMMAND cat ${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags
        OUTPUT_VARIABLE PAMADM_FLAGS
        )
else()
    # they do not exist: run pamadm -a
    set(PAMADM_FLAGS "-a")
endif()

execute_process(
    COMMAND python pamadm ${PAMADM_FLAGS}
    TIMEOUT 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_QUIET
)

message("-- Packaging sources and tests using pamadm ${PAMADM_FLAGS}")
