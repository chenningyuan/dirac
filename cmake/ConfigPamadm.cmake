if(EXISTS ${CMAKE_SOURCE_DIR}/pamadm)
    # check whether all pamadm-generated files exist
    if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/SourcesList.cmake"      AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/DefinitionsList.cmake"  AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/IncludeFilesList.cmake" AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/TestsList.cmake"        AND
       EXISTS "${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags")
       # they exist, read previously used flags from "pamadm-flags"
       file(READ ${CMAKE_SOURCE_DIR}/cmake/pamadm-generated/pamadm-flags PAMADM_FLAGS)
       message("-- Read previous pamadm flags from file, these are: ${PAMADM_FLAGS}")
       # miro: inserted check for Windows
       if (PAMADM_FLAGS STREQUAL "") 
           message(FATAL_ERROR "No pamadm flags in the read file !")
       endif()
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
    #miro: tweak null device for Windows & Unix
    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set(NULL_DEVICE "/dev/null")
    else()
        set(NULL_DEVICE "NUL")
    endif()
	
    add_custom_target(
        run_pamadm
        COMMAND python pamadm ${PAMADM_FLAGS} > ${NULL_DEVICE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )

    message("-- Packaging sources and tests using pamadm ${PAMADM_FLAGS}")
else()
    # pamadm does not exist
    # this must then be released code
    # in this case don't bother
endif()
