if(VERBOSE_OUTPUT)
    # get compilation info to cmake output
    execute_process(
        COMMAND python ${CMAKE_BINARY_DIR}/compilation_info.py
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # get last git revision to cmake output
    execute_process(
        COMMAND cat ${CMAKE_SOURCE_DIR}/main/config.h
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()
