if(VERBOSE_OUTPUT)
    execute_process(
        COMMAND python ${CMAKE_BINARY_DIR}/compilation_info.py
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND python ${CMAKE_BINARY_DIR}/git_revision_info.py
        TIMEOUT 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()
