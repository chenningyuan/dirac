execute_process(
    COMMAND cat ${CMAKE_SOURCE_DIR}/VERSION
    OUTPUT_VARIABLE DIRAC_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# reset GIT_REVISION
set(GIT_REVISION)

# if GIT_HASH exists then this is exported code
# in this case we read git hash from this file and set DEVELOPMENT_CODE to false
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/GIT_HASH")
    execute_process(
        COMMAND cat ${CMAKE_SOURCE_DIR}/cmake/GIT_HASH
        OUTPUT_VARIABLE GIT_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(DEVELOPMENT_CODE FALSE)
else()
    set(DEVELOPMENT_CODE TRUE)
endif()
