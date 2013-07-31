# only in the development code we get the hash from git
# in the released code we read it from file, in this case
# it is already set at this stage
if(DEVELOPMENT_CODE)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND          ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            OUTPUT_VARIABLE  GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND     ${GIT_EXECUTABLE} status -uno
            OUTPUT_FILE ${PROJECT_BINARY_DIR}/GIT_STATUS_AT_BUILD
            ERROR_QUIET
            )
        if(EXISTS "${PROJECT_BINARY_DIR}/GIT_STATUS_AT_BUILD")
            file(STRINGS ${PROJECT_BINARY_DIR}/GIT_STATUS_AT_BUILD GIT_STATUS LIMIT_COUNT 1)
            if(NOT ${GIT_STATUS} STREQUAL "")
                string(STRIP ${GIT_STATUS} GIT_STATUS)
            endif()
        endif()
    endif()
endif()
