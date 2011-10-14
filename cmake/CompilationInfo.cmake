set(COMPILATION_INFO_FILE ${CMAKE_SOURCE_DIR}/gp/compilation_info.F90)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/compilation_info.py.in
    ${CMAKE_BINARY_DIR}/compilation_info.py
    )
execute_process(
    COMMAND python ${CMAKE_BINARY_DIR}/compilation_info.py
    TIMEOUT 1
    OUTPUT_VARIABLE RABOOF
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
#message("testing return: ${RABOOF}")

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/git_revision_info.py.in
    ${CMAKE_BINARY_DIR}/git_revision_info.py
    )
execute_process(
    COMMAND python ${CMAKE_BINARY_DIR}/git_revision_info.py
    TIMEOUT 1
    OUTPUT_FILE ${CMAKE_SOURCE_DIR}/gp/git_revision_info.F90
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# get revision info for cmake output
# radovan: that's code doubling!!!!!! i don't understand why
execute_process(
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/git_revision_info_cmake.py
    TIMEOUT 1
    OUTPUT_VARIABLE GIT_REVISION
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

