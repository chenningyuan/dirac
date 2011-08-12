# doxygen documentation
configure_file(
    ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
    ${PROJECT_BINARY_DIR}/Doxyfile
)
add_custom_target(
    doc
    COMMAND doxygen
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
