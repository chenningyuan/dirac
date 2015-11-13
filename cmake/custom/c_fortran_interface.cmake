include_directories(${PROJECT_BINARY_DIR})

add_custom_command(
    OUTPUT
    ${PROJECT_BINARY_DIR}/fortran_interface.h
    DEPENDS
    gen_fortran_interface_h
    COMMAND
    gen_fortran_interface_h > ${PROJECT_BINARY_DIR}/fortran_interface.h
    )

add_executable(
    gen_fortran_interface_h
    ${PROJECT_SOURCE_DIR}/utils/gen_fortran_interface_h.c
    )

set(_generated_files
    ${PROJECT_BINARY_DIR}/fortran_interface.h
    )

add_custom_target(
    dirac_interfaces
    ALL
    DEPENDS
    ${_generated_files}
    )

add_dependencies(dirac dirac_interfaces)

unset(_generated_files)
