include_directories(${PROJECT_BINARY_DIR})

add_custom_command(
    OUTPUT
    ${PROJECT_BINARY_DIR}/fortran_interface.h
    DEPENDS
    gen_fortran_interface_h
    COMMAND
    gen_fortran_interface_h > ${PROJECT_BINARY_DIR}/fortran_interface.h
    )

#
#  miro: replace utils/gen_fortran_interface_h.c with utils/gen_fortran_interface_h.F90
#  to avoid possible Intel -mkl flag collision with the non-Intel C-compiler for gen_fortran_interface_h.c
#
add_library(gfih_c ${PROJECT_SOURCE_DIR}/utils/gen_fortran_interface_h.c) # create library first
add_executable(
    gen_fortran_interface_h 
    ${PROJECT_SOURCE_DIR}/utils/gen_fortran_interface_h.F90)
target_link_libraries(gen_fortran_interface_h gfih_c) # secure dependency


set(_generated_files
    ${PROJECT_BINARY_DIR}/fortran_interface.h
    )

add_custom_target(
    dirac_interfaces
    ALL
    DEPENDS
    ${_generated_files}
    )


unset(_generated_files)
