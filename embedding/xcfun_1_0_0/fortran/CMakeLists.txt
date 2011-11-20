
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(example_fortran_1_0_0
	example.F90
	xcfun_module.F90
	xcfun_autogen.F90
	)

set_target_properties(example_fortran_1_0_0
	PROPERTIES LINKER_LANGUAGE Fortran
	)

target_link_libraries(example_fortran_1_0_0 xcfun_1_0_0 m)
