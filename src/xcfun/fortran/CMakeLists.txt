
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(example_fortran 
	example.F90
	xcfun_module.F90
	xcfun_autogen.F90
	)

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#   miro: necessary for PGI
    target_link_libraries(example_fortran stdc++)
    message("added libstdc++ to linker flags of xcfun due to PGI ")
endif()

set_target_properties(example_fortran 
	PROPERTIES LINKER_LANGUAGE Fortran
	)

target_link_libraries(example_fortran xcfun m)
