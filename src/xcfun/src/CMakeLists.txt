
include_directories(taylor functionals)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/include)

set (BASE_LIB_SRCS 
	xcfun.cpp
	evaluators.cpp
	functional.cpp
	parameters.cpp
	potential.cpp
	fortran.c
	)

file (GLOB FUNCTIONAL_SRCS functionals/*.cpp) 

add_library(xcfun
    ${BASE_LIB_SRCS}
    ${FUNCTIONAL_SRCS}
    )

add_custom_command (
	OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/include
		#${CMAKE_CURRENT_BINARY_DIR}/fortran # cmake bug
	COMMAND ${CMAKE_COMMAND} -E make_directory 
		${CMAKE_CURRENT_BINARY_DIR}/include
	COMMAND ${CMAKE_COMMAND} -E make_directory 
		${CMAKE_CURRENT_BINARY_DIR}/fortran
	)

install(TARGETS xcfun ARCHIVE DESTINATION lib)
