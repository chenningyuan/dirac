# get compiler version with Python script
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/cmake/compilers/getFortranCompilerVersion.py ${CMAKE_Fortran_COMPILER_ID} 
				OUTPUT_VARIABLE PYTHON_VERSION)

if(DEFINED CMAKE_Fortran_COMPILER_VERSION)
	SET(CMAKE_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum reqiured version
	if(NOT ((CMAKE_VERSION VERSION_EQUAL PYTHON_VERSION) AND 
	(CMAKE_VERSION VERSION_GREATER MIN_Fortran_COMPILER_VERSION OR CMAKE_VERSION VERSION_EQUAL MIN_Fortran_COMPILER_VERSION)))
		
		MESSAGE(FATAL_ERROR "Not supported version of Fortran compiler.")
		
	endif(NOT ((CMAKE_VERSION VERSION_EQUAL PYTHON_VERSION) AND 
	(CMAKE_VERSION VERSION_GREATER MIN_Fortran_COMPILER_VERSION OR CMAKE_VERSION VERSION_EQUAL MIN_Fortran_COMPILER_VERSION)))
	
else(DEFINED CMAKE_Fortran_COMPILER_VERSION)
	
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	if(NOT (PYTHON_VERSION VERSION_GREATER MIN_Fortran_COMPILER_VERSION OR PYTHON_VERSION VERSION_EQUAL MIN_Fortran_COMPILER_VERSION))
	
		MESSAGE(FATAL_ERROR "Not supported version of Fortran compiler.")
	
	endif(NOT (PYTHON_VERSION VERSION_GREATER MIN_Fortran_COMPILER_VERSION OR PYTHON_VERSION VERSION_EQUAL MIN_Fortran_COMPILER_VERSION))
	
endif(DEFINED CMAKE_Fortran_COMPILER_VERSION)
