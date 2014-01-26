# get compiler version from Python script for C, C++ and Fortran (C;C++;Fortran -- form like CMake list)
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/cmake/compilers/getCompilerVersion.py 
				${CMAKE_C_COMPILER_ID} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_Fortran_COMPILER_ID}
				OUTPUT_VARIABLE PYTHON_OUTPUT)
				
# get versions from python output
list(GET PYTHON_OUTPUT 0 PYTHON_C_VERSION)
list(GET PYTHON_OUTPUT 1 PYTHON_CXX_VERSION)
list(GET PYTHON_OUTPUT 2 PYTHON_Fortran_VERSION)

# cut off \n added by Python print()
string(STRIP ${PYTHON_Fortran_VERSION} PYTHON_Fortran_VERSION)

# what is minimum required version?
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_C ${MIN_GCC_COMPILER_VERSION})
endif()
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_C ${MIN_Intel_C_COMPILER_VERSION})
endif()
if(CMAKE_C_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_C ${MIN_PGI_C_COMPILER_VERSION})
endif()
if(CMAKE_C_COMPILER_ID MATCHES "XL")
	set(REQUIRED_C ${MIN_XL_C_COMPILER_VERSION})
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_CXX ${MIN_GXX_COMPILER_VERSION})
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_CXX ${MIN_Intel_CXX_COMPILER_VERSION})
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_CXX ${MIN_PGI_CXX_COMPILER_VERSION})
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "XL")
	set(REQUIRED_CXX ${MIN_XL_CXX_COMPILER_VERSION})
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_Fortran ${MIN_GFORTRAN_COMPILER_VERSION})
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_Fortran ${MIN_Intel_FORTRAN_COMPILER_VERSION})
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_Fortran ${MIN_PGI_FORTRAN_COMPILER_VERSION})
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
	set(REQUIRED_Fortran ${MIN_XL_FORTRAN_COMPILER_VERSION})
endif()

# checking for correct version of C compiler
if(DEFINED CMAKE_C_COMPILER_VERSION)
	set(CMAKE_C_VERSION ${CMAKE_C_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT ((CMAKE_C_VERSION VERSION_EQUAL PYTHON_C_VERSION) AND 
	(CMAKE_C_VERSION VERSION_GREATER REQUIRED_C OR CMAKE_C_VERSION VERSION_EQUAL REQUIRED_C)))
		MESSAGE(FATAL_ERROR "Not supported version of C compiler. Found: ${CMAKE_C_VERSION}, Required: ${REQUIRED_C}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_C_VERSION VERSION_GREATER REQUIRED_C OR PYTHON_C_VERSION VERSION_EQUAL REQUIRED_C))
		MESSAGE(FATAL_ERROR "Not supported version of C compiler. Found: ${PYTHON_C_VERSION}, Required: ${REQUIRED_C}")
	endif()
endif()

# checking for correct version of C++ compiler
if(DEFINED CMAKE_CXX_COMPILER_VERSION)
	set(CMAKE_CXX_VERSION ${CMAKE_CXX_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT ((CMAKE_CXX_VERSION VERSION_EQUAL PYTHON_CXX_VERSION) AND 
	(CMAKE_CXX_VERSION VERSION_GREATER REQUIRED_CXX OR CMAKE_CXX_VERSION VERSION_EQUAL REQUIRED_CXX)))
		MESSAGE(FATAL_ERROR "Not supported version of C++ compiler. Found: ${CMAKE_CXX_VERSION}, Required: ${REQUIRED_CXX}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_CXX_VERSION VERSION_GREATER REQUIRED_CXX OR PYTHON_CXX_VERSION VERSION_EQUAL REQUIRED_CXX))
		MESSAGE(FATAL_ERROR "Not supported version of C++ compiler. Found: ${PYTHON_CXX_VERSION}, Required: ${REQUIRED_CXX}")
	endif()
endif()

# checking for correct version of Fortran compiler
if(DEFINED CMAKE_Fortran_COMPILER_VERSION)
	set(CMAKE_Fortran_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT ((CMAKE_Fortran_VERSION VERSION_EQUAL PYTHON_Fortran_VERSION) AND 
	(CMAKE_Fortran_VERSION VERSION_GREATER REQUIRED_Fortran OR CMAKE_Fortran_VERSION VERSION_EQUAL REQUIRED_Fortran)))
		MESSAGE(FATAL_ERROR 
		"Not supported version of Fortran compiler. Found: ${CMAKE_Fortran_VERSION}, Required: ${REQUIRED_Fortran}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_Fortran_VERSION VERSION_GREATER REQUIRED_Fortran OR PYTHON_Fortran_VERSION VERSION_EQUAL REQUIRED_Fortran))
		MESSAGE(FATAL_ERROR 
		"Not supported version of Fortran compiler. Found: ${PYTHON_Fortran_VERSION}, Required: ${REQUIRED_Fortran}")
	endif()
endif()
