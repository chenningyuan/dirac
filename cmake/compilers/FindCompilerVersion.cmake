########################################################################################
#
# CMake script for getting (from Python script and from CMake variables if defined)
# and comparing compilers versions
#
# At this moment the CMAKE_<LANG>_COMPILER_ID variable must be defined for C, C++ and Fortran: 
# see CFlags.cmake, CXXFlags.cmake and FortranFlags.cmake, otherwise quit with fatal error.
#
# Written by Ivan Hrasko, February 2014.
# 
#########################################################################################

# What is minimum required compiler version depending on vendor? 
# Please add according to your experience. But if you want to specify more numbers 
# or different format than you see here then you must also make changes in GetCompilerVersion.py.
# (becasuse that script matches only as many numbers as you see here)
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_C "4.1.2")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_C "10.1")
elseif(CMAKE_C_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_C "12.5")
elseif(CMAKE_C_COMPILER_ID MATCHES "XL")
	set(REQUIRED_C "11.01")
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
	set(REQUIRED_C "3.1")
else()
	message(FATAL_ERROR "Vendor of your C compiler is not supported")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_CXX "4.1.2")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_CXX "10.1")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_CXX "12.5")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")
	set(REQUIRED_CXX "11.01")
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
	set(REQUIRED_CXX "3.1")
else()
	message(FATAL_ERROR "Vendor of your C++ compiler is not supported")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(REQUIRED_Fortran "4.3.0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(REQUIRED_Fortran "10.1")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
	set(REQUIRED_Fortran "12.5")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
	set(REQUIRED_Fortran "11.01")
else()
	message(FATAL_ERROR "Vendor of your Fortran compiler is not supported")
endif()

# get compiler version from Python script for C, C++ and Fortran compilers (C;C++;Fortran -- in form like CMake list)
# some compilers have different behaviour on Windows (Intel) or ...
execute_process(COMMAND python ${PROJECT_SOURCE_DIR}/cmake/compilers/GetCompilerVersion.py 
				${CMAKE_C_COMPILER_ID} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_SYSTEM_NAME}
				OUTPUT_VARIABLE PYTHON_OUTPUT)
				
# get versions from python output
list(GET PYTHON_OUTPUT 0 PYTHON_C_VERSION)
list(GET PYTHON_OUTPUT 1 PYTHON_CXX_VERSION)
list(GET PYTHON_OUTPUT 2 PYTHON_Fortran_VERSION)

# cut off \n added by Python print()
string(STRIP ${PYTHON_Fortran_VERSION} PYTHON_Fortran_VERSION)

# checking for correct version of C compiler
if(DEFINED CMAKE_C_COMPILER_VERSION)
	set(CMAKE_C_VERSION ${CMAKE_C_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT (CMAKE_C_VERSION VERSION_EQUAL PYTHON_C_VERSION))
                # Ivan, please fix this:
		#message("Versions from CMake variable and Python script for C compiler are not equal")
		#message("CMAKE_C_VERSION=${CMAKE_C_VERSION}   PYTHON_C_VERSION=${PYTHON_C_VERSION}")
	endif()
	
	if(NOT (CMAKE_C_VERSION VERSION_GREATER REQUIRED_C OR CMAKE_C_VERSION VERSION_EQUAL REQUIRED_C))
		message(FATAL_ERROR "Not supported version of C compiler. Found: ${CMAKE_C_VERSION}, Required: ${REQUIRED_C}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_C_VERSION VERSION_GREATER REQUIRED_C OR PYTHON_C_VERSION VERSION_EQUAL REQUIRED_C))
		message(FATAL_ERROR "Not supported version of C compiler. Found: ${PYTHON_C_VERSION}, Required: ${REQUIRED_C}")
	endif()
endif()

# checking for correct version of C++ compiler
if(DEFINED CMAKE_CXX_COMPILER_VERSION)
	set(CMAKE_CXX_VERSION ${CMAKE_CXX_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT (CMAKE_CXX_VERSION VERSION_EQUAL PYTHON_CXX_VERSION))
		#message(FATAL_ERROR "Versions from CMake variable and Python script for C++ compiler are not equal")
	endif()
	
	if(NOT (CMAKE_CXX_VERSION VERSION_GREATER REQUIRED_CXX OR CMAKE_CXX_VERSION VERSION_EQUAL REQUIRED_CXX))
		message(FATAL_ERROR "Not supported version of C++ compiler. Found: ${CMAKE_CXX_VERSION}, Required: ${REQUIRED_CXX}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_CXX_VERSION VERSION_GREATER REQUIRED_CXX OR PYTHON_CXX_VERSION VERSION_EQUAL REQUIRED_CXX))
		message(FATAL_ERROR "Not supported version of C++ compiler. Found: ${PYTHON_CXX_VERSION}, Required: ${REQUIRED_CXX}")
	endif()
endif()

# checking for correct version of Fortran compiler
if(DEFINED CMAKE_Fortran_COMPILER_VERSION)
	set(CMAKE_Fortran_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
	
	# if CMake variable is available then must be equal to version obtained from Python script
	# and equal or greater than minimum required version
	
	if(NOT (CMAKE_Fortran_VERSION VERSION_EQUAL PYTHON_Fortran_VERSION))
		#message(FATAL_ERROR "Versions from CMake variable and Python script for Fortran compiler are not equal")
	endif()
	
	if(NOT (CMAKE_Fortran_VERSION VERSION_GREATER REQUIRED_Fortran OR CMAKE_Fortran_VERSION VERSION_EQUAL REQUIRED_Fortran))
		message(FATAL_ERROR 
		"Not supported version of Fortran compiler. Found: ${CMAKE_Fortran_VERSION}, Required: ${REQUIRED_Fortran}")
	endif()
else()
	# if CMake variable is not defined, only compiler version from Python script must be equal or greater
	# than minimum required version
	
	if(NOT (PYTHON_Fortran_VERSION VERSION_GREATER REQUIRED_Fortran OR PYTHON_Fortran_VERSION VERSION_EQUAL REQUIRED_Fortran))
		message(FATAL_ERROR 
		"Not supported version of Fortran compiler. Found: ${PYTHON_Fortran_VERSION}, Required: ${REQUIRED_Fortran}")
	endif()
endif()
