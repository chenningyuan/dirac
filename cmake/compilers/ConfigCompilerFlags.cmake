# This variables are not guaranteed to be defined for all compilers or languages
if(NOT DEFINED CMAKE_C_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED CMAKE_CXX_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_CXX_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED CMAKE_Fortran_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_Fortran_COMPILER_ID variable is not defined! (CMake Error)")
endif()

# save flags to cache macro
include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS AND CMAKE_CXX_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)   
    # check if support for xHost flag is available (Intel)
    if(CMAKE_C_COMPILER_ID MATCHES "Intel")
        # compile test program
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
            execute_process(
                COMMAND ifort /QxHost ${CMAKE_SOURCE_DIR}/cmake/compilers/Test_xHost_Flag.F90 /exe:${CMAKE_BINARY_DIR}/test_xhost_flag.exe
            )
        else()
            execute_process(
                COMMAND ifort -xHost ${CMAKE_SOURCE_DIR}/cmake/compilers/Test_xHost_Flag.F90 -o ${CMAKE_BINARY_DIR}/test_xhost_flag.out
            )
        endif()

        # try to run test program
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
            execute_process(
                COMMAND ${CMAKE_BINARY_DIR}/test_xhost_flag.exe

                OUTPUT_VARIABLE TEST
                ERROR_VARIABLE TEST
            )
        else()
            execute_process(
                COMMAND ${CMAKE_BINARY_DIR}/test_xhost_flag.out

                OUTPUT_VARIABLE TEST
                ERROR_VARIABLE TEST
            )
        endif()

        # depends on result define variable or not
        if(TEST MATCHES "Intel xHost flag is available on this computer.")
            set(XHOST_AVAILABLE ON)
        endif()
    endif()
    
    # set flags
    include(CFlags)
	include(CXXFlags)
	include(FortranFlags)
	
	# check for correct compiler version with Python script and CMake variable (if defined)
	include(FindCompilerVersion)
endif()
