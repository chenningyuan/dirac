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

# macro for saving flags to cache
include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS AND CMAKE_CXX_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)

    # check whether we can use -xHost on this machine
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

        # compile test program
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
            set(_xhost_binary ${CMAKE_BINARY_DIR}/test_xhost_flag.exe)
            execute_process(COMMAND ifort /QxHost ${CMAKE_SOURCE_DIR}/cmake/compilers/test_xhost_flag.f90 /exe:${_xhost_binary})
        else()
            set(_xhost_binary ${CMAKE_BINARY_DIR}/test_xhost_flag.x)
            execute_process(COMMAND ifort -xHost ${CMAKE_SOURCE_DIR}/cmake/compilers/test_xhost_flag.f90 -o ${_xhost_binary})
        endif()

        # run test program
        execute_process(
            COMMAND ${_xhost_binary}
            OUTPUT_VARIABLE _xhost_binary_output
            ERROR_VARIABLE _xhost_binary_output
        )

        # depending on result set XHOST_FLAG_AVAILABLE
        if(_xhost_binary_output MATCHES "Intel xHost flag is available on this computer.")
            set(XHOST_FLAG_AVAILABLE ON)
        else()
            set(XHOST_FLAG_AVAILABLE OFF)
        endif()
    endif()

    # set flags
    include(CFlags)
    include(CXXFlags)
    include(FortranFlags)

    # check for correct compiler version with Python script and CMake variable (if defined)
    include(FindCompilerVersion)
endif()
