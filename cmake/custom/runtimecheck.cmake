#.rst:
#
# Enables run-time checking by appending corresponding compiler flags.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS
#   CMAKE_C_FLAGS
#   CMAKE_CXX_FLAGS
#
# autocmake.cfg configuration::
#
#   docopt:  --check Enable runtime checking of compiled executables [default: False]
#   define: '-DENABLE_RUNTIMECHECK={0}'.format(arguments['--check'])

option(ENABLE_RUNTIMECHECK "Enable run-time checking" OFF)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS_RUNTIMECHECK "-fbounds-check")
    set(CMAKE_C_FLAGS_RUNTIMECHECK "-fbounds-check")
    set(CMAKE_CXX_FLAGS_RUNTIMECHECK "-fbounds-check ")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS_RUNTIMECHECK " -check all -fpstkchk")
    set(CMAKE_C_FLAGS_RUNTIMECHECK " -fp-stack-check -check-uninit")
    set(CMAKE_CXX_FLAGS_RUNTIMECHECK " -fp-stack-check -check-uninit -early-template-check")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS_RUNTIMECHECK " ")
    set(CMAKE_C_FLAGS_RUNTIMECHECK " ")
    set(CMAKE_CXX_FLAGS_RUNTIMECHECK " ")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS_RUNTIMECHECK " ")
    set(CMAKE_C_FLAGS_RUNTIMECHECK " ")
    set(CMAKE_CXX_FLAGS_RUNTIMECHECK " ")
endif()

if(ENABLE_RUNTIMECHECK)
    message(STATUS "Run-time checking through compiler flags ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RUNTIMECHECK}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RUNTIMECHECK}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RUNTIMECHECK}")
endif()
