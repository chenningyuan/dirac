#.rst:
#
# Enables profiling by appending corresponding compiler flags.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS
#   CMAKE_C_FLAGS
#   CMAKE_CXX_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--prof Enable profiling [default: False]"
#   define: "'-DENABLE_PROFILING={0}'.format(arguments['--check'])"

option(ENABLE_PROFILING "Enable profiling" OFF)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS_PROFILE " -pg")
    set(CMAKE_C_FLAGS_PRFILE " -pg")
    set(CMAKE_CXX_FLAGS_PROFILE " -pg")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS_PROFILE " -pg")
    set(CMAKE_C_FLAGS_PROFILE " -pg")
    set(CMAKE_CXX_FLAGS_PROFILE " -pg")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS_PROFILE " -pg")
    set(CMAKE_C_FLAGS_PROFILE " -pg")
    set(CMAKE_CXX_FLAGS_PROFILE " -pg")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS_PROFILE " -pg")
    set(CMAKE_C_FLAGS_PROFILE " -pg")
    set(CMAKE_CXX_FLAGS_PROFILE " -pg")
endif()

if(ENABLE_PROFILING)
    message(STATUS "Code profiling through compiler flags ENABLED")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_PROFILE}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_PROFILE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_PROFILE}")
endif()
