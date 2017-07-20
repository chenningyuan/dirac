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
# autocmake.yml configuration::
#
#   docopt: "--check Enable runtime checking of compiled executables [default: False]"
#   define: "'-DENABLE_RUNTIMECHECK={0}'.format(arguments['--check'])"

option(ENABLE_RUNTIMECHECK "Enable run-time checking" OFF)

message(STATUS "Enable run-time checking: ${ENABLE_RUNTIMECHECK}")

if(ENABLE_RUNTIMECHECK)
    if(DEFINED CMAKE_Fortran_COMPILER_ID)
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -fbounds-check")
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -check all -fpstkchk")
        endif()
    endif()

    if(DEFINED CMAKE_C_COMPILER_ID)
        if(CMAKE_C_COMPILER_ID MATCHES GNU)
            set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -fbounds-check")
        endif()
        if(CMAKE_C_COMPILER_ID MATCHES Intel)
            set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -fp-stack-check -check-uninit")
        endif()
    endif()

    if(DEFINED CMAKE_CXX_COMPILER_ID)
        if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
            set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fbounds-check")
        endif()
        if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
            set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fp-stack-check -check-uninit -early-template-check")
        endif()
    endif()
endif()
