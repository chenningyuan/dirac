if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcray-pointer -fbacktrace -fno-range-check -ffree-line-length-none -DVAR_GFORTRAN -DVAR_MFDS")

    if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "9.9.9")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()

    if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -maix64")
        message(STATUS "Added gfortran -maix64 flag due to IBM AIX and integer*8")
    endif()

    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -w")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -fbounds-check -finit-real=nan -finit-integer=-9999  -finit-logical=true  -Wall")
    #set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif()
