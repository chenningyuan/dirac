message("-- DIRAC: Fortran compiler is  \"${CMAKE_Fortran_COMPILER}\" (\"${CMAKE_Fortran_COMPILER_ID}\") ")
message("-- DIRAC:       C compiler is  \"${CMAKE_C_COMPILER}\" (\"${CMAKE_C_COMPILER_ID}\") ") 

# Fortran compilers

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    add_definitions(-DVAR_GFORTRAN)
                                     # this is for VAR_MPI2
    set(CMAKE_Fortran_FLAGS         "-fcray-pointer -fbacktrace")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    if(DIRAC_ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(DIRAC_ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check"
            )
    endif()
    if(DIRAC_ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -ftest-coverage"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    add_definitions(-DVAR_G95)
    set(CMAKE_Fortran_FLAGS         "-fno-second-underscore -ftrace=full")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char")
    if(DIRAC_ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(DIRAC_ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
            )
    endif()
    if(DIRAC_ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS}"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFC)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    if(DIRAC_ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(DIRAC_ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion -traceback"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    add_definitions(-DVAR_PGF90)
    set(CMAKE_Fortran_FLAGS         " ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(DIRAC_ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(DIRAC_ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(DIRAC_ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

    set_source_files_properties(${DIRAC_FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree"
        )
    set_source_files_properties(${DIRAC_FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
endif()

# C compilers

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-Wall")
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_C_FLAGS_RELEASE "-g -O2 -Wno-unused")
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-g -O2")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE " ")
endif()

