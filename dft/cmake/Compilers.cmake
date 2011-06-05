# Fortran compilers

set(FORTRAN_COMPILER_VERSION "")

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-g -fcray-pointer -fbacktrace -DVAR_GFORTRAN -fno-range-check")
                                   # -fcray-pointer is for VAR_MPI2
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -ftest-coverage"
            )
    endif()

    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        execute_process(
            COMMAND ./get-gfortran-version
            TIMEOUT 1
            OUTPUT_VARIABLE FORTRAN_COMPILER_VERSION
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    set(CMAKE_Fortran_FLAGS         "-fno-second-underscore -ftrace=full -g -DVAR_G95")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS}"
            )
    endif()

    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        execute_process(
            COMMAND ./get-g95-version
            TIMEOUT 1
            OUTPUT_VARIABLE FORTRAN_COMPILER_VERSION
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl -DVAR_IFORT -g -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion -ftrapuv"
            )
    endif()

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()

    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        execute_process(
            COMMAND ./get-ifort-version
            TIMEOUT 1
            OUTPUT_VARIABLE FORTRAN_COMPILER_VERSION
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()

    set_source_files_properties(${FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree"
        )
    set_source_files_properties(${FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
endif()