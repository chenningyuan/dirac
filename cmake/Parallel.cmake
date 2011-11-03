if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
#miro: print out OpenMPI integer setting first, do sanity check aftewards
        configure_file(
            ${CMAKE_SOURCE_DIR}/cmake/ompi_integer_info.py.in
            ${CMAKE_BINARY_DIR}/ompi_integer_info.py)
        execute_process(
            COMMAND python ${CMAKE_BINARY_DIR}/ompi_integer_info.py
            TIMEOUT 1
            OUTPUT_VARIABLE OMPI_INTEGER_INFO_OUTPUT
            ERROR_VARIABLE  OMPI_INTEGER_INFO_STDERR
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (VERBOSE_OUTPUT)
            if (OMPI_INTEGER_INFO_OUTPUT)
                message("-- Open MPI Fortran ${OMPI_INTEGER_INFO_OUTPUT}")
            endif()
            if (OMPI_INTEGER_INFO_STDERR)
                message("stderr from ompi_integer_info.py.in:${OMPI_INTEGER_INFO_STDERR}")
            endif()
        endif()

        if(ENABLE_SANITY_CHECKS)
            # test whether MPI integer type matches
            exec_program(
                ${cat_program}
                ARGS ${CMAKE_SOURCE_DIR}/cmake/test-MPI-itype-compatibility.F90
                OUTPUT_VARIABLE _source
                )
            check_fortran_source_compiles(
                ${_source}
                MPI_ITYPE_MATCHES
                )
            if(MPI_ITYPE_MATCHES)
            else()
                if(ENABLE_64BIT_INTEGERS)
                   message(
                       FATAL_ERROR
                       "ERROR: You want to compile DIRAC with 64bit integers but your MPI is built with 32bit integers"
                       )
                else()
                   message(
                       FATAL_ERROR
                       "ERROR: You want to compile DIRAC with 32bit integers but your MPI is built with 64bit integers"
                       )
                endif()
            endif()
        endif()

        if(ENABLE_MPI2_DETECTION)
            # test whether MPI-2 is available
            exec_program(
                ${cat_program}
                ARGS ${CMAKE_SOURCE_DIR}/cmake/test-MPI-2-compatibility.F90
                OUTPUT_VARIABLE _source
                )
            check_fortran_source_compiles(
                ${_source}
                MPI_2_COMPATIBLE
                )
            if(MPI_2_COMPATIBLE)
                add_definitions(-DVAR_MPI2)
                message("-- MPI-2 support found")
            else()
                message("-- no MPI-2 support found, will try with MPI-1")
            endif()
        else()
            message("-- MPI-2 support disabled, will try with MPI-1")
        endif()

    endif()
endif()

if(ENABLE_OMP)
    find_package(OpenMP)
    if(OpenMP_FOUND)
        add_definitions(-DVAR_OMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_COMPILE_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()
