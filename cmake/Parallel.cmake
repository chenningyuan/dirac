if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})

        if(ENABLE_MPI2_DETECTION)
            # test whether MPI-2 is available
            exec_program(
                ${cat_program}
                ARGS ${CMAKE_SOURCE_DIR}/cmake/MPI-2-test.F90
                OUTPUT_VARIABLE mpi2test
                )
            check_fortran_source_compiles(
                ${mpi2test}
                MPI_2_SUPPORT_OK
                )
            if(MPI_2_SUPPORT_OK)
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
