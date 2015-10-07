if(MPI_FOUND)
    include(CheckFortranSourceCompiles)

    add_definitions(-DVAR_MPI)
    add_definitions(-DVAR_MPI2)

    # test whether MPI module is compatible with compiler
    file(READ "${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-compiler-compatibility.F90" _source)
    check_fortran_source_compiles(
        ${_source}
        MPI_COMPILER_MATCHES
        )
    if(MPI_COMPILER_MATCHES)
        message(STATUS "mpi.mod matches current compiler, setting -DUSE_MPI_MOD_F90")
        add_definitions(-DUSE_MPI_MOD_F90)
    else()
        message(STATUS "WARNING: mpi.mod compiled with different compiler")
    endif()

    # test whether MPI integer type matches
    file(READ "${PROJECT_SOURCE_DIR}/cmake/custom/test-MPI-itype-compatibility.F90" _source)
    check_fortran_source_compiles(
        ${_source}
        MPI_ITYPE_MATCHES
        )
    if(NOT MPI_ITYPE_MATCHES)
        if(ENABLE_64BIT_INTEGERS)
            message(STATUS "No 64-bit integer MPI interface found, will use 32-bit integer MPI interface")
            add_definitions(-DVAR_MPI_32BIT_INT)
            set(USE_32BIT_MPI_INTERFACE TRUE)
        else()
            message(STATUS "WARNING: Cannot determine whether MPI is built for 32bit integers")
        endif()
    endif()

    if(ENABLE_64BIT_INTEGERS AND FORCE_32BIT_MPI_INTERFACE)
        message(STATUS "32-bit integer MPI interface activated by the user")
        add_definitions(-DVAR_MPI_32BIT_INT)
        set(USE_32BIT_MPI_INTERFACE TRUE)
    endif()
endif()
