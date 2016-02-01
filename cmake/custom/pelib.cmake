option(ENABLE_PELIB "Enable polarizabile embedding library" ON)

if(ENABLE_PELIB)
    include_directories(${PROJECT_BINARY_DIR}/external/pelib-build)

    add_definitions(-DBUILD_PELIB)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	-DINTEGRAL_LIBRARY=GEN1INT
	-DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    )

    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
	src/peqm/pe_dirac_interfaces.F
        )

    add_external(pelib)

    add_dependencies(pelib gen1int)

    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libpelib.a
        ${EXTERNAL_LIBS})
endif()

message(STATUS "PElib module: ${ENABLE_PELIB}")
