set(MATH_FOUND  FALSE)
set(MATH_LANG   "Fortran")

if(DIRAC_ENABLE_EXTERNAL_MATH)

    # user sets
    set(DIRAC_USERDEFINED_LIBS
        "${DIRAC_USERDEFINED_LIBS}"
        CACHE STRING
        "User set math libraries"
        FORCE
        )
    if(NOT "${DIRAC_USERDEFINED_LIBS}" STREQUAL "")
        set(MATH_LIBS
            "${DIRAC_USERDEFINED_LIBS}"
            CACHE STRING
            "User set math libraries"
            FORCE
            )
        message("-- DIRAC: User set math libraries: ${MATH_LIBS}")
        set(MATH_FOUND TRUE)
    endif()

    # try to find the best library using environment variables MATH_ROOT or
    # BLAS_ROOT and LAPACK_ROOT
    if(NOT MATH_FOUND)
        find_package(BLAS)
        find_package(LAPACK)
        if(BLAS_FOUND AND LAPACK_FOUND)
            set(MATH_LIBS
                ${BLAS_LIBRARIES}
                ${LAPACK_LIBRARIES}
                )
            set(MATH_FOUND TRUE)
        endif()
    endif()

    if(MATH_FOUND)
        set(LIBS
            ${LIBS}
            ${MATH_LIBS}
            )
    else()
        message("-- No external math libraries found - defaulting to DIRAC own math libraries")
    endif()
if (ENABLE_CUDA)
	find_package(CUDA)
endif ()

endif()

if(NOT MATH_FOUND)
    message("-- Using own BLAS and LAPACK math")
    add_definitions(-DUSE_DIRAC_BLAS)
    add_definitions(-DUSE_DIRAC_DLAMCH)
    add_definitions(-DUSE_DIRAC_DLAPACK)
    add_definitions(-DUSE_DIRAC_ZLAPACK)
endif()
