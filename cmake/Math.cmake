set(MATH_LANG "Fortran")

set(BLAS_LIBRARIES)
set(LAPACK_LIBRARIES)

option(FORCE_OWN_BLAS   "Force use of own BLAS"   OFF)
option(FORCE_OWN_LAPACK "Force use of own LAPACK" OFF)
mark_as_advanced(FORCE_OWN_BLAS)
mark_as_advanced(FORCE_OWN_LAPACK)

if(NOT ENABLE_INTERNAL_MATH)

    # user sets
    set(USERDEFINED_MATH
        "${USERDEFINED_MATH}"
        CACHE STRING
        "User set math libraries"
        FORCE
        )
    if(NOT "${USERDEFINED_MATH}" STREQUAL "")
        set(MATH_LIBS
            "${USERDEFINED_MATH}"
            CACHE STRING
            "User set math libraries"
            FORCE
            )
        message("-- User set math libraries: ${MATH_LIBS}")
        set(BLAS_LIBRARIES   TRUE)
        set(LAPACK_LIBRARIES TRUE)
        set(LIBS
            ${LIBS}
            ${MATH_LIBS}
            )
    endif()

    # try to find the best libraries using environment variables
    if(NOT FORCE_OWN_BLAS)
        if(NOT BLAS_LIBRARIES)
            find_package(BLAS)
            if(BLAS_LIBRARIES)
                set(LIBS
                    ${LIBS}
                    ${BLAS_LIBRARIES}
                    )
            endif()
        endif()
        if(NOT BLAS_LIBRARIES)
            message("-- No external BLAS library found")
        endif()
    endif()
    if(NOT FORCE_OWN_LAPACK)
        if(NOT LAPACK_LIBRARIES)
            find_package(LAPACK)
            if(LAPACK_LIBRARIES)
                set(LIBS
                    ${LAPACK_LIBRARIES}
                    ${LIBS}
                    )
            endif()
        endif()
        if(NOT LAPACK_LIBRARIES)
            message("-- No external LAPACK library found")
        endif()
    endif()
endif()

if(NOT BLAS_LIBRARIES OR FORCE_OWN_BLAS)
    message("-- Using own BLAS implementation (slow)")
    add_definitions(-DUSE_OWN_BLAS)
endif()
if(NOT LAPACK_LIBRARIES OR FORCE_OWN_LAPACK)
    message("-- Using own LAPACK implementation (slow)")
    add_definitions(-DUSE_OWN_DLAMCH)
    add_definitions(-DUSE_OWN_DLAPACK)
    add_definitions(-DUSE_OWN_ZLAPACK)
endif()
