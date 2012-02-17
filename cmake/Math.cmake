set(MATH_LANG "Fortran")

set(EXPLICIT_LIBS
    "${EXPLICIT_LIBS}"
    CACHE STRING
    "User set math libraries"
    FORCE
    )

set(BLAS_FOUND FALSE)
option(FORCE_OWN_BLAS "Force use of own BLAS" OFF)
mark_as_advanced(FORCE_OWN_BLAS)

if(ENABLE_BLAS AND NOT EXPLICIT_LIBS)
    if(FORCE_OWN_BLAS)
        add_definitions(-DUSE_OWN_BLAS)
    else()
        find_package(BLAS)
        if(BLAS_FOUND)
            set(LIBS
                ${LIBS}
                ${BLAS_LIBRARIES}
                )
        else()
            message("-- No external BLAS library found")
            message("-- Using own BLAS implementation (slow)")
        endif()
    endif()
endif()

set(LAPACK_FOUND FALSE)
option(FORCE_OWN_LAPACK "Force use of own LAPACK" OFF)
mark_as_advanced(FORCE_OWN_LAPACK)

if(ENABLE_LAPACK AND NOT EXPLICIT_LIBS)
    if(FORCE_OWN_LAPACK)
        add_definitions(-DUSE_OWN_LAPACK)
    else()
        find_package(LAPACK)
        if(LAPACK_FOUND)
            set(LIBS
                ${LIBS}
                ${LAPACK_LIBRARIES}
                )
        else()
            message("-- No external LAPACK library found")
            message("-- Using own LAPACK implementation (slow)")
        endif()
    endif()
endif()

if(EXPLICIT_LIBS)
    set(LIBS
        ${LIBS}
        ${EXPLICIT_LIBS}
        )
    message("-- User set math libraries: ${EXPLICIT_LIBS}")
endif()
