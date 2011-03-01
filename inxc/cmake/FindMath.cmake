
# radovan: this routine is very primitive
#          feel free to refactor/generalize

macro(Check_Fortran_Libraries LIBRARIES _prefix _list)
    set(${LIBRARIES})
    set(_libraries_work TRUE)
    set(_paths_ok TRUE)
    foreach(_library ${_list})
        if(_libraries_work)

            set(TRY_PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH)

            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
            find_library(${_prefix}_${_library}_LIBRARY
                NAMES ${_library}
                PATHS ${TRY_PATHS}
                )
            set(_libraries_work ${${_prefix}_${_library}_LIBRARY})

            find_path(DIR
                NAMES lib${_library}.a
                PATHS ${TRY_PATHS}
                )
            if(${DIR} MATCHES "NOTFOUND")
                set(_paths_ok FALSE)
            endif()
            if("${${LIBRARIES}}" STREQUAL "")
                set(LEADING "")
            else()
                set(LEADING " ")
            endif()
            set(${LIBRARIES} "${${LIBRARIES}}${LEADING}-L${DIR} -l${_library}")

        endif()
    endforeach(_library ${_list})
    if(NOT _libraries_work)
        set(${LIBRARIES} FALSE)
    endif()
    if(NOT _paths_ok)
        set(${LIBRARIES} FALSE)
    endif()
endmacro()

set(MATH_FOUND FALSE)

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

# try MKL
# if(NOT MATH_FOUND)
#     if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
#         check_fortran_libraries(
#             MKL_LIBRARIES
#             MKL
#             "mkl_core;mkl_intel_lp64;mkl_sequential;mkl_lapack;guide;pthread"
#             )
#     else()
#         check_fortran_libraries(
#             MKL_LIBRARIES
#             MKL
#             "mkl_core;mkl_intel;mkl_sequential;mkl_lapack;guide;pthread"
#             )
#     endif()
#     if(MKL_LIBRARIES)
#         set(MATH_LIBS
#             ${MKL_LIBRARIES}
#             )
# #       message("found MKL libraries: ${MKL_LIBRARIES}")
#         set(MATH_FOUND TRUE)
#     endif()
# endif()
# mkl code above is too stupid for the moment deactivate
# i will come back to it asap

# try default
if(NOT MATH_FOUND)
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND LAPACK_FOUND)
        set(MATH_LIBS
            blas
            lapack
            )
        message("-- Set system default math libraries: ${MATH_LIBS} (${LAPACK_lapack_LIBRARY} ${BLAS_f77blas_LIBRARY})")
        set(MATH_FOUND TRUE)
    endif()
endif()

# nothing found
if(NOT MATH_FOUND)
    message("-- No math libraries found - defaulting to DIRAC own math libraries")
endif()
