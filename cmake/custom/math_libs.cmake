#.rst:
#
# Detects and links to BLAS and LAPACK libraries using CMake's own
# FindBLAS and FindLAPACK modules. Replaces the deprecated Autocmake
# math_libs.cmake module and keeps its interface towards the rest of DIRAC.
#
# On Cray systems the compiler wrappers (ftn, cc, CC) link Cray LibSci
# implicitly; FindBLAS/FindLAPACK detect this ("implicitly linked") and
# MATH_LIBS then stays empty.
#
# Variables used::
#
#   ENABLE_BLAS             (auto or off)
#   ENABLE_LAPACK           (auto or off)
#   ENABLE_BUILTIN_BLAS
#   ENABLE_BUILTIN_LAPACK
#   ENABLE_64BIT_INTEGERS
#   ENABLE_STATIC_LINKING
#   BLA_VENDOR              (see CMake's FindBLAS documentation; All searches all)
#   MKL_FLAG
#   MATH_LIBS               (if set by the user, detection is skipped;
#                            BLAS_TYPE may then also be set by the user)
#
# Variables defined::
#
#   MATH_LIBS
#   BLAS_FOUND, LAPACK_FOUND
#   BLAS_TYPE, LAPACK_TYPE  (MKL, OPENBLAS, LIBSCI, ATLAS, ESSL, ...)
#
# Preprocessor definitions::
#
#   HAVE_<TYPE>_BLAS, HAVE_<TYPE>_LAPACK  if a library is found
#   USE_BUILTIN_BLAS, USE_BUILTIN_LAPACK  if none is found
#
# Environment variables used::
#
#   MATH_ROOT
#   BLAS_ROOT, LAPACK_ROOT  (used directly by find_package)
#   MKL_ROOT, MKLROOT
#   CRAY_LIBSCI_PREFIX_DIR
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--blas=<BLAS> Detect and link BLAS library (auto or off) [default: auto]."
#     - "--lapack=<LAPACK> Detect and link LAPACK library (auto or off) [default: auto]."
#     - "--math-vendor=<BLA_VENDOR> BLAS/LAPACK vendor for CMake's FindBLAS, e.g. Intel10_64lp, OpenBLAS, FlexiBLAS, AOCL, Arm, Generic (All searches all) [default: All]."
#     - "--qmkl=<MKL> Pass MKL flag to the Intel compiler and linker and skip BLAS/LAPACK detection (sequential, parallel, cluster, or off) [default: off]."
#   define:
#     - "'-DENABLE_BLAS={0}'.format(arguments['--blas'])"
#     - "'-DENABLE_LAPACK={0}'.format(arguments['--lapack'])"
#     - "'-DBLA_VENDOR={0}'.format(arguments['--math-vendor'])"
#     - "'-DMKL_FLAG={0}'.format(arguments['--qmkl'])"

# guess the library type (as used elsewhere in DIRAC) from the library names
function(_dirac_math_type _out)
    set(_type)
    foreach(_lib ${ARGN})
        get_filename_component(_name "${_lib}" NAME)
        string(TOLOWER "${_name}" _name)
        if(_name MATCHES "mkl")
            set(_type MKL)
        elseif(_name MATCHES "flexiblas")
            set(_type FLEXIBLAS)
        elseif(_name MATCHES "openblas")
            set(_type OPENBLAS)
        elseif(_name MATCHES "sci_")
            set(_type LIBSCI)
        elseif(_name MATCHES "essl")
            set(_type ESSL)
        elseif(_name MATCHES "acml")
            set(_type ACML)
        elseif(_name MATCHES "armpl")
            set(_type ARMPL)
        elseif(_name MATCHES "atlas")
            set(_type ATLAS)
        elseif(_name MATCHES "blis|flame")
            set(_type BLIS)
        elseif(_name MATCHES "ssl2|fjlapack")
            set(_type FUJITSU)
        elseif(_name MATCHES "accelerate")
            set(_type APPLE)
        endif()
        if(_type)
            break()
        endif()
    endforeach()
    if(NOT _type)
        if("${ARGN}" STREQUAL "")
            # linked implicitly by the compiler (wrapper)
            if(DEFINED ENV{CRAY_LIBSCI_PREFIX_DIR})
                set(_type LIBSCI)
            else()
                set(_type IMPLICIT)
            endif()
        elseif(BLA_VENDOR AND NOT BLA_VENDOR STREQUAL "All")
            string(TOUPPER "${BLA_VENDOR}" _type)
        else()
            set(_type SYSTEM_NATIVE)
        endif()
    endif()
    set(${_out} ${_type} PARENT_SCOPE)
endfunction()

if(NOT BLA_VENDOR)
    set(BLA_VENDOR All)
endif()

# avoid warnings about unused variables when detection is skipped
set(_unused ${ENABLE_BLAS} ${ENABLE_LAPACK})
unset(_unused)

if(ENABLE_STATIC_LINKING)
    set(BLA_STATIC ON)
endif()

# FindBLAS/FindLAPACK check the integer size (CMake 3.22 or later)
if(NOT DEFINED BLA_SIZEOF_INTEGER)
    if(ENABLE_64BIT_INTEGERS)
        set(BLA_SIZEOF_INTEGER 8)
    else()
        set(BLA_SIZEOF_INTEGER 4)
    endif()
endif()

# BLAS_ROOT and LAPACK_ROOT (variables or environment) are honored by find_package itself
if(DEFINED ENV{MATH_ROOT})
    foreach(_service BLAS LAPACK)
        if(NOT ${_service}_ROOT AND NOT DEFINED ENV{${_service}_ROOT})
            set(${_service}_ROOT $ENV{MATH_ROOT})
            message(STATUS "${_service} will be searched for based on MATH_ROOT=${${_service}_ROOT}")
        endif()
    endforeach()
endif()
if(DEFINED ENV{MKL_ROOT} AND NOT DEFINED ENV{MKLROOT})
    set(ENV{MKLROOT} $ENV{MKL_ROOT})
endif()

if(MKL_FLAG)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} -qmkl=${MKL_FLAG})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qmkl=${MKL_FLAG}")
    message(STATUS "User set explicit MKL flag which is passed to the compiler and linker: -qmkl=${MKL_FLAG}")
    message(STATUS "This disables math detection and builtin math libraries")
    foreach(_service BLAS LAPACK)
        set(${_service}_FOUND TRUE)
        set(${_service}_TYPE MKL CACHE STRING "${_service} type" FORCE)
    endforeach()
# MATH_LIBS set by the user (and not by an earlier run of this module): no detection
elseif(NOT "${MATH_LIBS}" STREQUAL "" AND NOT "${MATH_LIBS}" STREQUAL "${_DIRAC_DETECTED_MATH_LIBS}")
    message(STATUS "Using user-defined MATH_LIBS, BLAS/LAPACK detection skipped")
    foreach(_service BLAS LAPACK)
        set(${_service}_FOUND TRUE)
    endforeach()
else()
    message(STATUS "Searching for BLAS and LAPACK (BLA_VENDOR=${BLA_VENDOR})")
    foreach(_service BLAS LAPACK)
        unset(${_service}_TYPE CACHE)
        if(ENABLE_${_service} STREQUAL "auto" AND NOT ENABLE_BUILTIN_${_service})
            # FindLAPACK also searches for BLAS
            find_package(${_service})
            if(${_service}_FOUND)
                _dirac_math_type(_type ${${_service}_LIBRARIES})
                set(${_service}_TYPE ${_type} CACHE STRING "${_service} type" FORCE)
                mark_as_advanced(${_service}_TYPE)
                message(STATUS "${_service} type: ${${_service}_TYPE}")
            else()
                add_definitions(-DUSE_BUILTIN_${_service})
                set(USE_BUILTIN_${_service} TRUE)
                message(STATUS "No ${_service} library found, using builtin ${_service} (slow)")
            endif()
        endif()
    endforeach()

    # LAPACK_LIBRARIES already contains BLAS_LIBRARIES
    if(LAPACK_FOUND AND NOT ENABLE_BUILTIN_LAPACK)
        set(_math_libs ${LAPACK_LIBRARIES})
        separate_arguments(_flags UNIX_COMMAND "${LAPACK_LINKER_FLAGS}")
    elseif(BLAS_FOUND AND NOT ENABLE_BUILTIN_BLAS)
        set(_math_libs ${BLAS_LIBRARIES})
        separate_arguments(_flags UNIX_COMMAND "${BLAS_LINKER_FLAGS}")
    endif()
    set(_math_libs ${_math_libs} ${_flags})

    # further adaptation for the static linking
    if(ENABLE_STATIC_LINKING AND _math_libs)
        if(BLAS_TYPE MATCHES "ATLAS|SYSTEM_NATIVE|OPENBLAS" OR LAPACK_TYPE MATCHES "ATLAS|SYSTEM_NATIVE|OPENBLAS")
            set(_math_libs ${_math_libs} -Wl,--whole-archive -lpthread -ldl -Wl,--no-whole-archive -lm)
        endif()
        if(BLAS_TYPE MATCHES MKL OR LAPACK_TYPE MATCHES MKL)
            set(_math_libs ${_math_libs} -ldl)
        endif()
    endif()

    set(MATH_LIBS "${_math_libs}" CACHE STRING "Math libraries" FORCE)
    set(_DIRAC_DETECTED_MATH_LIBS "${_math_libs}" CACHE INTERNAL "MATH_LIBS as detected by math_libs.cmake")
    unset(_math_libs)
    unset(_flags)
    unset(_type)
endif()

foreach(_service BLAS LAPACK)
    if(${_service}_TYPE)
        add_definitions(-DHAVE_${${_service}_TYPE}_${_service})
        message(STATUS "Setting -DHAVE_${${_service}_TYPE}_${_service}")
        set(HAVE_${_service} ON CACHE INTERNAL "Defined if ${_service} is available")
    endif()
    set(${_service}_FOUND ${${_service}_FOUND} CACHE BOOL "${_service} found" FORCE)
endforeach()
unset(_service)
