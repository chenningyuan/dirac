#.rst:
#
# Creates the build_info.h file in the build directory.
# This file can be included into sources to print out
# build information variables to the executable program output.
#

#miro: major part this code is repeated in autocmake's git_info_sub.cmake
# TODO: refeactor the code at autocmake's side
include(FindGit)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 --abbrev-commit HEAD
        OUTPUT_VARIABLE _git_last_commit_hash
        ERROR_QUIET
        )
    if(_git_last_commit_hash)
        string(STRIP ${_git_last_commit_hash} _git_last_commit_hash)
    endif()
    execute_process(
        COMMAND ${GIT_EXECUTABLE} log --max-count=1 HEAD
        OUTPUT_VARIABLE _git_last_commit
        ERROR_QUIET
        )
    if(_git_last_commit)
        string(REGEX MATCH "Author:[ ]*(.+)<.*@.*>" temp "${_git_last_commit}")
        set(_git_last_commit_author ${CMAKE_MATCH_1})
        string(STRIP ${_git_last_commit_author} _git_last_commit_author)
        string(REGEX MATCH "Date:[ ]*(.+(\\+|-)[0-9][0-9][0-9][0-9])" temp "${_git_last_commit}")
        set(_git_last_commit_date ${CMAKE_MATCH_1})
    endif()
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE _git_branch
        ERROR_QUIET
        )
    string(STRIP ${_git_branch} _git_branch)
# printout  
    message(STATUS "-------------------------------------------------------------")
    message(STATUS "Git branch               : ${_git_branch}")
    message(STATUS "Git last commit hash     : ${_git_last_commit_hash}")
    message(STATUS "Git last commit date     : ${_git_last_commit_date}")
  # some authors have characters with non-standard encoding which can make setup choke
  # message(STATUS "Git last commit author   : ${_git_last_commit_author}")
    message(STATUS "-------------------------------------------------------------")
endif()

set(_configuration_time     "unknown")
# get configuration time in UTC
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import datetime; print(datetime.datetime.utcnow())"
    OUTPUT_VARIABLE _configuration_time
    )
string(STRIP ${_configuration_time} _configuration_time) # delete newline

set(_user_name              "unknown")
execute_process(
    COMMAND whoami
    TIMEOUT 1
    OUTPUT_VARIABLE _user_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
message(STATUS "User name                : ${_user_name}")

set(_host_name              "unknown")
execute_process(
    COMMAND hostname
    TIMEOUT 1
    OUTPUT_VARIABLE _host_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
message(STATUS "Host name                : ${_host_name}")

set(_system_name            "unknown")
if (CMAKE_SYSTEM_NAME)
    set(_system_name ${CMAKE_SYSTEM_NAME})
endif()

set(_system                 "unknown")
if (CMAKE_SYSTEM)
    set(_system ${CMAKE_SYSTEM})
    message(STATUS "Operating system         : ${CMAKE_SYSTEM}")
endif()

# set python version variable when it is not defined automatically
# this can happen when cmake version is equal or less than 2.8.5
if(NOT _python_version)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -V
        OUTPUT_VARIABLE _python_version
        ERROR_VARIABLE _python_version
        )
    string(REGEX MATCH "Python ([0-9].[0-9].[0-9])" temp "${_python_version}")
    set(_python_version ${CMAKE_MATCH_1})
endif()
message(STATUS "Python version           : ${_python_version}")

set(_cmake_version          "unknown")
if (CMAKE_VERSION)
    set(_cmake_version     ${CMAKE_VERSION})
    message(STATUS "CMake version            : ${CMAKE_VERSION}")
endif()

set(_cmake_generator       "unknown")
if (CMAKE_GENERATOR)
    set(_cmake_generator     ${CMAKE_GENERATOR})
    message(STATUS "CMake generator          : ${CMAKE_GENERATOR}")
endif()

set(_cmake_build_type         "unknown")
if (CMAKE_BUILD_TYPE)
    set(_cmake_build_type ${CMAKE_BUILD_TYPE})
    message(STATUS "CMake build type         : ${CMAKE_BUILD_TYPE}")
endif()

set(_enable_mpi ${ENABLE_MPI})
message(STATUS "MPI parallelization      : ${ENABLE_MPI}")

set(_enable_64bit_integers ${ENABLE_64BIT_INTEGERS})
message(STATUS "64-bit integers          : ${ENABLE_64BIT_INTEGERS}")

set(_fortran_compiler_id       "unknown")
if (CMAKE_Fortran_COMPILER_ID)
    set(_fortran_compiler_id ${CMAKE_Fortran_COMPILER_ID})
    #message("-- Fortran compiler ID      : ${CMAKE_Fortran_COMPILER_ID}")
endif()

set(_fortran_compiler       "unknown")
if (CMAKE_Fortran_COMPILER)
    set(_fortran_compiler ${CMAKE_Fortran_COMPILER})
    message(STATUS "Fortran compiler         : ${CMAKE_Fortran_COMPILER}")
endif()

set(_python_fortran_version "unknown")
if (PYTHON_Fortran_VERSION)
    set(_python_fortran_version ${PYTHON_Fortran_VERSION})
    message(STATUS "Fortran compiler version : ${CMAKE_Fortran_COMPILER_ID}  ${PYTHON_Fortran_VERSION}")
endif()

set(_fortran_compiler_flags "unknown")
if (CMAKE_Fortran_FLAGS)
    set(_fortran_compiler_flags ${CMAKE_Fortran_FLAGS})
    message(STATUS "Fortran compiler flags   : ${CMAKE_Fortran_FLAGS}")
endif()

set(_c_compiler_id       "unknown")
if (CMAKE_C_COMPILER_ID)
    set(_c_compiler_id ${CMAKE_C_COMPILER_ID})
endif()

set(_c_compiler       "unknown")
if (CMAKE_C_COMPILER)
    set(_c_compiler ${CMAKE_C_COMPILER})
    message(STATUS "C compiler               : ${CMAKE_C_COMPILER}")
endif()

set(_python_c_version "unknown")
if (PYTHON_C_VERSION)
    set(_python_c_version ${PYTHON_C_VERSION})
    message(STATUS "C compiler version       : ${CMAKE_C_COMPILER_ID} ${PYTHON_C_VERSION}")
endif()

set(_c_compiler_flags "unknown")
if (CMAKE_C_FLAGS)
    set(_c_compiler_flags ${CMAKE_C_FLAGS})
    message(STATUS "C compiler flags         : ${CMAKE_C_FLAGS}")
endif()

set(_cxx_compiler_id       "unknown")
if (CMAKE_CXX_COMPILER_ID)
    set(_cxx_compiler_id ${CMAKE_CXX_COMPILER_ID})
endif()

set(_cxx_compiler           "unknown")
if (CMAKE_CXX_COMPILER)
    set(_cxx_compiler ${CMAKE_CXX_COMPILER})
    message(STATUS "C++ compiler             : ${CMAKE_CXX_COMPILER}")
endif()

set(_python_cxx_version "unknown")
if (PYTHON_CXX_VERSION)
    set(_python_cxx_version ${PYTHON_CXX_VERSION})
    message(STATUS "C++ compiler version     : ${CMAKE_CXX_COMPILER_ID}  ${PYTHON_CXX_VERSION}")
endif()

set(_cxx_compiler_flags     "unknown")
if (CMAKE_CXX_FLAGS)
    set(_cxx_compiler_flags ${CMAKE_CXX_FLAGS})
    message(STATUS "C++ compiler flags       : ${CMAKE_CXX_FLAGS}")
endif()

set(_static_linking         "unknown")
if (CMAKE_SYSTEM_NAME)
    set(_static_linking ${ENABLE_STATIC_LINKING})
    message(STATUS "Static linking           : ${ENABLE_STATIC_LINKING}")
endif()

set(_enable_builtin_blas ${ENABLE_BUILTIN_BLAS})
if (ENABLE_BUILTIN_BLAS)
    message(STATUS "Builtin BLAS (slow)      : ${ENABLE_BUILTIN_BLAS}")
endif()

set(_enable_builtin_lapack ${ENABLE_BUILTIN_LAPACK})
if (ENABLE_BUILTIN_LAPACK)
    message(STATUS "Builtin LAPACK (slow)    : ${ENABLE_BUILTIN_LAPACK}")
endif()

set(_math_libs         "unknown")
if (MATH_LIBS)
    set(_math_libs ${MATH_LIBS})
    message(STATUS "Mathematical libraries   : ${MATH_LIBS}")
endif()

set(_explicit_libs       "unknown")
if (EXPLICIT_LIBS)
    set(_explicit_libs ${EXPLICIT_LIBS})
    message(STATUS " Explicit libraries      : ${EXPLICIT_LIBS}")
endif()

set(_mkl_flag         "unknown")
if (MKL_FLAG)
    set(_mkl_flag ${MKL_FLAG})
    message(STATUS "Intel MKL math libraries : ${MKL_FLAG}")
endif()

get_directory_property(_list_of_definitions DIRECTORY ${PROJECT_SOURCE_DIR} COMPILE_DEFINITIONS)
message(STATUS "Definitions              : ${_list_of_definitions}")

message(STATUS "Configuration time       : ${_configuration_time}")

# generate the build_info.h include file with defined set of variables
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/build_info.h.in
    ${PROJECT_BINARY_DIR}/build_info.h
    @ONLY
)

add_custom_target(
    build_info
    ALL DEPENDS ${PROJECT_BINARY_DIR}/build_info.h
    )
