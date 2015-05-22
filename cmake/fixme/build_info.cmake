#.rst:
#
# Generates source code that echoes configuration, version, and build
# information to the program output.

find_package(PythonInterp REQUIRED)

# generated cmake files will be written to this path
# only build info is generated
# everything else is static for the user
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated_by_cmake)
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_BINARY_DIR}/generated_by_cmake
    )

#-------------------------------------------------------------------------------

# get configuration time in UTC
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import datetime; print(datetime.datetime.utcnow())"
    OUTPUT_VARIABLE _configuration_time
    )
string(STRIP ${_configuration_time} _configuration_time) # delete newline

execute_process(
    COMMAND whoami
    TIMEOUT 1
    OUTPUT_VARIABLE _user_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND hostname
    TIMEOUT 1
    OUTPUT_VARIABLE _host_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

get_directory_property(_list_of_definitions DIRECTORY ${PROJECT_SOURCE_DIR} COMPILE_DEFINITIONS)

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

configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/lib/config_info.in.py
    ${PROJECT_BINARY_DIR}/config_info.py
    )

unset(_configuration_time)
unset(_user_name)
unset(_host_name)
unset(_list_of_definitions)
unset(_python_version)

exec_program(
    ${PYTHON_EXECUTABLE}
    ${PROJECT_BINARY_DIR}
    ARGS config_info.py Fortran > ${PROJECT_BINARY_DIR}/config_info.F90
    OUTPUT_VARIABLE _discard # we don't care about the output
    )

exec_program(
    ${PYTHON_EXECUTABLE}
    ${PROJECT_BINARY_DIR}
    ARGS config_info.py CMake > ${PROJECT_BINARY_DIR}/generated_by_cmake/config_info_generated.cmake
    OUTPUT_VARIABLE _discard # we don't care about the output
    )

include(config_info_generated)
