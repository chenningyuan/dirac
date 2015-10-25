option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)

if(ENABLE_UNIT_TESTS)
    message(STATUS "Unit control tests ENABLED")

    execute_process(COMMAND mkdir -p ${PROJECT_BINARY_DIR}/unit)
    execute_process(
        COMMAND python ${PROJECT_SOURCE_DIR}/cmake/custom/unit/gather.py
                       ${PROJECT_SOURCE_DIR}/external/unit-tests
                       ${PROJECT_BINARY_DIR}/unit/list_of_unit_tests.cmake
        )
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR}/unit)

    # set unit test libs
    set(UNIT_TEST_LIBS
        dirac
        ${EXTERNAL_LIBS}
        )

    macro(add_unit_test _path _name)
        add_custom_target(
            ${_name}_files_copied
            COMMAND cp ${_path}/* ${PROJECT_BINARY_DIR}/unit-tests/${_name}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            )
        add_executable(
            ${_name}.x
            ${_path}/test.F90
            )
        set_target_properties(
            ${_name}.x
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY
            "unit-tests/${_name}"
            )
        target_link_libraries(
            ${_name}.x
            ${UNIT_TEST_LIBS}
            )
        set_property(
            TARGET ${_name}.x
            PROPERTY
            LINKER_LANGUAGE Fortran
            )
        add_dependencies(
            ${_name}.x
            ${_name}_files_copied
            )
        add_test(
            ${_name}
            python ${PROJECT_SOURCE_DIR}/cmake/custom/unit-test-wrap.py --build-dir=${PROJECT_BINARY_DIR} --name=${_name}
            )
    endmacro()

    include(list_of_unit_tests)
else()
    message(STATUS "Unit control tests DISABLED")
endif()
