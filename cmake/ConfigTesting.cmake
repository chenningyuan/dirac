# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${cmake_build_type_tolower}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest's own timeout in secs
set(DART_TESTING_TIMEOUT
    "1500"
    CACHE STRING
    "Set default timeout in seconds for every single test"
    )

macro(dirac_test _name _labels)
    add_test(${_name} python ${CMAKE_SOURCE_DIR}/test/${_name}/test --binary-dir=${PROJECT_BINARY_DIR} --work-dir=${PROJECT_BINARY_DIR}/test/${_name})
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

include(Tests)
include(CTest)
enable_testing()
