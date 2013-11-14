# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${cmake_build_type_tolower}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest's own timeout in secs
set(DART_TESTING_TIMEOUT
    "600"
    CACHE STRING
    "Set default timeout in seconds for every single test"
    )

include(Tests)
include(CTest)
enable_testing()
