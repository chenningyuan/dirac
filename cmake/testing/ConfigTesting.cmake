# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest's own timeout in secs
set(DART_TESTING_TIMEOUT
    "600"
    CACHE STRING
    "Set default timeout in seconds for every single test"
    )

# set timeout (this is not DART_TESTING_TIMEOUT)
if(NOT "${RUNTEST_TIMEOUT}" STREQUAL "")
    set(RUNTEST_FLAG_TIMEOUT "--timeout=${RUNTEST_TIMEOUT}")
endif()

include(Tests)
include(CTest)
enable_testing()
