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

# set number of CPU for test
if(NOT "${RUNTEST_NR_PROC}" STREQUAL "")
    set(RUNTEST_FLAG_NR_PROC "--mpi=${RUNTEST_NR_PROC}")
endif()
# set timeout (this is not DART_TESTING_TIMEOUT)
if(NOT "${RUNTEST_TIMEOUT}" STREQUAL "")
    set(RUNTEST_FLAG_TIMEOUT "--timeout=${RUNTEST_TIMEOUT}")
endif()
# choose full path binary - override diracrc setting
if(NOT "${RUNTEST_BINARY}" STREQUAL "")
    set(RUNTEST_FLAG_BINARY "--binary=${RUNTEST_BINARY}")
endif()

include(Tests)
include(CTest)
enable_testing()
