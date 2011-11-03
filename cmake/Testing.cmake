# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest's own timeout
set(DART_TESTING_TIMEOUT
    "600"
    CACHE STRING
    "Set default timeout in seconds for every single test"
    )

if(NOT "${RUNTEST_NR_PROC}" STREQUAL "")
    set(RUNTEST_FLAG_NR_PROC "--mpi=${RUNTEST_NR_PROC}")
    message("-- runtest extra flag active, RUNTEST_FLAG_NR_PROC=\"${RUNTEST_FLAG_NR_PROC}\" ")
endif()
if(NOT "${RUNTEST_TIMEOUT}" STREQUAL "")
    set(RUNTEST_FLAG_TIMEOUT "--timeout=${RUNTEST_TIMEOUT}")
    message("-- runtest extra flag active, RUNTEST_FLAG_TIMEOUT=\"${RUNTEST_FLAG_TIMEOUT}\" ")
endif()

include(Tests)
include(CTest)
enable_testing()
