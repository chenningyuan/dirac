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

include(Tests)
include(CTest)
enable_testing()
