#
# Script for regulating the CTest part of the build
#
set(CTEST_PROJECT_NAME       "DIRAC")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CEST")
set(CTEST_DROP_METHOD        "http")
set(CTEST_DROP_SITE          "repo.ctcc.no")
set(CTEST_DROP_LOCATION      "/CDash/submit.php?project=DIRAC")
set(CTEST_DROP_SITE_CDASH    TRUE)
