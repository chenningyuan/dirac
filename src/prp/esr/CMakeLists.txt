add_library(esr)

target_sources(esr
  PRIVATE
    pamesr.F
    pamesrci.F
)

target_include_directories(esr
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

dirac_test(esr "medium;gascip;ci" "")
