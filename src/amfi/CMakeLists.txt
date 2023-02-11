add_library(amfi)

target_sources(amfi
  PRIVATE
    amfi.F
    amfi1.F
    amfi2.F
    amfi3.F
    relscf.F
  )

target_include_directories(amfi
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )
