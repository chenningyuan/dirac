add_library(moltra)

target_sources(moltra
  PRIVATE
    tradr1.F
    tradr2.F
    tradr4.F
    tradr6.F
    trainp.F
    traone.F
    traout.F
    trapam.F
    traprp.F
    dcbtra_h.F90
    tradr3.F90
    tradr5.F90
  )

target_include_directories(moltra
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )

target_link_libraries(moltra
  PRIVATE
    gp
    relccsd
  )
