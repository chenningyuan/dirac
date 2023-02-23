add_library(input)

target_sources(input
  PRIVATE
    input_reader.F90
    input_reader_sections.F90
    parse_input.F90
  )

target_link_libraries(input
  PRIVATE
    reladc
    gp
    dirac
    main
    visual
    relccsd
    grid
    qcorr
  )

if(ENABLE_EXATENSOR)
  target_link_libraries(input
    PRIVATE
      exacorr
    )
endif()

target_include_directories(input
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
  )
